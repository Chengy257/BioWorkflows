###############################################
# Shared definitions: paths, sample parsing, FASTQ discovery, and helpers.
# This file is included first by workflow/Snakefile.
###############################################
import os
import re

SCRIPTS = os.path.join(WORKFLOW_DIR, "scripts")


def runtime_env(name, default=""):
    """Return a runtime value exported by run.sh/software.yaml."""
    return os.environ.get(name, default)


def tool(name, default):
    """Resolve a named executable from the unified runtime context."""
    key = f"RNASEQ_TOOL_{name.upper().replace('-', '_')}"
    return runtime_env(key, default)


RSCRIPT = runtime_env("RNASEQ_RSCRIPT", tool("rscript", "Rscript"))
PYTHON = runtime_env("RNASEQ_PYTHON", tool("python", "python3"))

# Root directory for workflow outputs, relative to the project directory.
RD = config.get("results_dir", "results").rstrip(os.sep) + os.sep


def R(path=""):
    """Return a path under the configured results directory."""
    return f"{RD}{path}"


def res(key, default=None):
    """Read a resource value after species presets have been merged."""
    return config.get(key, default)


def lres(key):
    """Read an lncRNA-specific setting, falling back to the top-level config."""
    value = (config.get("lncrna") or {}).get(key)
    if value in (None, "") and key in config:
        value = config[key]
    if value in (None, ""):
        raise ValueError(
            f"Missing required config value 'lncrna.{key}'. See the lncRNA section in config.yaml."
        )
    return value


def get_samples():
    """Read sample ids from SampleListFile; the header must contain id and group."""
    path = config["SampleListFile"]
    if not os.path.isfile(path):
        raise ValueError(
            f"SampleListFile {path!r} does not exist (configure it relative to the "
            "working directory, or as an absolute path; see config/samples.csv)"
        )
    ids = []
    with open(path, "r") as samples_list:
        next(samples_list)
        for line in samples_list:
            line = line.strip().split(",")
            ids.append(line[0])
    return ids


def validate_config(cfg):
    """Parse-time hard validation (mirrors chip_cuttag_atac_faire): required
    keys, types and value ranges aggregated into a single report. Reference
    existence only warns (--lint/dry-run often run where references are
    absent)."""
    errors, warnings = [], []
    for key in ("SampleListFile", "control_group", "results_dir", "FoldChange", "padj"):
        if key not in cfg:
            errors.append(f"missing required config key: {key}")
    if not isinstance(cfg.get("results_dir", "results"), str):
        errors.append(f"results_dir must be a string, got {cfg.get('results_dir')!r}")
    for key in ("FoldChange", "padj"):
        if key in cfg:
            try:
                float(cfg[key])
            except (TypeError, ValueError):
                errors.append(f"{key} must be numeric, got {cfg[key]!r}")
    try:
        if not 0 < float(cfg.get("padj", 0.05)) <= 1:
            errors.append(f"padj must be in (0, 1], got {cfg.get('padj')!r}")
    except (TypeError, ValueError):
        pass  # already reported by the numeric check above
    for key in ("pca_ntop", "threads"):
        if key in cfg:
            v = cfg[key]
            if isinstance(v, bool) or not isinstance(v, int) or v < 1:
                errors.append(f"{key} must be an integer >= 1, got {v!r}")
    if cfg.get("batch_correction") not in ("T", "F"):
        errors.append(f"batch_correction must be 'T' or 'F', got {cfg.get('batch_correction')!r}")
    trim = cfg.get("trim")
    if trim is not None:
        if not isinstance(trim, dict):
            errors.append(f"trim must be a mapping of {{quality, stringency, error_rate, extra}}, got {trim!r}")
        else:
            for key, lo in (("quality", 0), ("stringency", 1)):
                if key in trim:
                    v = trim[key]
                    if isinstance(v, bool) or not isinstance(v, int) or v < lo:
                        errors.append(f"trim.{key} must be an integer >= {lo}, got {v!r}")
            if "error_rate" in trim:
                try:
                    er = float(trim["error_rate"])
                    if not 0 < er <= 1:
                        errors.append(f"trim.error_rate must be in (0, 1], got {er!r}")
                except (TypeError, ValueError):
                    errors.append(f"trim.error_rate must be numeric, got {trim['error_rate']!r}")
            if "extra" in trim and not isinstance(trim["extra"], str):
                errors.append(f"trim.extra must be a string, got {trim['extra']!r}")
    umi = cfg.get("umi")
    if umi is not None:
        if not isinstance(umi, dict):
            errors.append(f"umi must be a mapping of {{enabled, read1_len, read2_len}}, got {umi!r}")
        else:
            if "enabled" in umi and not isinstance(umi["enabled"], bool):
                errors.append(f"umi.enabled must be true/false, got {umi['enabled']!r}")
            for key in ("read1_len", "read2_len"):
                if key in umi:
                    v = umi[key]
                    if isinstance(v, bool) or not isinstance(v, int) or v < 0:
                        errors.append(f"umi.{key} must be an integer >= 0, got {v!r}")
    resources = cfg.get("resources") or {}
    if not isinstance(resources, dict):
        errors.append(f"resources must be a mapping of rule -> {{threads, mem_mb, runtime_min}}, got {resources!r}")
    else:
        for rule, entry in resources.items():
            if not isinstance(entry, dict):
                errors.append(f"resources.{rule} must be a mapping, got {entry!r}")
                continue
            for field in entry:
                if field not in ("threads", "mem_mb", "runtime_min"):
                    errors.append(f"resources.{rule}.{field}: unknown field (supported: threads/mem_mb/runtime_min)")
    if isinstance(umi, dict) and umi.get("enabled"):
        # A malformed read1_len is already flagged by the umi type checks
        # above; the conflict probe must not crash on it before the
        # aggregated report is raised.
        read1_len = umi.get("read1_len")
        if (
            isinstance(read1_len, int)
            and not isinstance(read1_len, bool)
            and read1_len > 0
            and "--clip5pNbases" in str(cfg.get("star_extra_args") or "")
        ):
            errors.append(
                "conflicting UMI clipping: umi.read*_len and --clip5pNbases in star_extra_args "
                "would both be appended to the STAR command; keep only one of them"
            )
    if errors:
        raise ValueError(
            f"config validation failed ({len(errors)} issues):\n  " + "\n  ".join(errors)
        )
    for key in ("genome", "gtf", "bed"):
        value = cfg.get(key)
        if not isinstance(value, str) or not value:
            continue
        if "/path/to/" in value:
            warnings.append(f"{key} is still a placeholder: {value} (fill real paths in the project config)")
        elif not os.path.exists(value) and not os.path.isabs(value):
            warnings.append(f"{key} is relative and not found from the working directory: {value}")
    for w in warnings:
        print(f"[config warning] {w}")


SAMPLES = get_samples()
SAMPLE_WILDCARD = "(?:" + "|".join(re.escape(sample) for sample in SAMPLES) + ")"

validate_config(config)


def _raw_reads(sample):
    """Detect paired-end first, then single-end FASTQ files using exact names.

    Besides the project-native {id}_1/_2 names, the common sequencer-delivery
    names {id}_R1/_R2 (.fastq.gz / .fq.gz) are accepted; the _1/_2 patterns
    keep priority so existing projects resolve identically."""
    for pattern_1, pattern_2 in (
        ("1.rawdata/{0}_1.fastq.gz", "1.rawdata/{0}_2.fastq.gz"),
        ("1.rawdata/{0}_1.fq.gz", "1.rawdata/{0}_2.fq.gz"),
        ("1.rawdata/{0}_R1.fastq.gz", "1.rawdata/{0}_R2.fastq.gz"),
        ("1.rawdata/{0}_R1.fq.gz", "1.rawdata/{0}_R2.fq.gz"),
    ):
        if os.path.exists(pattern_1.format(sample)) and os.path.exists(pattern_2.format(sample)):
            return "PE", [pattern_1.format(sample), pattern_2.format(sample)]
    for pattern in ("1.rawdata/{0}.fastq.gz", "1.rawdata/{0}.fq.gz"):
        if os.path.exists(pattern.format(sample)):
            return "SE", [pattern.format(sample)]
    raise ValueError(
        "sample {0}: no raw FASTQ files were found under 1.rawdata/. "
        "Supported names are {0}_1.fastq.gz + {0}_2.fastq.gz, "
        "{0}_1.fq.gz + {0}_2.fq.gz, {0}_R1.fastq.gz + {0}_R2.fastq.gz, or "
        "{0}_R1.fq.gz + {0}_R2.fq.gz for paired-end data, or "
        "{0}.fastq.gz / {0}.fq.gz for single-end data.".format(sample)
    )


def get_fastq(wildcards):
    """Return trimmed FASTQ inputs for STAR according to the detected layout."""
    layout, _ = _raw_reads(wildcards.sample)
    sample = wildcards.sample
    if layout == "PE":
        return [
            R(f"2.cleandata/trim/{sample}_1_val_1.fq.gz"),
            R(f"2.cleandata/trim/{sample}_2_val_2.fq.gz"),
        ]
    return [R(f"2.cleandata/trim/{sample}_trimmed.fq.gz")]


def trimmed_reads(sample):
    """Return the STAR --readFilesIn argument for one sample."""
    layout, _ = _raw_reads(sample)
    if layout == "PE":
        return (
            f"{RD}2.cleandata/trim/{sample}_1_val_1.fq.gz "
            f"{RD}2.cleandata/trim/{sample}_2_val_2.fq.gz"
        )
    return f"{RD}2.cleandata/trim/{sample}_trimmed.fq.gz"


def is_paired_end(wildcards):
    return "True" if _raw_reads(wildcards.sample)[0] == "PE" else "False"


def trim_reports(sample):
    """Return normalized Trim Galore report paths for one sample."""
    layout, _ = _raw_reads(sample)
    if layout == "PE":
        return [
            R(f"2.cleandata/trim/{sample}_1_trimming_report.txt"),
            R(f"2.cleandata/trim/{sample}_2_trimming_report.txt"),
        ]
    return [R(f"2.cleandata/trim/{sample}_trimming_report.txt")]


def all_trim_reports(wildcards=None):
    """Return Trim Galore reports for all samples."""
    reports = []
    for sample in SAMPLES:
        reports += trim_reports(sample)
    return reports


def raw_se(wildcards):
    """Return the raw single-end FASTQ for one declared sample."""
    layout, files = _raw_reads(wildcards.sample)
    if layout != "SE":
        raise ValueError(
            f"sample {wildcards.sample}: detected {layout} data but the SE trim rule was selected"
        )
    return files[0]


def raw_pe(wildcards):
    """Return the two raw paired-end FASTQs for one declared sample."""
    layout, files = _raw_reads(wildcards.sample)
    if layout != "PE":
        raise ValueError(
            f"sample {wildcards.sample}: detected {layout} data but the PE trim rule was selected"
        )
    return {"fq1": files[0], "fq2": files[1]}


def multiqc_inputs(wildcards=None):
    """Return FastQC, trimming, STAR, and featureCounts inputs for MultiQC."""
    files = []
    for sample in SAMPLES:
        layout, _ = _raw_reads(sample)
        if layout == "PE":
            files += [
                R(f"2.cleandata/trim/fastqc/{sample}_1_val_1_fastqc.zip"),
                R(f"2.cleandata/trim/fastqc/{sample}_2_val_2_fastqc.zip"),
            ]
        else:
            files.append(R(f"2.cleandata/trim/fastqc/{sample}_trimmed_fastqc.zip"))
        files += trim_reports(sample)
        files.append(R(f"3.align/{sample}_Log.final.out"))
    if PIPELINE in ("upstream", "deg", "lncrna"):
        files += [R(f"4.expression/{sample}.log") for sample in SAMPLES]
    if PIPELINE == "lncrna":
        files += [R(f"4.lncrna/expression/{sample}.log") for sample in SAMPLES]
    return files


# Per-rule resource defaults. Project config can override any field under
# config["resources"][<name>]. The legacy top-level "threads" value is kept as
# a global thread cap for backward compatibility.
RESOURCE_DEFAULTS = {
    "software_versions": {"threads": 1, "mem_mb": 1024, "runtime_min": 10},
    "trim": {"threads": 4, "mem_mb": 8000, "runtime_min": 180},
    "multiqc": {"threads": 2, "mem_mb": 4000, "runtime_min": 60},
    "star_index": {"threads": 12, "mem_mb": 48000, "runtime_min": 480},
    "star_align": {"threads": 12, "mem_mb": 32000, "runtime_min": 360},
    "mapping_stat": {"threads": 1, "mem_mb": 2000, "runtime_min": 30},
    "strandedness": {"threads": 1, "mem_mb": 4000, "runtime_min": 60},
    "featurecounts": {"threads": 8, "mem_mb": 16000, "runtime_min": 240},
    "count_merge": {"threads": 1, "mem_mb": 4000, "runtime_min": 60},
    "deseq2": {"threads": 4, "mem_mb": 16000, "runtime_min": 240},
    "enrichment": {"threads": 2, "mem_mb": 12000, "runtime_min": 360},
    "deg_compare": {"threads": 4, "mem_mb": 12000, "runtime_min": 360},
    "stringtie": {"threads": 8, "mem_mb": 16000, "runtime_min": 360},
    "gtf_merge": {"threads": 4, "mem_mb": 12000, "runtime_min": 240},
    "isoform_expr": {"threads": 8, "mem_mb": 16000, "runtime_min": 360},
    "lnc_coding": {"threads": 8, "mem_mb": 16000, "runtime_min": 720},
    "pfam": {"threads": 16, "mem_mb": 24000, "runtime_min": 720},
    "nr": {"threads": 16, "mem_mb": 32000, "runtime_min": 1440},
    "lnc_final": {"threads": 2, "mem_mb": 8000, "runtime_min": 120},
    "lnc_featurecounts": {"threads": 8, "mem_mb": 16000, "runtime_min": 240},
    "lnc_count_merge": {"threads": 1, "mem_mb": 4000, "runtime_min": 60},
}


def _rule_resource(name, field):
    defaults = RESOURCE_DEFAULTS[name]
    overrides = (config.get("resources") or {}).get(name) or {}
    return int(overrides.get(field, defaults[field]))


def rthreads(name):
    """Return rule threads with backward-compatible legacy fallbacks."""
    overrides = (config.get("resources") or {}).get(name) or {}
    if "threads" in overrides:
        return max(1, int(overrides["threads"]))
    if name in ("pfam", "nr"):
        legacy_lnc_threads = (config.get("lncrna") or {}).get("threads")
        if legacy_lnc_threads not in (None, ""):
            return max(1, int(legacy_lnc_threads))
    requested = max(1, RESOURCE_DEFAULTS[name]["threads"])
    legacy_threads = config.get("threads")
    if legacy_threads not in (None, ""):
        return min(requested, max(1, int(legacy_threads)))
    return requested


def rmem(name):
    """Return memory requested by one rule in MB."""
    return max(1, _rule_resource(name, "mem_mb"))


def rruntime(name):
    """Return walltime requested by one rule in minutes."""
    return max(1, _rule_resource(name, "runtime_min"))


def rruntime_sec(name):
    """Return walltime requested by one rule in seconds (for SGE h_rt)."""
    return rruntime(name) * 60


# STAR parameters. The AS pipeline uses assembly-oriented splice filtering;
# other pipelines use the standard mapping configuration. star_extra_args is
# appended in both modes.
STAR_ARGS_STANDARD = (
    "--twopassMode Basic --genomeLoad NoSharedMemory "
    "--quantMode GeneCounts --outSAMattrIHstart 0 "
)
STAR_ARGS_ASSEMBLY = (
    "--twopassMode Basic --outFilterType BySJout --alignIntronMin 20 "
    "--alignIntronMax 5000 --alignMatesGapMax 5000 "
    "--outFilterMatchNminOverLread 0.66 --outFilterScoreMinOverLread 0.66 "
    "--winAnchorMultimapNmax 70 --seedSearchStartLmax 45 --outSAMattrIHstart 0 "
    "--outSAMstrandField intronMotif --genomeLoad NoSharedMemory "
    "--quantMode TranscriptomeSAM GeneCounts "
)
STAR_ARGS = (STAR_ARGS_ASSEMBLY if PIPELINE == "as" else STAR_ARGS_STANDARD) + config.get(
    "star_extra_args", ""
)


# Trim Galore arguments for both trim rules. Defaults reproduce the historical
# hardcoded command (-q 30 --stringency 3 -e 0.1); the optional trim: section
# overrides them (mirroring the chip workflow), with trim.extra appended last.
def trim_args():
    trim = config.get("trim") or {}
    args = (
        f"-q {trim.get('quality', 30)} "
        f"--stringency {trim.get('stringency', 3)} "
        f"-e {trim.get('error_rate', 0.1)}"
    )
    extra = (trim.get("extra") or "").strip()
    if extra:
        args += f" {extra}"
    return args


# Per-sample STAR arguments. Identical to STAR_ARGS unless the optional umi:
# section enables 5' UMI clipping; the --clip5pNbases string depends on the
# sample's detected PE/SE layout, so it must stay out of the global constant.
def star_args_for(sample):
    args = STAR_ARGS
    umi = config.get("umi") or {}
    read1_len = int(umi.get("read1_len", 0) or 0)
    if umi.get("enabled") and read1_len > 0:
        read2_len = int(umi.get("read2_len", 0) or 0)
        if _raw_reads(sample)[0] == "PE":
            # STAR 2.7.10b requires two comma-separated values for PE data.
            args += f" --clip5pNbases {read1_len},{read2_len or read1_len}"
        else:
            args += f" --clip5pNbases {read1_len}"
    return args
