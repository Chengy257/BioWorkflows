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
    ids = []
    with open(config["SampleListFile"], "r") as samples_list:
        next(samples_list)
        for line in samples_list:
            line = line.strip().split(",")
            ids.append(line[0])
    return ids


SAMPLES = get_samples()
SAMPLE_WILDCARD = "(?:" + "|".join(re.escape(sample) for sample in SAMPLES) + ")"


def _raw_reads(sample):
    """Detect paired-end first, then single-end FASTQ files using exact names."""
    for pattern_1, pattern_2 in (
        ("1.rawdata/{0}_1.fastq.gz", "1.rawdata/{0}_2.fastq.gz"),
        ("1.rawdata/{0}_1.fq.gz", "1.rawdata/{0}_2.fq.gz"),
    ):
        if os.path.exists(pattern_1.format(sample)) and os.path.exists(pattern_2.format(sample)):
            return "PE", [pattern_1.format(sample), pattern_2.format(sample)]
    for pattern in ("1.rawdata/{0}.fastq.gz", "1.rawdata/{0}.fq.gz"):
        if os.path.exists(pattern.format(sample)):
            return "SE", [pattern.format(sample)]
    raise ValueError(
        "sample {0}: no raw FASTQ files were found under 1.rawdata/. "
        "Supported names are {0}_1.fastq.gz + {0}_2.fastq.gz, "
        "{0}_1.fq.gz + {0}_2.fq.gz for paired-end data, or "
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
        files += [R(f"5.expression/lncRNA/{sample}.log") for sample in SAMPLES]
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
