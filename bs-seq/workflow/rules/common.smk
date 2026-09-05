# ---------------------------------------------------------------------
# Shared definitions: sample-table parsing, config validation, PE/SE
# raw-read detection, trim-aware alignment inputs, resource helpers, and
# target aggregation. Included first by workflow/Snakefile; every
# rules/*.smk uses the names defined here. BASE_DIR / WORKFLOW_DIR come
# from the Snakefile.
# ---------------------------------------------------------------------
import csv
import os
import re

from snakemake.exceptions import WorkflowError

# Sample names allow alphanumerics . _ - (no "__", no leading '-'); they
# become file names and bismark arguments.
_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def _resolve_sample_table(path):
    """Resolve the sample table: absolute > relative to the working
    directory > relative to the repository root."""
    p = os.path.expanduser(str(path))
    if os.path.isabs(p):
        return p
    if os.path.exists(p):
        return os.path.abspath(p)
    return os.path.join(BASE_DIR, p)


def load_sample_table(path):
    """Parse the single-column sample table (header must be exactly
    sample_id); returns the ordered list of unique sample ids."""
    samples = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames != ["sample_id"]:
            raise WorkflowError(
                f"Sample table {path} must have exactly one column with header "
                f"'sample_id' (got {reader.fieldnames}); see config/samples.csv"
            )
        for lineno, row in enumerate(reader, start=2):
            sid = (row["sample_id"] or "").strip()
            if not sid:
                raise WorkflowError(f"Sample table line {lineno}: sample_id must not be empty")
            if not _NAME_RE.match(sid) or "__" in sid:
                raise WorkflowError(
                    f"Sample table line {lineno}: sample_id={sid!r} contains illegal "
                    "characters; only alphanumerics and . _ - are allowed "
                    "(no leading '-', no '__')"
                )
            if sid in samples:
                raise WorkflowError(f"Sample table line {lineno}: duplicate sample_id {sid!r}")
            samples.append(sid)
    if not samples:
        raise WorkflowError(f"Sample table {path} has no data rows")
    return samples


SAMPLES = load_sample_table(_resolve_sample_table(config["SampleListFile"]))
SAMPLE_WILDCARD = "(?:" + "|".join(re.escape(s) for s in SAMPLES) + ")"


def validate_config(cfg):
    """Parse-time hard validation (mirrors srna-seq/seclip-seq): required
    keys, types and value ranges aggregated into a single report. Reference
    existence only warns (--lint/dry-run often run where references are
    absent)."""
    errors, warnings = [], []
    for key in ("SampleListFile", "results_dir", "threads", "trim",
                "bismark", "methylation_extractor", "genome"):
        if key not in cfg:
            errors.append(f"missing required config key: {key}")
    if not isinstance(cfg.get("results_dir", "results"), str):
        errors.append(f"results_dir must be a string, got {cfg.get('results_dir')!r}")
    v = cfg.get("threads")
    if isinstance(v, bool) or not isinstance(v, int) or v < 1:
        errors.append(f"threads must be an integer >= 1, got {v!r}")
    if isinstance(cfg.get("trim"), dict):
        t = cfg["trim"]
        if not isinstance(t.get("enabled"), bool):
            errors.append(f"trim.enabled must be true/false, got {t.get('enabled')!r}")
        for key, lo in (("quality", 0), ("min_len", 1), ("stringency", 1)):
            v = t.get(key)
            if isinstance(v, bool) or not isinstance(v, int) or v < lo:
                errors.append(f"trim.{key} must be an integer >= {lo}, got {v!r}")
        try:
            er = float(t.get("error_rate"))
            if not 0 < er <= 1:
                errors.append(f"trim.error_rate must be in (0, 1], got {t.get('error_rate')!r}")
        except (TypeError, ValueError):
            errors.append(f"trim.error_rate must be numeric, got {t.get('error_rate')!r}")
        adapter = t.get("adapter")
        if not isinstance(adapter, str) or not re.fullmatch(r"[ACGTN]+", adapter or ""):
            errors.append(f"trim.adapter must be a non-empty ACGTN string, got {adapter!r}")
        if not isinstance(t.get("extra", ""), str):
            errors.append("trim.extra must be a string")
    else:
        errors.append("trim must be a mapping with sub-keys")
    bismark = cfg.get("bismark")
    if not isinstance(bismark, dict) or not isinstance(bismark.get("align_extra", ""), str):
        errors.append("bismark must be a mapping and bismark.align_extra a string")
    mex = cfg.get("methylation_extractor")
    if isinstance(mex, dict):
        for key in ("cx_report", "merge_cpg"):
            if not isinstance(mex.get(key), bool):
                errors.append(f"methylation_extractor.{key} must be true/false, got {mex.get(key)!r}")
        v = mex.get("buffer_frac")
        if isinstance(v, bool) or not isinstance(v, int) or v < 1:
            errors.append(f"methylation_extractor.buffer_frac must be an integer >= 1, got {v!r}")
    else:
        errors.append("methylation_extractor must be a mapping with sub-keys")
    genome = cfg.get("genome")
    if not isinstance(genome, str) or not genome.strip():
        errors.append(f"genome must be a non-empty fasta path string, got {genome!r}")
    resources = cfg.get("resources") or {}
    if not isinstance(resources, dict):
        errors.append(f"resources must be a mapping of rule -> {{threads, mem_mb, runtime_min}}, got {resources!r}")
    else:
        for rule, entry in resources.items():
            if not isinstance(entry, dict):
                errors.append(f"resources.{rule} must be a mapping, got {entry!r}")
                continue
            for field, value in entry.items():
                if field not in ("threads", "mem_mb", "runtime_min"):
                    errors.append(f"resources.{rule}.{field}: unknown field (supported: threads/mem_mb/runtime_min)")
                elif isinstance(value, bool) or not isinstance(value, int) or value < 1:
                    errors.append(f"resources.{rule}.{field} must be an integer >= 1, got {value!r}")
    if errors:
        raise WorkflowError(
            f"config validation failed ({len(errors)} issues):\n  " + "\n  ".join(errors)
        )
    if "/path/to/" in genome:
        warnings.append(f"genome is still a placeholder: {genome} (fill real paths in the project config)")
    for w in warnings:
        print(f"[config warning] {w}")


validate_config(config)

TRIM_ENABLED = bool(config["trim"]["enabled"])

# ---------------------------------------------------------------------
# Output redirection (aligned with rna-seq/chip/seclip-seq/srna-seq):
# raw inputs (1.rawdata/) stay at the project working-directory root;
# every derived artifact lives under the configurable results root
# (default "results").
# ---------------------------------------------------------------------
RD = str(config.get("results_dir", "results")).rstrip(os.sep) + os.sep


def R(path=""):
    """Return a path under the configured results directory."""
    return f"{RD}{path}"


def raw_reads(sample):
    """Detect paired-end first, then single-end raw FASTQ files for one
    declared sample, using exact names under 1.rawdata/. Returns
    ("PE"|"SE", [files])."""
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


def layout_of(sample):
    """Library layout of one declared sample: "PE" or "SE"."""
    return raw_reads(sample)[0]


def raw_se(wildcards):
    """Input function: the raw single-end FASTQ of one declared sample."""
    layout, files = raw_reads(wildcards.sample)
    if layout != "SE":
        raise ValueError(
            f"sample {wildcards.sample}: detected {layout} data but the SE trim rule was selected"
        )
    return files[0]


def raw_pe(wildcards):
    """Input function (use with unpack): the two raw paired-end FASTQs of
    one declared sample."""
    layout, files = raw_reads(wildcards.sample)
    if layout != "PE":
        raise ValueError(
            f"sample {wildcards.sample}: detected {layout} data but the PE trim rule was selected"
        )
    return {"fq1": files[0], "fq2": files[1]}


def align_input(sample):
    """Alignment input FASTQ list for one sample: the Trim Galore outputs
    when trim.enabled is true, otherwise the raw reads."""
    layout, raw_files = raw_reads(sample)
    if not TRIM_ENABLED:
        return list(raw_files)
    if layout == "PE":
        return [
            R(f"2.cleandata/{sample}_1_val_1.fq.gz"),
            R(f"2.cleandata/{sample}_2_val_2.fq.gz"),
        ]
    return [R(f"2.cleandata/{sample}_trimmed.fq.gz")]


def align_input_arg(sample):
    """Bismark command-line argument string for one sample's FASTQ input:
    --1/--2 pairs for PE, a positional path for SE."""
    files = align_input(sample)
    if len(files) == 2:
        return f"--1 {files[0]} --2 {files[1]}"
    return files[0]


# ---------------------------------------------------------------------
# Per-rule scheduler resources (same model as the sibling projects).
# ---------------------------------------------------------------------
RESOURCE_DEFAULTS = {
    "software_versions": {"threads": 1, "mem_mb": 1024, "runtime_min": 10},
    "trim": {"threads": 4, "mem_mb": 4096, "runtime_min": 120},
    "fastqc": {"threads": 2, "mem_mb": 2048, "runtime_min": 30},
    "multiqc": {"threads": 2, "mem_mb": 4096, "runtime_min": 30},
    "bismark_genome_prep": {"threads": 12, "mem_mb": 16000, "runtime_min": 240},
    "bismark_align": {"threads": 12, "mem_mb": 32000, "runtime_min": 720},
    "deduplicate": {"threads": 4, "mem_mb": 16000, "runtime_min": 240},
    "bam2nuc_genome": {"threads": 2, "mem_mb": 8000, "runtime_min": 60},
    "bam2nuc_sample": {"threads": 2, "mem_mb": 8000, "runtime_min": 60},
    "methylation_extractor": {"threads": 4, "mem_mb": 32000, "runtime_min": 1440},
    "coverage2cytosine": {"threads": 4, "mem_mb": 16000, "runtime_min": 720},
    "bismark2report": {"threads": 1, "mem_mb": 4096, "runtime_min": 30},
    "bismark2summary": {"threads": 1, "mem_mb": 4096, "runtime_min": 30},
}


def _rule_resource(name, field):
    if name not in RESOURCE_DEFAULTS:
        raise WorkflowError(f"no resource defaults for rule {name!r}")
    defaults = RESOURCE_DEFAULTS[name]
    overrides = (config.get("resources") or {}).get(name) or {}
    return int(overrides.get(field, defaults[field]))


def rthreads(name):
    """Rule threads, capped by the legacy global config["threads"]."""
    requested = max(1, _rule_resource(name, "threads"))
    legacy = config.get("threads")
    if legacy not in (None, ""):
        return min(requested, max(1, int(legacy)))
    return requested


def rmem(name):
    """Memory requested by one rule, in MB."""
    return max(1, _rule_resource(name, "mem_mb"))


def rruntime(name):
    """Walltime requested by one rule, in minutes."""
    return max(1, _rule_resource(name, "runtime_min"))


def rruntime_sec(name):
    """Walltime requested by one rule, in seconds (for PBS/SGE)."""
    return rruntime(name) * 60


# ---------------------------------------------------------------------
# Target aggregation.
#
# Derived Bismark naming convention (keep the rule modules in sync):
# deduplicate_bismark inherits the full alignment basename, so the dedup
# BAM is {sample}.bam.deduplicated.bam; the methylation extractor output
# inherits that deduplicated basename again, i.e.
# {sample}.bam.deduplicated.bismark.cov.gz (+ the _seq_context.html
# report), and coverage2cytosine --merge_CpG writes
# {sample}.CpG_merged.tsv.gz. The C3-C5 modules must emit exactly these
# paths.
# ---------------------------------------------------------------------
TARGETS = [
    R("5.QC/bismark2summary.html"),
    R("5.QC/multiqc/multiqc_report.html"),
    R("5.QC/software_versions.yaml"),
]
TARGETS += [expand(R("4.dedup/{sample}.bam.deduplicated.bam"), sample=SAMPLES)]
TARGETS += [expand(R("5.methylation/{sample}/{sample}.bam.deduplicated.bismark.cov.gz"), sample=SAMPLES)]
TARGETS += [expand(R("5.methylation/{sample}/{sample}_seq_context.html"), sample=SAMPLES)]
if config["methylation_extractor"]["merge_cpg"]:
    TARGETS += [expand(R("5.methylation/{sample}/{sample}.CpG_merged.tsv.gz"), sample=SAMPLES)]
