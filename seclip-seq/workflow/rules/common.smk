# ---------------------------------------------------------------------
# Shared definitions: sample-table parsing, config validation, path and
# resource helpers, CLIPper resolution, and target aggregation.
# Included first by workflow/Snakefile; every rules/*.smk uses the names
# defined here. BASE_DIR / WORKFLOW_DIR come from the Snakefile.
# ---------------------------------------------------------------------
import csv
import os
import re

from snakemake.exceptions import WorkflowError

# Sample names allow alphanumerics . _ - (no "__", no leading '-'); they
# become file names and STAR/bowtie arguments.
_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")

# CLIPper is an external legacy install resolved from software.yaml via
# run.sh; empty means the callpeak_clipper rule is skipped entirely.
CLIPPER = os.environ.get("SECLIP_TOOL_CLIPPER", "")
FILTER_REPEATS = None  # set after validate_config below


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
    """Parse-time hard validation (mirrors chip_cuttag_atac_faire): required
    keys, types and value ranges aggregated into a single report. Reference
    existence only warns (--lint/dry-run often run where references are
    absent)."""
    errors, warnings = [], []
    required = ("SampleListFile", "results_dir", "threads", "filter_repeats",
                "umi", "cutadapt", "star", "callpeak")
    for key in required:
        if key not in cfg:
            errors.append(f"missing required config key: {key}")
    if not isinstance(cfg.get("results_dir", "results"), str):
        errors.append(f"results_dir must be a string, got {cfg.get('results_dir')!r}")
    for key in ("threads",):
        v = cfg.get(key)
        if isinstance(v, bool) or not isinstance(v, int) or v < 1:
            errors.append(f"{key} must be an integer >= 1, got {v!r}")
    if not isinstance(cfg.get("filter_repeats"), bool):
        errors.append(f"filter_repeats must be true/false, got {cfg.get('filter_repeats')!r}")
    if cfg.get("filter_repeats") and not str(cfg.get("repeats_fa") or "").strip():
        errors.append("filter_repeats=true requires a non-empty repeats_fa (species preset or project config)")
    if isinstance(cfg.get("umi"), dict):
        pattern = cfg["umi"].get("pattern", "")
        if not isinstance(pattern, str) or not re.fullmatch(r"[ACGTN]+", pattern or ""):
            errors.append(f"umi.pattern must be a non-empty string of ACGTN, got {pattern!r}")
    else:
        errors.append("umi must be a mapping with sub-keys")
    if isinstance(cfg.get("cutadapt"), dict):
        c = cfg["cutadapt"]
        for key, lo in (("min_len", 1), ("quality_cutoff", 0)):
            v = c.get(key)
            if isinstance(v, bool) or not isinstance(v, int) or v < lo:
                errors.append(f"cutadapt.{key} must be an integer >= {lo}, got {v!r}")
        try:
            er = float(c.get("error_rate"))
            if not 0 < er <= 1:
                errors.append(f"cutadapt.error_rate must be in (0, 1], got {c.get('error_rate')!r}")
        except (TypeError, ValueError):
            errors.append(f"cutadapt.error_rate must be numeric, got {c.get('error_rate')!r}")
        adapters = c.get("adapters")
        if not isinstance(adapters, list) or not adapters:
            errors.append("cutadapt.adapters must be a non-empty list")
        else:
            for i, a in enumerate(adapters):
                if not isinstance(a, str) or not re.fullmatch(r"[ACGT]+", a or ""):
                    errors.append(f"cutadapt.adapters[{i}] must be an ACGT string, got {a!r}")
    else:
        errors.append("cutadapt must be a mapping with sub-keys")
    if isinstance(cfg.get("star"), dict):
        s = cfg["star"]
        for key in ("sjdb_overhang", "genome_sa_index_nbases", "repeats_sa_index_nbases",
                    "filter_multimap_nmax", "align_multimap_nmax"):
            v = s.get(key)
            if isinstance(v, bool) or not isinstance(v, int) or v < 1:
                errors.append(f"star.{key} must be an integer >= 1, got {v!r}")
        v = s.get("repeats_limit_ram")
        if isinstance(v, bool) or not isinstance(v, int) or v <= 0:
            errors.append(f"star.repeats_limit_ram must be a positive integer, got {v!r}")
    else:
        errors.append("star must be a mapping with sub-keys")
    if isinstance(cfg.get("callpeak"), dict):
        p = cfg["callpeak"]
        for key in ("pureclip", "clipper"):
            if not isinstance(p.get(key), bool):
                errors.append(f"callpeak.{key} must be true/false, got {p.get(key)!r}")
        if p.get("clipper") and not str(p.get("clipper_species") or "").strip():
            errors.append("callpeak.clipper=true requires a non-empty callpeak.clipper_species")
    else:
        errors.append("callpeak must be a mapping with sub-keys")
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
    if errors:
        raise WorkflowError(
            f"config validation failed ({len(errors)} issues):\n  " + "\n  ".join(errors)
        )
    for key in ("genome", "gtf", "repeats_fa"):
        value = cfg.get(key)
        if not isinstance(value, str) or not value:
            warnings.append(f"{key} is not configured (species preset should provide it)")
        elif "/path/to/" in value:
            warnings.append(f"{key} is still a placeholder: {value} (fill real paths in the project config)")
    for w in warnings:
        print(f"[config warning] {w}")


validate_config(config)

FILTER_REPEATS = bool(config["filter_repeats"])
CLIPPER_ENABLED = bool(config["callpeak"]["clipper"]) and bool(CLIPPER)
if bool(config["callpeak"]["clipper"]) and not CLIPPER:
    print("[config warning] callpeak.clipper is enabled but no CLIPper executable "
          "is configured (software.yaml tools.clipper); CLIPper peak calling is "
          "skipped and only PureCLIP runs")

# ---------------------------------------------------------------------
# Output redirection (aligned with rna-seq/chip): raw inputs (1.rawdata/)
# stay at the project working-directory root; every derived artifact
# lives under the configurable results root (default "results").
# ---------------------------------------------------------------------
RD = str(config.get("results_dir", "results")).rstrip(os.sep) + os.sep


def R(path=""):
    """Return a path under the configured results directory."""
    return f"{RD}{path}"


def raw_fastq(sample):
    """Locate the raw SE FASTQ for a declared sample (exact names)."""
    for pattern in ("1.rawdata/{0}_R1.fastq.gz", "1.rawdata/{0}_R1.fq.gz",
                    "1.rawdata/{0}.fastq.gz", "1.rawdata/{0}.fq.gz"):
        if os.path.exists(pattern.format(sample)):
            return pattern.format(sample)
    raise ValueError(
        f"sample {sample}: no raw FASTQ found under 1.rawdata/. Supported names: "
        f"{sample}_R1.fastq.gz, {sample}_R1.fq.gz, {sample}.fastq.gz, {sample}.fq.gz"
    )


def _align_input(wc):
    """Genome alignment input: repeats-unmapped reads, or the trimmed reads
    when the repeats filter is disabled (config filter_repeats=false)."""
    if FILTER_REPEATS:
        return R(f"3.align/repeats/{wc.sample}_Unmapped.out.mate1")
    return R(f"2.cleandata/{wc.sample}_clean.fqTrTr.sorted.fq.gz")


# ---------------------------------------------------------------------
# Scheduler resources per rule; precedence is a project-local resources.yaml
# (or a "resources" block in the project config), then the defaults below.
# The legacy top-level "threads" value still acts as a global thread cap.
# ---------------------------------------------------------------------
RESOURCE_DEFAULTS = {
    "software_versions": {"threads": 1, "mem_mb": 1024, "runtime_min": 10},
    "fastqc": {"threads": 2, "mem_mb": 2048, "runtime_min": 30},
    "multiqc": {"threads": 2, "mem_mb": 4096, "runtime_min": 60},
    "star_index_genome": {"threads": 12, "mem_mb": 48000, "runtime_min": 480},
    "star_index_repeats": {"threads": 12, "mem_mb": 32000, "runtime_min": 240},
    "umi_extract": {"threads": 4, "mem_mb": 8000, "runtime_min": 180},
    "cutadapt_trim": {"threads": 8, "mem_mb": 8000, "runtime_min": 180},
    "fastq_sort": {"threads": 4, "mem_mb": 8000, "runtime_min": 120},
    "star_filter_repeats": {"threads": 12, "mem_mb": 32000, "runtime_min": 360},
    "star_align": {"threads": 12, "mem_mb": 32000, "runtime_min": 360},
    "umi_dedup": {"threads": 4, "mem_mb": 16000, "runtime_min": 360},
    "read_count": {"threads": 1, "mem_mb": 2000, "runtime_min": 10},
    "callpeak_clipper": {"threads": 4, "mem_mb": 16000, "runtime_min": 240},
    "callpeak_pureclip": {"threads": 8, "mem_mb": 16000, "runtime_min": 360},
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
# Target aggregation
# ---------------------------------------------------------------------
TARGETS = [
    R("5.QC/multiqc/multiqc_report.html"),
    R("5.QC/software_versions.yaml"),
]
TARGETS += [expand(R("4.rmdup/{sample}_readnum.txt"), sample=SAMPLES)]
if config["callpeak"]["pureclip"]:
    TARGETS += [expand(R("5.callpeak/{sample}.pureclip.bed"), sample=SAMPLES)]
if CLIPPER_ENABLED:
    TARGETS += [expand(R("5.callpeak/{sample}.clipper.peakClusters.bed"), sample=SAMPLES)]
