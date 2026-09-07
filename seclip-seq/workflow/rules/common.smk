# ---------------------------------------------------------------------
# Shared definitions: sample-table parsing (optional condition/role
# columns), config validation, path and resource helpers, CLIPper
# resolution, and target aggregation.
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


_GROUPED_HEADER = ["sample_id", "condition", "role"]


def _read_sample_table(path):
    """Parse the sample table; returns (samples, conditions, roles).

    Accepted headers are exactly ["sample_id"] (classic single-column table)
    or exactly ["sample_id", "condition", "role"] (the optional grouping
    columns that drive the reproducible_peaks stage). With the grouping
    columns present, "condition" follows the sample-id name rules (it becomes
    a file name) and "role" is restricted to ip|input; every row must fill
    both. When the columns are absent, the conditions/roles dicts come back
    empty (so lookups via SAMPLE_CONDITIONS.get(sample) yield None)."""
    samples, conditions, roles = [], {}, {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        grouped = reader.fieldnames == _GROUPED_HEADER
        if reader.fieldnames != ["sample_id"] and not grouped:
            raise WorkflowError(
                f"Sample table {path} header must be exactly ['sample_id'] or exactly "
                f"['sample_id', 'condition', 'role'] (got {reader.fieldnames}); "
                "see config/samples.csv"
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
            if grouped:
                cond = (row["condition"] or "").strip()
                role = (row["role"] or "").strip()
                if not cond:
                    raise WorkflowError(f"Sample table line {lineno}: condition must not be empty")
                if not _NAME_RE.match(cond) or "__" in cond:
                    raise WorkflowError(
                        f"Sample table line {lineno}: condition={cond!r} contains illegal "
                        "characters; only alphanumerics and . _ - are allowed "
                        "(no leading '-', no '__')"
                    )
                if role not in ("ip", "input"):
                    raise WorkflowError(
                        f"Sample table line {lineno}: role={role!r} must be 'ip' or 'input'"
                    )
                conditions[sid] = cond
                roles[sid] = role
    if not samples:
        raise WorkflowError(f"Sample table {path} has no data rows")
    return samples, conditions, roles


def load_sample_table(path):
    """Parse the sample table; returns the ordered list of unique sample ids
    (see _read_sample_table for the accepted header shapes)."""
    return _read_sample_table(path)[0]


SAMPLES, SAMPLE_CONDITIONS, SAMPLE_ROLES = _read_sample_table(
    _resolve_sample_table(config["SampleListFile"]))
SAMPLE_WILDCARD = "(?:" + "|".join(re.escape(s) for s in SAMPLES) + ")"


def validate_config(cfg):
    """Parse-time hard validation (mirrors chip_cuttag_atac_faire): required
    keys, types and value ranges aggregated into a single report. Reference
    existence only warns (--lint/dry-run often run where references are
    absent)."""
    errors, warnings = [], []
    required = ("SampleListFile", "results_dir", "threads", "filter_repeats",
                "genome", "gtf", "umi", "cutadapt", "star", "callpeak")
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
    # Optional v0.2 stage switches (both default off; absent sections are
    # valid and mean "disabled", which keeps older project configs working).
    reproducible = cfg.get("reproducible_peaks")
    annotate = cfg.get("annotate_peaks")
    for name, section in (("reproducible_peaks", reproducible), ("annotate_peaks", annotate)):
        if section is not None and not isinstance(section, dict):
            errors.append(f"{name} must be a mapping with sub-keys, got {section!r}")
    if isinstance(reproducible, dict):
        if not isinstance(reproducible.get("enabled"), bool):
            errors.append(f"reproducible_peaks.enabled must be true/false, got {reproducible.get('enabled')!r}")
        v = reproducible.get("min_replicates", 2)
        if isinstance(v, bool) or not isinstance(v, int) or v < 1:
            errors.append(f"reproducible_peaks.min_replicates must be an integer >= 1, got {v!r}")
        for key in ("input_control", "filter_by_input"):
            v = reproducible.get(key, False)
            if not isinstance(v, bool):
                errors.append(f"reproducible_peaks.{key} must be true/false, got {v!r}")
    if isinstance(annotate, dict):
        if not isinstance(annotate.get("enabled"), bool):
            errors.append(f"annotate_peaks.enabled must be true/false, got {annotate.get('enabled')!r}")
    pureclip_on = isinstance(cfg.get("callpeak"), dict) and bool(cfg["callpeak"].get("pureclip"))
    reproducible_on = isinstance(reproducible, dict) and bool(reproducible.get("enabled"))
    annotate_on = isinstance(annotate, dict) and bool(annotate.get("enabled"))
    input_control = isinstance(reproducible, dict) and bool(reproducible.get("input_control", False))
    filter_by_input = isinstance(reproducible, dict) and bool(reproducible.get("filter_by_input", False))
    if input_control and not reproducible_on:
        warnings.append(
            "reproducible_peaks.input_control=true has no effect while "
            "reproducible_peaks.enabled=false (the consensus stage is off)"
        )
    if filter_by_input and not input_control:
        errors.append(
            "reproducible_peaks.filter_by_input=true requires "
            "reproducible_peaks.input_control=true (there is no input background "
            "to filter against)"
        )
    if reproducible_on and SAMPLE_CONDITIONS:
        # Cross-check the condition/role table: every declared condition needs
        # at least one ip sample (a consensus must exist to flag), and with
        # input_control every ip condition should ideally also carry inputs.
        conds_all = sorted(set(SAMPLE_CONDITIONS.values()))
        conds_with_ip = set()
        for cond in conds_all:
            has_ip = any(SAMPLE_CONDITIONS.get(s) == cond and SAMPLE_ROLES.get(s) == "ip"
                         for s in SAMPLES)
            has_input = any(SAMPLE_CONDITIONS.get(s) == cond and SAMPLE_ROLES.get(s) == "input"
                            for s in SAMPLES)
            if has_ip:
                conds_with_ip.add(cond)
            elif has_input:
                errors.append(
                    f"reproducible_peaks: condition {cond!r} has input sample(s) but no ip "
                    "sample; there is no ip consensus to build or flag"
                )
        if input_control and reproducible_on:
            for cond in sorted(conds_with_ip):
                has_input = any(SAMPLE_CONDITIONS.get(s) == cond and SAMPLE_ROLES.get(s) == "input"
                                for s in SAMPLES)
                if not has_input:
                    warnings.append(
                        f"reproducible_peaks: condition {cond!r} has ip sample(s) but no input "
                        "sample; the input background is absent and its consensus stays "
                        "unflagged (plain BED4)"
                    )
    if input_control and reproducible_on and not SAMPLE_CONDITIONS:
        errors.append(
            "reproducible_peaks.input_control=true requires a sample table with the optional "
            "condition and role columns (header sample_id,condition,role); got the "
            "single-column table"
        )
    if reproducible_on:
        if not SAMPLE_CONDITIONS:
            errors.append(
                "reproducible_peaks.enabled=true requires a sample table with the optional "
                "condition and role columns (header sample_id,condition,role); got the "
                "single-column table"
            )
        if not pureclip_on:
            errors.append(
                "reproducible_peaks.enabled=true requires callpeak.pureclip=true "
                "(the consensus groups the per-sample PureCLIP beds)"
            )
        if SAMPLE_CONDITIONS and isinstance(reproducible.get("min_replicates"), int) \
                and not isinstance(reproducible.get("min_replicates"), bool):
            ip_counts = {}
            for sid, cond in SAMPLE_CONDITIONS.items():
                if SAMPLE_ROLES.get(sid) == "ip":
                    ip_counts[cond] = ip_counts.get(cond, 0) + 1
            for cond in sorted(set(SAMPLE_CONDITIONS.values())):
                n = ip_counts.get(cond, 0)
                if n < reproducible["min_replicates"]:
                    errors.append(
                        f"reproducible_peaks: condition {cond!r} has {n} ip sample(s) but "
                        f"min_replicates={reproducible['min_replicates']}; consensus support "
                        "can never reach the threshold"
                    )
    if annotate_on and not pureclip_on:
        errors.append(
            "annotate_peaks.enabled=true requires callpeak.pureclip=true "
            "(annotation consumes the per-sample PureCLIP peak sets)"
        )
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
          "is configured (software.yaml paths.clipper); CLIPper peak calling is "
          "skipped and only PureCLIP runs")

# Optional v0.2 stages (both default off; see config.yaml reproducible_peaks /
# annotate_peaks). The consensus stage groups the ip-role samples of one
# condition; its condition list drives the {condition} wildcard.
REPRODUCIBLE_PEAKS_ENABLED = bool((config.get("reproducible_peaks") or {}).get("enabled", False))
ANNOTATE_PEAKS_ENABLED = bool((config.get("annotate_peaks") or {}).get("enabled", False))
MIN_REPLICATES = int((config.get("reproducible_peaks") or {}).get("min_replicates", 2))
CONDITIONS = sorted({c for s, c in SAMPLE_CONDITIONS.items() if SAMPLE_ROLES.get(s) == "ip"})
CONDITION_WILDCARD = (
    "(?:" + "|".join(re.escape(c) for c in CONDITIONS) + ")" if CONDITIONS else "(?!x)x"
)
# Input-control extension of the consensus stage (TODO item 3): with
# input_control, role=input samples run the peak-calling chain too and their
# PureCLIP beds form a per-condition background union against which the ip
# consensus is flagged (and, with filter_by_input, filtered). Both switches
# default to false, so W7-era configs keep their exact DAG.
INPUT_CONTROL_ENABLED = bool((config.get("reproducible_peaks") or {}).get("input_control", False))
FILTER_BY_INPUT = bool((config.get("reproducible_peaks") or {}).get("filter_by_input", False))

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


def condition_ip_samples(condition):
    """Ordered (sample-table order) list of the ip-role samples of one
    condition; used by the reproducible_peaks consensus rule. Helper lives
    here to keep the rule modules rules-only."""
    return [s for s in SAMPLES
            if SAMPLE_CONDITIONS.get(s) == condition and SAMPLE_ROLES.get(s) == "ip"]


def condition_input_samples(condition):
    """Ordered (sample-table order) list of the input-role samples of one
    condition; feeds the input_background union consensus (empty for
    single-column tables, where no roles are declared)."""
    return [s for s in SAMPLES
            if SAMPLE_CONDITIONS.get(s) == condition and SAMPLE_ROLES.get(s) == "input"]


# Conditions of the ip-role samples, split by whether they also carry input
# samples (only those get a background consensus and the flagged consensus
# shape; the others keep the plain BED4 consensus).
CONDITIONS_WITH_INPUT = [c for c in CONDITIONS if condition_input_samples(c)]
CONDITIONS_WITHOUT_INPUT = [c for c in CONDITIONS if not condition_input_samples(c)]

# Samples that receive per-sample PureCLIP peak calling. With a single-column
# table (or while reproducible_peaks is disabled) this is every declared
# sample; with the grouping columns and reproducible_peaks enabled, the
# ip-role samples only — unless input_control is enabled, in which case the
# input controls run peak calling too so their beds can form the background.
if SAMPLE_CONDITIONS and REPRODUCIBLE_PEAKS_ENABLED and not INPUT_CONTROL_ENABLED:
    PEAKCALL_SAMPLES = [s for s in SAMPLES if SAMPLE_ROLES.get(s) == "ip"]
else:
    PEAKCALL_SAMPLES = SAMPLES


def consensus_annotated_bed(condition):
    """Peak BED consumed by the consensus annotation rule: the filtered
    consensus when filter_by_input is active and the condition carries input
    samples, the (flagged or plain) consensus BED otherwise."""
    if INPUT_CONTROL_ENABLED and FILTER_BY_INPUT and condition in CONDITIONS_WITH_INPUT:
        return R(f"6.reproducible_peaks/{condition}.consensus.filtered.bed")
    return R(f"6.reproducible_peaks/{condition}.consensus.bed")


def consensus_peak_cols(condition):
    """--peak-cols for the consensus annotation: the flagged consensus keeps
    the W7 BED4 columns and appends the binary in_input_background flag as
    column 5; every other consensus shape is plain BED4."""
    if INPUT_CONTROL_ENABLED and not FILTER_BY_INPUT and condition in CONDITIONS_WITH_INPUT:
        return 5
    return 4


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
    "consensus_peaks": {"threads": 1, "mem_mb": 2048, "runtime_min": 30},
    "input_background": {"threads": 1, "mem_mb": 2048, "runtime_min": 30},
    "flag_input_background": {"threads": 1, "mem_mb": 2048, "runtime_min": 30},
    "filter_input_background": {"threads": 1, "mem_mb": 2048, "runtime_min": 30},
    "gtf_gene_regions": {"threads": 1, "mem_mb": 4096, "runtime_min": 60},
    "annotate_peaks": {"threads": 1, "mem_mb": 4096, "runtime_min": 60},
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
    # PEAKCALL_SAMPLES excludes role=input samples unless input_control is
    # enabled (see the definition in common.smk above).
    TARGETS += [expand(R("5.callpeak/{sample}.pureclip.bed"), sample=PEAKCALL_SAMPLES)]
if CLIPPER_ENABLED:
    TARGETS += [expand(R("5.callpeak/{sample}.clipper.peakClusters.bed"), sample=SAMPLES)]
if REPRODUCIBLE_PEAKS_ENABLED:
    # The final consensus BED is produced by the flagging rule for conditions
    # with input samples (input_control on) and by the consensus rule itself
    # for every other condition (rules/consensus.smk owns the split).
    TARGETS += [expand(R("6.reproducible_peaks/{condition}.consensus.bed"), condition=CONDITIONS)]
    if INPUT_CONTROL_ENABLED:
        TARGETS += [expand(R("6.reproducible_peaks/{condition}.input_background.bed"),
                           condition=CONDITIONS_WITH_INPUT)]
        if FILTER_BY_INPUT:
            TARGETS += [expand(R("6.reproducible_peaks/{condition}.consensus.filtered.bed"),
                               condition=CONDITIONS_WITH_INPUT)]
    if ANNOTATE_PEAKS_ENABLED:
        TARGETS += [expand(R("6.annotation/{condition}.consensus.annotation.tsv"), condition=CONDITIONS)]
if ANNOTATE_PEAKS_ENABLED and bool(config["callpeak"]["pureclip"]):
    TARGETS += [expand(R("6.annotation/{sample}.annotation.tsv"), sample=PEAKCALL_SAMPLES)]
