# ---------------------------------------------------------------------
# Shared definitions: path constants, sample-table parsing, config
# validation, rule query helpers, and target aggregation.
# Included first by workflow/Snakefile; every rules/*.smk uses the names
# defined here. BASE_DIR / WORKFLOW_DIR come from the Snakefile.
# ---------------------------------------------------------------------
import csv
import os
import re

from snakemake.exceptions import WorkflowError

ASSAYS = ("chip", "cuttag", "atac", "faire")

# Sample and group names allow only alphanumerics . _ - : commas break the
# peak-list concatenation and MACS2 multi-file arguments, "__" is the FRiP
# output separator, and spaces/tabs break shell expansion.
_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")

# ---------------------------------------------------------------------
# Sample table parsing and validation
# ---------------------------------------------------------------------
REQUIRED_COLUMNS = ("sample_id", "role", "group", "seqtype", "layout", "peak_type")


def _resolve_sample_table(path):
    """Resolve the sample table path: absolute > relative to the working
    directory > relative to the repository root."""
    p = os.path.expanduser(str(path))
    if os.path.isabs(p):
        return p
    if os.path.exists(p):
        return os.path.abspath(p)
    return os.path.join(BASE_DIR, p)


def load_sample_table(path):
    """Parse the sample table; returns (sample ids, group dict, sample->seqtype map).

    Group dict layout: group -> {seqtype, peak_type, layout, treat: [], control: []}
    """
    samples = []          # deduplicated sample ids (controls are aligned too)
    seqtype_of = {}
    groups = {}

    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        missing = [c for c in REQUIRED_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise WorkflowError(
                f"Sample table {path} is missing columns: {missing}; "
                f"it must contain {list(REQUIRED_COLUMNS)}; see config/samples.csv"
            )
        for lineno, row in enumerate(reader, start=2):
            sid = (row["sample_id"] or "").strip()
            role = (row["role"] or "").strip().lower()
            grp = (row["group"] or "").strip()
            seqtype = (row["seqtype"] or "").strip().lower()
            layout = (row["layout"] or "").strip().upper()
            peak_type = (row["peak_type"] or "").strip().lower()

            if not sid or not grp:
                raise WorkflowError(f"Sample table line {lineno}: sample_id/group must not be empty")
            for label, value in (("sample_id", sid), ("group", grp)):
                if not _NAME_RE.match(value) or "__" in value:
                    raise WorkflowError(
                        f"Sample table line {lineno}: {label}={value!r} contains illegal "
                        "characters; only alphanumerics and . _ - are allowed (must not "
                        "start with '-', no consecutive '__')"
                    )
            if role not in ("treat", "control"):
                raise WorkflowError(f"Sample table line {lineno}: role must be treat or control, got {role!r}")
            if seqtype not in ASSAYS:
                raise WorkflowError(f"Sample table line {lineno}: seqtype must be one of {'/'.join(ASSAYS)}, got {seqtype!r}")
            if layout != "PE":
                raise WorkflowError(
                    f"Sample table line {lineno}: only PE (paired-end) data is supported, layout={layout!r}"
                )
            if peak_type not in ("narrow", "broad", "none"):
                raise WorkflowError(
                    f"Sample table line {lineno}: peak_type must be narrow/broad/none, got {peak_type!r}"
                )
            if seqtype in ("atac", "faire") and peak_type != "none":
                raise WorkflowError(
                    f"Sample table line {lineno}: peak_type for atac/faire must be none (the peak mode is fixed by the workflow)"
                )
            if sid in seqtype_of and seqtype_of[sid] != seqtype:
                raise WorkflowError(f"Sample table line {lineno}: sample {sid} appears in groups of different seqtypes")
            if role == "treat" and peak_type == "none" and seqtype in ("chip", "cuttag"):
                raise WorkflowError(
                    f"Sample table line {lineno}: chip/cuttag treat sample {sid} must declare narrow or broad"
                )

            if sid not in seqtype_of:
                seqtype_of[sid] = seqtype
                samples.append(sid)

            g = groups.setdefault(
                grp,
                {"seqtype": seqtype, "peak_type": peak_type, "layout": layout,
                 "treat": [], "control": []},
            )
            if (g["seqtype"], g["peak_type"], g["layout"]) != (seqtype, peak_type, layout):
                raise WorkflowError(
                    f"Sample table line {lineno}: all rows of group {grp} must share the same seqtype/peak_type/layout"
                )
            g[role].append(sid)

    bad_groups = [g for g, v in groups.items() if not v["treat"]]
    if bad_groups:
        raise WorkflowError(f"These groups have no role=treat sample: {bad_groups}")
    if not samples:
        raise WorkflowError(f"Sample table {path} has no data rows")

    return samples, groups, seqtype_of


SAMPLES, GROUPS, SEQTYPE_OF = load_sample_table(_resolve_sample_table(config["grouplist"]))

try:
    config["threads"] = int(config["threads"])
except (TypeError, ValueError):
    raise WorkflowError(
        f"config threads must be an integer (do not quote it), got {config['threads']!r}"
    )

# ATAC peak-calling mode whitelist: mode only accepts bampe / shifted; a typo
# must be caught at parse time instead of silently falling into shifted.
if config["peak"]["atac"]["mode"] not in ("bampe", "shifted"):
    raise WorkflowError(
        "config peak.atac.mode must be bampe or shifted, got "
        f"{config['peak']['atac']['mode']!r}"
    )


def validate_config(cfg):
    """Central config validation: required keys, sub-keys, types and value
    ranges, aggregated into a single report. Reference-file existence only
    warns (--lint/dry-run often run where reference files are absent)."""
    errors, warnings = [], []
    required = ("genome_fa", "gtf", "bed", "chromsize", "genome_size",
                "grouplist", "threads", "bowtie2_extra", "min_mapq",
                "region_flank", "dedup", "peak", "qc", "trim")
    for key in required:
        if key not in cfg:
            errors.append(f"missing required config key: {key}")
    if "results_dir" in cfg and not isinstance(cfg["results_dir"], str):
        errors.append(f"results_dir must be a string, got {cfg['results_dir']!r}")
    for key in ("dedup", "qc", "peak", "trim"):
        if key in cfg and not isinstance(cfg[key], dict):
            errors.append(f"{key} must be a mapping (with sub-keys), got {cfg[key]!r}")
    if isinstance(cfg.get("dedup"), dict):
        for assay in ASSAYS:
            v = cfg["dedup"].get(assay)
            if not isinstance(v, bool):
                errors.append(f"dedup.{assay} must be true/false, got {v!r}")
    if isinstance(cfg.get("qc"), dict):
        for key in ("nsc_rsc", "frip", "deeptools"):
            v = cfg["qc"].get(key)
            if not isinstance(v, bool):
                errors.append(f"qc.{key} must be true/false, got {v!r}")
    for key in ("min_mapq", "region_flank"):
        v = cfg.get(key)
        if isinstance(v, bool) or not isinstance(v, int) or v < 0:
            errors.append(f"{key} must be an integer >= 0, got {v!r}")
    if isinstance(cfg.get("peak"), dict):
        p = cfg["peak"]
        for key in ("qvalue", "broad_cutoff"):
            try:
                if not 0 < float(p[key]) <= 1:
                    errors.append(f"peak.{key} must be in (0, 1], got {p[key]!r}")
            except (TypeError, ValueError):
                errors.append(f"peak.{key} must be numeric, got {p[key]!r}")
    if isinstance(cfg.get("trim"), dict):
        t = cfg["trim"]
        for key, lo in (("quality", 0), ("stringency", 1)):
            v = t.get(key)
            if isinstance(v, bool) or not isinstance(v, int) or v < lo:
                errors.append(f"trim.{key} must be an integer >= {lo}, got {v!r}")
        try:
            er = float(t.get("error_rate"))
            if not 0 < er <= 1:
                errors.append(f"trim.error_rate must be in (0, 1], got {er!r}")
        except (TypeError, ValueError):
            errors.append(f"trim.error_rate must be numeric, got {t.get('error_rate')!r}")
        if not isinstance(t.get("extra", ""), str):
            errors.append("trim.extra must be a string")
    if errors:
        raise WorkflowError(
            f"config validation failed ({len(errors)} issues):\n  " + "\n  ".join(errors)
        )
    for key in ("genome_fa", "gtf", "bed", "chromsize"):
        p = str(cfg[key])
        if not os.path.isabs(p) and not os.path.exists(p):
            p = os.path.join(BASE_DIR, p)
        if not os.path.exists(p):
            warnings.append(f"reference file does not exist (verify before running): {key} = {cfg[key]}")
    for w in warnings:
        print(f"[config warning] {w}")


validate_config(config)

# ---------------------------------------------------------------------
# Output redirection (aligned with rna-seq): raw inputs (1.rawdata/) stay
# at the project working-directory root; every derived artifact lives
# under the configurable results root (default "results").
# ---------------------------------------------------------------------
RD = str(config.get("results_dir", "results")).rstrip(os.sep) + os.sep


def R(path=""):
    """Return a path under the configured results directory."""
    return f"{RD}{path}"


# ---------------------------------------------------------------------
# Shared query helpers (used by the rules/*.smk modules)
# ---------------------------------------------------------------------

def assay_needs_dedup(seqtype):
    """Whether picard dedup runs for an assay; CUT&Tag keeps PCR duplicates."""
    return bool(config["dedup"][seqtype])


def sample_bam(sample):
    suffix = "rmdup.bam" if assay_needs_dedup(SEQTYPE_OF[sample]) else "sorted.bam"
    return f"{RD}3.align/bowtie2/{sample}_{suffix}"


def group_bams(group, role):
    return [sample_bam(s) for s in GROUPS[group][role]]


def group_control_arg(wc):
    bams = group_bams(wc.group, "control")
    return "-c " + ",".join(bams) if bams else ""


def group_peak_file(group):
    if GROUPS[group]["peak_type"] == "broad":
        return f"{RD}4.peak/{group}_peaks.broadPeak"
    return f"{RD}4.peak/{group}_peaks.narrowPeak"


def _groups_of(assay, peak_type=None):
    return [g for g, v in GROUPS.items()
            if v["seqtype"] == assay and (peak_type is None or v["peak_type"] == peak_type)]


def _group_regex(groups):
    """Compile a group-name list into a wildcard constraint; an empty list
    yields a never-matching regex."""
    if not groups:
        return "(?!x)x"
    return "(?:" + "|".join(re.escape(g) for g in groups) + ")"


# ---------------------------------------------------------------------
# Per-rule scheduler resources (unified model, aligned with rna-seq).
# Precedence: config/resources.yaml (or a project-local resources.yaml /
# a "resources" block in the project config) > RESOURCE_DEFAULTS below.
# The legacy top-level "threads" value still acts as a global thread cap.
# ---------------------------------------------------------------------
RESOURCE_DEFAULTS = {
    "software_versions": {"threads": 1, "mem_mb": 1024, "runtime_min": 10},
    "trim_adapter": {"threads": 4, "mem_mb": 4096, "runtime_min": 60},
    "fastqc": {"threads": 2, "mem_mb": 2048, "runtime_min": 30},
    "multiqc": {"threads": 1, "mem_mb": 4096, "runtime_min": 30},
    "bowtie2_index": {"threads": 8, "mem_mb": 8192, "runtime_min": 120},
    "bowtie2_mapping": {"threads": 8, "mem_mb": 16384, "runtime_min": 240},
    "dedup": {"threads": 2, "mem_mb": 8192, "runtime_min": 120},
    "callpeak_narrow": {"threads": 1, "mem_mb": 8192, "runtime_min": 180},
    "callpeak_broad": {"threads": 1, "mem_mb": 8192, "runtime_min": 180},
    "callpeak_atac": {"threads": 1, "mem_mb": 8192, "runtime_min": 180},
    "bigwig": {"threads": 1, "mem_mb": 4096, "runtime_min": 60},
    "peak_annotation": {"threads": 1, "mem_mb": 8192, "runtime_min": 120},
    "frip": {"threads": 1, "mem_mb": 4096, "runtime_min": 60},
    "frip_summary": {"threads": 1, "mem_mb": 1024, "runtime_min": 10},
    "deeptools_multibamsummary": {"threads": 2, "mem_mb": 8192, "runtime_min": 120},
    "deeptools_correlation": {"threads": 1, "mem_mb": 4096, "runtime_min": 30},
    "deeptools_pca": {"threads": 1, "mem_mb": 4096, "runtime_min": 30},
    "deeptools_fingerprint": {"threads": 2, "mem_mb": 8192, "runtime_min": 60},
    "deeptools_fragmentsize": {"threads": 2, "mem_mb": 8192, "runtime_min": 60},
    "deeptools_profile": {"threads": 2, "mem_mb": 8192, "runtime_min": 120},
    "spp_crosscorr": {"threads": 2, "mem_mb": 8192, "runtime_min": 180},
    "spp_summary": {"threads": 1, "mem_mb": 1024, "runtime_min": 10},
}


def _rule_resource(name, field):
    if name not in RESOURCE_DEFAULTS:
        raise WorkflowError(f"unknown rule name in resource lookup: {name!r}")
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
    """Walltime requested by one rule, in seconds (for SGE h_rt)."""
    return rruntime(name) * 60


# ---------------------------------------------------------------------
# Target aggregation
# ---------------------------------------------------------------------

BAM_TARGETS = [f"{RD}3.align/bowtie2/{s}_sorted.bam" for s in SAMPLES]
BAM_TARGETS += [f"{RD}3.align/bowtie2/{s}_rmdup.bam"
                for s in SAMPLES if assay_needs_dedup(SEQTYPE_OF[s])]

PEAK_TARGETS = [group_peak_file(g) for g in GROUPS]
BW_TARGETS = [f"{RD}4.peak/{g}_FE.bw" for g in GROUPS]

QC_TARGETS = [R("2.cleandata/fastqc/multiqc/multiqc_report.html")]
if config["qc"]["nsc_rsc"]:
    QC_TARGETS += [R("5.QC/spp/NSC_RSC_mqc.tsv")]
if config["qc"]["frip"]:
    QC_TARGETS += [R("5.QC/frip/FRiP_summary.tsv")]
if config["qc"]["deeptools"]:
    QC_TARGETS += [
        R("5.QC/deeptools/multiBamSummary.npz"),
        R("5.QC/deeptools/heatmap_SpearmanCorr_readCounts.png"),
        R("5.QC/deeptools/PCA_readCounts.png"),
        R("5.QC/deeptools/fingerprints.png"),
        R("5.QC/deeptools/fragmentsize.png"),
        R("5.QC/deeptools/profile_scaled.png"),
    ]

# Software-version record: generated by the software_versions rule in
# meta.smk (no input dependency, scheduled freely within the DAG); always
# part of rule all.
VERSION_TARGETS = [R("5.QC/software_versions.yaml")]
