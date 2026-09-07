# ---------------------------------------------------------------------
# Shared definitions: sample-table parsing, config validation, cascade
# helpers, and target aggregation. Included first by workflow/Snakefile;
# every rules/*.smk uses the names defined here. BASE_DIR / WORKFLOW_DIR
# come from the Snakefile.
# ---------------------------------------------------------------------
import csv
import os
import re

from snakemake.exceptions import WorkflowError

# Sample and cascade class names allow alphanumerics . _ - (no "__",
# no leading '-'); both become file/directory names and bowtie arguments.
_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")

SCRIPTS = os.path.join(WORKFLOW_DIR, "scripts")


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
    """Parse-time hard validation: required keys, cascade entry shapes,
    types and value ranges aggregated into a single report. Reference
    existence only warns (lint/dry-run often runs without references)."""
    errors, warnings = [], []
    for key in ("SampleListFile", "results_dir", "threads", "trim", "cascade",
                "genome", "bowtie"):
        if key not in cfg:
            errors.append(f"missing required config key: {key}")
    if not isinstance(cfg.get("results_dir", "results"), str):
        errors.append(f"results_dir must be a string, got {cfg.get('results_dir')!r}")
    v = cfg.get("threads")
    if isinstance(v, bool) or not isinstance(v, int) or v < 1:
        errors.append(f"threads must be an integer >= 1, got {v!r}")
    cascade = cfg.get("cascade")
    if not isinstance(cascade, list) or not cascade:
        errors.append("cascade must be a non-empty ordered list of {name, fasta} entries")
    else:
        names = []
        for i, entry in enumerate(cascade):
            if not isinstance(entry, dict):
                errors.append(f"cascade[{i}] must be a mapping, got {entry!r}")
                continue
            name = str(entry.get("name") or "").strip()
            if not _NAME_RE.match(name) or "__" in name:
                errors.append(f"cascade[{i}].name={name!r} is not a legal class name")
            if name in names:
                errors.append(f"cascade[{i}].name={name!r} duplicates an earlier entry")
            names.append(name)
            fasta = entry.get("fasta", "")
            if not isinstance(fasta, str):
                errors.append(f"cascade[{i}].fasta must be a string, got {fasta!r}")
        configured = [str(e.get("fasta") or "").strip() for e in cascade
                      if isinstance(e, dict)]
        genome_fasta = str((cfg.get("genome") or {}).get("fasta") or "").strip()
        if not any(configured) and not genome_fasta:
            errors.append("neither any cascade class nor the genome has a configured fasta — nothing to align")
    if isinstance(cfg.get("trim"), dict):
        t = cfg["trim"]
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
    bowtie = cfg.get("bowtie")
    if not isinstance(bowtie, dict) or not isinstance(bowtie.get("extra", ""), str):
        errors.append("bowtie must be a mapping and bowtie.extra a string")
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
    for entry in cascade or []:
        if isinstance(entry, dict):
            fasta = str(entry.get("fasta") or "").strip()
            if not fasta:
                warnings.append(f"cascade class {entry.get('name')!r} has no fasta and is skipped")
            elif "/path/to/" in fasta:
                warnings.append(f"cascade.{entry.get('name')}.fasta is still a placeholder: {fasta}")
    genome_fasta = str((cfg.get("genome") or {}).get("fasta") or "").strip()
    if genome_fasta and "/path/to/" in genome_fasta:
        warnings.append(f"genome.fasta is still a placeholder: {genome_fasta}")
    for w in warnings:
        print(f"[config warning] {w}")


validate_config(config)


def apply_species_presets(cfg, preset):
    """Fill empty reference paths from one species preset, in place.

    A cascade entry whose fasta is empty inherits preset["cascade"][name];
    genome.fasta falls back to preset["genome"]. An empty preset (species
    "none") is a no-op. This mirrors the inline species merge in
    workflow/Snakefile -- which must stay there because it runs before
    this file is included (the parse-time validation above consumes the
    merged references); tests/test_common.py pins the two together.
    """
    for entry in cfg.get("cascade") or []:
        filled = (preset.get("cascade") or {}).get(entry.get("name"), "")
        if not str(entry.get("fasta") or "").strip() and filled:
            entry["fasta"] = filled
    if not str((cfg.get("genome") or {}).get("fasta") or "").strip() and preset.get("genome"):
        cfg["genome"] = dict(cfg.get("genome") or {}, fasta=preset["genome"])


def cascade_classes():
    """Ordered list of cascade class names with a configured fasta."""
    return [str(c["name"]).strip() for c in (config.get("cascade") or [])
            if str(c.get("fasta") or "").strip()]


def cascade_fasta(name):
    """Reference fasta of one cascade class."""
    for entry in config.get("cascade") or []:
        if str(entry["name"]).strip() == name:
            return entry["fasta"]
    raise WorkflowError(f"cascade class {name!r} is not configured")


def cascade_prev(name):
    """Previous class in the cascade, or None for the first."""
    classes = cascade_classes()
    i = classes.index(name)
    return classes[i - 1] if i > 0 else None


def cascade_input(wc):
    """Input FASTQ of one cascade stage: previous stage's unmapped reads,
    or the trimmed reads for the first stage."""
    prev = cascade_prev(wc.klass)
    if prev is None:
        return R(f"2.cleandata/{wc.sample}_trimmed.fq.gz")
    return R(f"3.align/filter/{prev}/{wc.sample}_unmapped.fq")


_CLASSES = cascade_classes()
# Wildcard over the classes that actually run; an empty cascade yields a
# never-matching constraint (no filter jobs).
KLASS_WILDCARD = "(?:" + "|".join(re.escape(c) for c in _CLASSES) + ")" if _CLASSES else "(?!x)x"
_GENOME_CONFIGURED = bool(str((config.get("genome") or {}).get("fasta") or "").strip())
# bowtie index targets cover the running classes plus the genome.
_INDEX_WILDCARD = "(?:" + "|".join(re.escape(c) for c in _CLASSES + (["genome"] if _GENOME_CONFIGURED else [])) + ")" \
    if (_CLASSES or _GENOME_CONFIGURED) else "(?!x)x"

# ---------------------------------------------------------------------
# Output redirection (aligned with rna-seq/chip/seclip-seq).
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


# ---------------------------------------------------------------------
# Per-rule scheduler resources (same model as the sibling projects).
# ---------------------------------------------------------------------
RESOURCE_DEFAULTS = {
    "software_versions": {"threads": 1, "mem_mb": 1024, "runtime_min": 10},
    "trim": {"threads": 4, "mem_mb": 4096, "runtime_min": 60},
    "fastqc": {"threads": 2, "mem_mb": 2048, "runtime_min": 30},
    "multiqc": {"threads": 2, "mem_mb": 4096, "runtime_min": 30},
    "bowtie_index": {"threads": 8, "mem_mb": 8192, "runtime_min": 120},
    "cascade_stage": {"threads": 4, "mem_mb": 8192, "runtime_min": 120},
    "genome_align": {"threads": 8, "mem_mb": 8192, "runtime_min": 180},
    "count_stage": {"threads": 1, "mem_mb": 2048, "runtime_min": 30},
    "merge_counts": {"threads": 1, "mem_mb": 4096, "runtime_min": 30},
    "cascade_summary": {"threads": 1, "mem_mb": 2048, "runtime_min": 30},
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
    R("4.expression/all_classes_counts.tsv"),
    R("5.QC/cascade_summary.tsv"),
    R("5.QC/multiqc/multiqc_report.html"),
    R("5.QC/software_versions.yaml"),
]
if _CLASSES:
    TARGETS += [expand(R("4.expression/{klass}/{klass}_counts.tsv"), klass=_CLASSES)]
if _GENOME_CONFIGURED:
    TARGETS += [expand(R("3.align/genome/{sample}.sam"), sample=SAMPLES)]
