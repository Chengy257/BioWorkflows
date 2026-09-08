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
    """Parse the sample table; returns (sample ids, group dict, sample->seqtype
    map, treat sample->condition, treat sample->batch).

    Group dict layout: group -> {seqtype, peak_type, layout, condition, treat: [], control: []}
    condition/batch are OPTIONAL columns (differential binding, v0.5 Phase 3):
    absent columns leave the maps empty and nothing changes for legacy tables.
    """
    samples = []          # deduplicated sample ids (controls are aligned too)
    seqtype_of = {}
    groups = {}
    treat_condition_of = {}
    treat_batch_of = {}

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
            condition = (row.get("condition") or "").strip()
            batch = (row.get("batch") or "").strip()

            if not sid or not grp:
                raise WorkflowError(f"Sample table line {lineno}: sample_id/group must not be empty")
            for label, value in (("sample_id", sid), ("group", grp),
                                 ("condition", condition), ("batch", batch)):
                if value and (not _NAME_RE.match(value) or "__" in value):
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
                 "condition": None, "treat": [], "control": []},
            )
            if (g["seqtype"], g["peak_type"], g["layout"]) != (seqtype, peak_type, layout):
                raise WorkflowError(
                    f"Sample table line {lineno}: all rows of group {grp} must share the same seqtype/peak_type/layout"
                )
            g[role].append(sid)
            if role == "treat":
                if condition:
                    if g["condition"] is None:
                        g["condition"] = condition
                    elif g["condition"] != condition:
                        raise WorkflowError(
                            f"Sample table line {lineno}: treat samples of group {grp} declare "
                            f"different condition values ({g['condition']!r} vs {condition!r}); "
                            "a group's treats must share one condition"
                        )
                treat_condition_of[sid] = condition
                treat_batch_of[sid] = batch

    bad_groups = [g for g, v in groups.items() if not v["treat"]]
    if bad_groups:
        raise WorkflowError(f"These groups have no role=treat sample: {bad_groups}")
    if not samples:
        raise WorkflowError(f"Sample table {path} has no data rows")

    return samples, groups, seqtype_of, treat_condition_of, treat_batch_of


SAMPLES, GROUPS, SEQTYPE_OF, TREAT_CONDITION_OF, TREAT_BATCH_OF = \
    load_sample_table(_resolve_sample_table(config["grouplist"]))

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
        # Replicate-aware peak stage (optional; defaults to disabled)
        rep = p.get("replicate")
        if rep is not None:
            if not isinstance(rep, dict):
                errors.append(f"peak.replicate must be a mapping, got {rep!r}")
            else:
                if not isinstance(rep.get("enabled", False), bool):
                    errors.append(f"peak.replicate.enabled must be true/false, got {rep.get('enabled')!r}")
                for key in ("qvalue", "idr_threshold"):
                    try:
                        if not 0 < float(rep[key]) <= 1:
                            errors.append(f"peak.replicate.{key} must be in (0, 1], got {rep[key]!r}")
                    except (TypeError, ValueError):
                        errors.append(f"peak.replicate.{key} must be numeric, got {rep[key]!r}")
                if rep.get("idr_rank", "p.value") not in ("p.value", "signal.value"):
                    errors.append(
                        "peak.replicate.idr_rank must be p.value or signal.value, got "
                        f"{rep.get('idr_rank')!r}")
                mr = rep.get("consensus_min_replicates", 2)
                if isinstance(mr, bool) or not isinstance(mr, int) or mr < 2:
                    errors.append(
                        f"peak.replicate.consensus_min_replicates must be an integer >= 2, got {mr!r}")
                if rep.get("frip_on", "pooled") not in ("pooled", "consensus"):
                    errors.append(
                        "peak.replicate.frip_on must be pooled or consensus, got "
                        f"{rep.get('frip_on')!r}")
                if (rep.get("frip_on", "pooled") == "consensus"
                        and not rep.get("enabled", False)):
                    errors.append(
                        "peak.replicate.frip_on=consensus requires peak.replicate.enabled=true")
        # Alternative pooled peak caller (v0.6): macs2 (default) | seacr. The
        # SEACR knobs are validated whenever present; caller=seacr together
        # with the replicate stage is an aggregated error (the replicate/IDR
        # machinery is MACS2-based, see workflow/rules/callpeak_replicate.smk).
        caller = p.get("caller", "macs2")
        if caller not in ("macs2", "seacr"):
            errors.append(f"peak.caller must be macs2 or seacr, got {caller!r}")
        sc = p.get("seacr")
        if sc is not None:
            if not isinstance(sc, dict):
                errors.append(f"peak.seacr must be a mapping, got {sc!r}")
            else:
                if sc.get("mode", "stringent") not in ("stringent", "relaxed"):
                    errors.append("peak.seacr.mode must be stringent or relaxed, got "
                                  f"{sc.get('mode')!r}")
                if sc.get("normalize", "norm") not in ("norm", "non"):
                    errors.append("peak.seacr.normalize must be norm or non, got "
                                  f"{sc.get('normalize')!r}")
                try:
                    if not 0 < float(sc.get("fdr_threshold", 0.01)) < 1:
                        errors.append("peak.seacr.fdr_threshold must be in (0, 1), got "
                                      f"{sc.get('fdr_threshold')!r}")
                except (TypeError, ValueError):
                    errors.append("peak.seacr.fdr_threshold must be numeric, got "
                                  f"{sc.get('fdr_threshold')!r}")
        _rep_for_caller = p.get("replicate")
        if (caller == "seacr" and isinstance(_rep_for_caller, dict)
                and _rep_for_caller.get("enabled", False)):
            errors.append(
                "peak.caller=seacr and peak.replicate.enabled=true are mutually "
                "exclusive: the replicate/IDR stage is MACS2-based; enable one "
                "of the two switches, not both")
    bl = cfg.get("blacklist", "")
    if not isinstance(bl, str):
        errors.append(f"blacklist must be a string path (or empty to disable), got {bl!r}")
    # Signal tracks (v0.5 Phase 4)
    bw = cfg.get("bigwig")
    if bw is not None:
        if not isinstance(bw, dict):
            errors.append(f"bigwig must be a mapping, got {bw!r}")
        else:
            if not isinstance(bw.get("per_sample", False), bool):
                errors.append(f"bigwig.per_sample must be true/false, got {bw.get('per_sample')!r}")
            if bw.get("normalize", "RPGC") not in ("RPGC", "CPM"):
                errors.append(
                    f"bigwig.normalize must be RPGC or CPM, got {bw.get('normalize')!r}")
            bsz = bw.get("bin", 25)
            if isinstance(bsz, bool) or not isinstance(bsz, int) or bsz < 1:
                errors.append(f"bigwig.bin must be an integer >= 1, got {bsz!r}")
    if isinstance(cfg.get("peak"), dict) and cfg["peak"].get("bigwig_measure", "FE") not in ("FE", "logFE"):
        errors.append(
            "peak.bigwig_measure must be FE or logFE, got "
            f"{cfg['peak'].get('bigwig_measure')!r}")
    # HOMER motif enrichment (v0.5 Phase 5)
    mo = cfg.get("motif")
    if mo is not None:
        if not isinstance(mo, dict):
            errors.append(f"motif must be a mapping, got {mo!r}")
        else:
            for key in ("homer_genome", "size", "background", "extra"):
                if key in mo and not isinstance(mo[key], str):
                    errors.append(f"motif.{key} must be a string, got {mo[key]!r}")
            if not isinstance(mo.get("enabled", False), bool):
                errors.append(f"motif.enabled must be true/false, got {mo.get('enabled')!r}")
            if mo.get("enabled", False) and not str(mo.get("homer_genome", "")).strip():
                errors.append("motif.homer_genome is required when motif.enabled is true "
                              "(a HOMER genome tag such as hg38, or custom:/path/to/genome)")
    # DiffBind differential binding (v0.5 Phase 3)
    db = cfg.get("diffbind")
    if db is not None:
        if not isinstance(db, dict):
            errors.append(f"diffbind must be a mapping, got {db!r}")
        else:
            if not isinstance(db.get("enabled", False), bool):
                errors.append(f"diffbind.enabled must be true/false, got {db.get('enabled')!r}")
            if db.get("analysis", "DESeq2") not in ("DESeq2", "edgeR"):
                errors.append(
                    f"diffbind.analysis must be DESeq2 or edgeR, got {db.get('analysis')!r}")
            sf = db.get("summit_flank", 250)
            if isinstance(sf, bool) or not isinstance(sf, int) or sf < 0:
                errors.append(f"diffbind.summit_flank must be an integer >= 0, got {sf!r}")
            for key in ("use_controls", "batch_correction"):
                if key in db and not isinstance(db[key], bool):
                    errors.append(f"diffbind.{key} must be true/false, got {db[key]!r}")
            try:
                if not 0 < float(db.get("fdr", 0.05)) <= 1:
                    errors.append(f"diffbind.fdr must be in (0, 1], got {db.get('fdr')!r}")
            except (TypeError, ValueError):
                errors.append(f"diffbind.fdr must be numeric, got {db.get('fdr')!r}")
            try:
                if float(db.get("foldchange", 1.0)) < 1:
                    errors.append(f"diffbind.foldchange must be >= 1, got {db.get('foldchange')!r}")
            except (TypeError, ValueError):
                errors.append(f"diffbind.foldchange must be numeric, got {db.get('foldchange')!r}")
    if isinstance(cfg.get("qc"), dict):
        for key in ("tss", "organelle"):
            if key in cfg["qc"] and not isinstance(cfg["qc"][key], bool):
                errors.append(f"qc.{key} must be true/false, got {cfg['qc'][key]!r}")
        pats = cfg["qc"].get("organelle_patterns")
        if pats is not None and (not isinstance(pats, list)
                                 or not all(isinstance(x, str) and x for x in pats)):
            errors.append(f"qc.organelle_patterns must be a list of non-empty strings, got {pats!r}")
        # QC gate summary (v0.6): qc.gates block, default off; thresholds are
        # informational and validated as numbers in [0, 1] except the NSC/RSC/
        # TSS floors (any positive number).
        gates = cfg["qc"].get("gates")
        if gates is not None:
            if not isinstance(gates, dict):
                errors.append(f"qc.gates must be a mapping, got {gates!r}")
            else:
                if not isinstance(gates.get("enabled", False), bool):
                    errors.append(
                        f"qc.gates.enabled must be true/false, got {gates.get('enabled')!r}")
                thr = gates.get("thresholds")
                if thr is not None:
                    if not isinstance(thr, dict):
                        errors.append(f"qc.gates.thresholds must be a mapping, got {thr!r}")
                    else:
                        for key in ("mapping_rate_min", "dup_rate_max",
                                    "frip_min", "organelle_max"):
                            if key in thr:
                                try:
                                    if not 0 <= float(thr[key]) <= 1:
                                        errors.append(
                                            f"qc.gates.thresholds.{key} must be in [0, 1], "
                                            f"got {thr[key]!r}")
                                except (TypeError, ValueError):
                                    errors.append(
                                        f"qc.gates.thresholds.{key} must be numeric, "
                                        f"got {thr[key]!r}")
                        for key in ("nsc_min", "rsc_min", "tss_min"):
                            if key in thr:
                                try:
                                    if not float(thr[key]) > 0:
                                        errors.append(
                                            f"qc.gates.thresholds.{key} must be a positive "
                                            f"number, got {thr[key]!r}")
                                except (TypeError, ValueError):
                                    errors.append(
                                        f"qc.gates.thresholds.{key} must be numeric, "
                                        f"got {thr[key]!r}")
    # Spike-in normalization (v0.6): spike_in block, default off. The FASTA
    # is a required non-empty string when enabled; a missing file (or a
    # "/path/to/"-style placeholder, which simply does not exist) only warns
    # like every other reference.
    sp = cfg.get("spike_in")
    if sp is not None:
        if not isinstance(sp, dict):
            errors.append(f"spike_in must be a mapping, got {sp!r}")
        else:
            if not isinstance(sp.get("enabled", False), bool):
                errors.append(
                    f"spike_in.enabled must be true/false, got {sp.get('enabled')!r}")
            if not isinstance(sp.get("scale_bigwigs", False), bool):
                errors.append(
                    f"spike_in.scale_bigwigs must be true/false, got {sp.get('scale_bigwigs')!r}")
            fasta = sp.get("fasta", "")
            if not isinstance(fasta, str):
                errors.append(f"spike_in.fasta must be a string path, got {fasta!r}")
            elif sp.get("enabled", False) and not fasta.strip():
                errors.append("spike_in.fasta is required when spike_in.enabled is true")
            name = sp.get("name", "spike")
            if not isinstance(name, str) or not name.strip():
                errors.append(f"spike_in.name must be a non-empty string, got {name!r}")
    # TOBIAS footprinting (v0.6): footprint block, default off. motifs is a
    # required non-empty string when enabled (the TOBIAS footprinting recipes
    # scan a fixed motif set); a "/path/to/"-style placeholder, which simply
    # does not exist, warns only like every other reference file.
    fp = cfg.get("footprint")
    if fp is not None:
        if not isinstance(fp, dict):
            errors.append(f"footprint must be a mapping, got {fp!r}")
        else:
            if not isinstance(fp.get("enabled", False), bool):
                errors.append(
                    f"footprint.enabled must be true/false, got {fp.get('enabled')!r}")
            if not isinstance(fp.get("bindetect", True), bool):
                errors.append(
                    f"footprint.bindetect must be true/false, got {fp.get('bindetect')!r}")
            motifs = fp.get("motifs", "")
            if not isinstance(motifs, str):
                errors.append(f"footprint.motifs must be a string path, got {motifs!r}")
            elif fp.get("enabled", False) and not motifs.strip():
                errors.append(
                    "footprint.motifs is required when footprint.enabled is true "
                    "(a JASPAR/HOMER/MEME motif PFM file)")
            mp = fp.get("motif_pvalue")
            if mp not in (None, ""):
                try:
                    if not 0 < float(mp) <= 1:
                        errors.append(
                            f"footprint.motif_pvalue must be in (0, 1], got {mp!r}")
                except (TypeError, ValueError):
                    errors.append(
                        f"footprint.motif_pvalue must be numeric, got {mp!r}")
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
    bl = str(cfg.get("blacklist") or "").strip()
    if bl and not os.path.exists(bl):
        warnings.append(f"blacklist file does not exist (verify before running): {bl}")
    _sp_cfg = cfg.get("spike_in") if isinstance(cfg.get("spike_in"), dict) else {}
    if _sp_cfg.get("enabled", False):
        _sp_fa = str(_sp_cfg.get("fasta") or "").strip()
        _sp_path = _sp_fa
        if _sp_path and not os.path.isabs(_sp_path) and not os.path.exists(_sp_path):
            _sp_path = os.path.join(BASE_DIR, _sp_path)
        if _sp_path and not os.path.exists(_sp_path):
            warnings.append(
                f"spike-in FASTA does not exist (verify before running): spike_in.fasta = {_sp_fa}")
    _fp_cfg = cfg.get("footprint") if isinstance(cfg.get("footprint"), dict) else {}
    if _fp_cfg.get("enabled", False):
        _fp_motifs = str(_fp_cfg.get("motifs") or "").strip()
        _fp_path = _fp_motifs
        if _fp_path and not os.path.isabs(_fp_path) and not os.path.exists(_fp_path):
            _fp_path = os.path.join(BASE_DIR, _fp_path)
        if _fp_path and not os.path.exists(_fp_path):
            warnings.append(
                f"footprint motif file does not exist (verify before running): "
                f"footprint.motifs = {_fp_motifs}")
    for w in warnings:
        print(f"[config warning] {w}")


validate_config(config)


# ---------------------------------------------------------------------
# Pure helpers for the replicate-aware peak stage and the extended QC
# (no config/GROUPS access; unit-tested via tests/run_tests.py extraction)
# ---------------------------------------------------------------------

def _unordered_pairs(items):
    """All unordered pairs of a list, in order (the replicate pairs IDR runs on)."""
    return [(items[i], items[j])
            for i in range(len(items)) for j in range(i + 1, len(items))]


def idr_pair_slug(a, b):
    """Filesystem token for one IDR comparison. Sample names cannot contain
    '__' (table validation), so the slug splits back unambiguously."""
    return f"{a}__vs__{b}"


def parse_idr_pair_slug(slug):
    a, b = slug.split("__vs__")
    return a, b


def _group_calls_narrow(v):
    """Whether a group's peaks are narrow: atac/faire always are (their
    peak_type is fixed to none by the table validation), chip/cuttag by
    peak_type. Decides IDR (narrow) vs overlap consensus (broad)."""
    return v["seqtype"] in ("atac", "faire") or v["peak_type"] == "narrow"


def _is_organelle_contig(name, patterns):
    """Classify a reference contig as organelle (chloroplast/mitochondrion).
    Short patterns (<= 3 chars, e.g. pt/mt/chrc/chrm) must equal the contig
    name case-insensitively — substring matching would swallow names like
    human ALT contigs; longer patterns (chloroplast, mitochondr...) match as
    substrings."""
    lowered = name.lower()
    for pattern in patterns:
        p = str(pattern).lower()
        if len(p) <= 3:
            if lowered == p:
                return True
        elif p in lowered:
            return True
    return False


# ---------------------------------------------------------------------
# Replicate-aware peak analysis and extended QC (v0.5). Every switch
# defaults to off, so existing project configs and the default dry-run
# baseline keep working unchanged.
# ---------------------------------------------------------------------
REPLICATE_DEFAULTS = {
    "enabled": False,
    "qvalue": 0.01,              # relaxed per-replicate narrow cutoff feeding IDR (ENCODE)
    "idr_threshold": 0.05,
    "idr_rank": "p.value",       # narrowPeak rank column: p.value | signal.value
    "consensus_min_replicates": 2,
    "frip_on": "pooled",         # peak set FRiP is computed against
}
REPLICATE = dict(REPLICATE_DEFAULTS)
REPLICATE.update(((config.get("peak") or {}).get("replicate") or {}))

BLACKLIST = str(config.get("blacklist") or "").strip()
QC_TSS = bool((config.get("qc") or {}).get("tss", False))
QC_ORGANELLE = bool((config.get("qc") or {}).get("organelle", False))
ORGANELLE_PATTERNS = list((config.get("qc") or {}).get("organelle_patterns")
                          or ["chrc", "chrm", "pt", "mt", "pltd",
                              "chloroplast", "mitochondr", "plastid"])

# QC gate summary (v0.6, default off): one PASS/WARN/FAIL row per sample
# aggregating the existing QC metrics (mapping rate, duplication, FRiP,
# NSC/RSC, TSS enrichment, organelle fraction). Strictly informational —
# the pipeline never hard-fails on a gate; aggregation lives in
# workflow/scripts/gates_summary.py.
GATES_DEFAULTS = {
    "enabled": False,
    "thresholds": {
        "mapping_rate_min": 0.70,   # flagstat mapped/total of the analysis BAM
        "dup_rate_max": 0.50,       # picard PERCENT_DUPLICATION (NA when dedup is off)
        "frip_min": 0.01,           # ENCODE TF floor; loosen for histone marks
        "nsc_min": 1.05,
        "rsc_min": 0.8,
        "tss_min": 6.0,
        "organelle_max": 0.20,
    },
}
_gates_cfg = (config.get("qc") or {}).get("gates") or {}
_gates_user_thr = _gates_cfg.get("thresholds") or {}
_unknown_gates_thr = [k for k in _gates_user_thr if k not in GATES_DEFAULTS["thresholds"]]
if _unknown_gates_thr:
    print(f"[config warning] unknown qc.gates.thresholds keys are ignored: {_unknown_gates_thr}")
_gates_thr = dict(GATES_DEFAULTS["thresholds"])
_gates_thr.update({k: v for k, v in _gates_user_thr.items()
                   if k in GATES_DEFAULTS["thresholds"]})
GATES = {"enabled": bool(_gates_cfg.get("enabled", False)), "thresholds": _gates_thr}

# Spike-in normalization (v0.6, default off): second-pass alignment of the
# unmapped read pairs against a spike-in genome, per-sample scale factors
# (1e6 / spike-in mapped reads) and a QC summary; scale_bigwigs optionally
# multiplies the bigWig signal tracks by the per-sample factors (group tracks
# use the mean over their treat samples).
SPIKE_IN_DEFAULTS = {"enabled": False, "fasta": "", "name": "spike",
                     "scale_bigwigs": False}
SPIKE_IN = dict(SPIKE_IN_DEFAULTS)
SPIKE_IN.update((config.get("spike_in") or {}))

# Alternative pooled peak caller (v0.6, default off via caller=macs2): with
# peak.caller=seacr the pooled MACS2 narrow/broad/atac group calls
# (callpeak.smk) are replaced by SEACR (workflow/rules/seacr.smk) over
# raw-depth bedGraph coverage, with SEACR's 6-column result converted back to
# the standard narrowPeak/broadPeak contract at the same {group}_peaks.* paths,
# so every downstream consumer (FRiP, annotation, blacklist, motif, DiffBind)
# is unchanged. The caller=seacr + peak.replicate.enabled combination is
# rejected by validate_config (the replicate/IDR machinery is MACS2-based).
PEAK_CALLER_DEFAULTS = {"caller": "macs2",
                        "seacr": {"mode": "stringent", "normalize": "norm",
                                  "fdr_threshold": 0.01}}
_peak_caller_cfg = config.get("peak") or {}
_seacr_user_cfg = _peak_caller_cfg.get("seacr") or {}
_unknown_seacr_keys = [k for k in _seacr_user_cfg if k not in PEAK_CALLER_DEFAULTS["seacr"]]
if _unknown_seacr_keys:
    print(f"[config warning] unknown peak.seacr keys are ignored: {_unknown_seacr_keys}")
_seacr_cfg = dict(PEAK_CALLER_DEFAULTS["seacr"])
_seacr_cfg.update({k: v for k, v in _seacr_user_cfg.items()
                   if k in PEAK_CALLER_DEFAULTS["seacr"]})
PEAK_CALLER = {"caller": _peak_caller_cfg.get("caller", "macs2"), "seacr": _seacr_cfg}

# TOBIAS footprinting (v0.6, default off): ATACorrect Tn5-bias correction over
# the pooled analysis BAM of every atac/faire group, ScoreBigwig footprint
# scores against a fixed motif PFM set, and (optional) BINDetect TF binding
# detection. Footprinting models cut-site/Tn5 bias, so it is meaningful only
# for the open-chromatin assays — chip/cuttag enrichment groups never get
# footprint jobs (see FOOTPRINT_GROUPS below). motif_pvalue is passed to
# BINDetect --motif-pvalue only when set (empty/None keeps the tool default).
FOOTPRINT_DEFAULTS = {"enabled": False, "motifs": "", "bindetect": True,
                      "motif_pvalue": 1e-4}
FOOTPRINT = dict(FOOTPRINT_DEFAULTS)
FOOTPRINT.update((config.get("footprint") or {}))
FOOTPRINT_MOTIF_PVALUE_ARG = (
    "" if FOOTPRINT.get("motif_pvalue") in (None, "")
    else f"--motif-pvalue {FOOTPRINT['motif_pvalue']}")

# Footprinting covers the ATAC/FAIRE groups only.
FOOTPRINT_GROUPS = [g for g, v in GROUPS.items() if v["seqtype"] in ("atac", "faire")]

# samtools flagstat count-line patterns (keep in sync with the compiled
# copies in workflow/scripts/gates_summary.py and
# workflow/scripts/spikein_summary.py — the standalone scripts duplicate
# them on purpose). The "primary mapped" and "with itself and mate mapped"
# lines must not match.
_FLAGSTAT_TOTAL_RE = re.compile(r"^\s*(\d+)\s+\+\s+\d+\s+in total\b")
_FLAGSTAT_MAPPED_RE = re.compile(r"^\s*(\d+)\s+\+\s+\d+\s+mapped(?:\s|\(|$)")


def flagstat_counts(path):
    """(total, mapped) from a samtools flagstat report; (0, 0) when the file
    is unreadable. Run-time helper — call it from params lambdas only, the
    flagstat files do not exist at parse time."""
    total = mapped = 0
    try:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                m = _FLAGSTAT_TOTAL_RE.match(line)
                if m:
                    total = int(m.group(1))
                    continue
                m = _FLAGSTAT_MAPPED_RE.match(line)
                if m:
                    mapped = int(m.group(1))
    except OSError:
        return 0, 0
    return total, mapped


def spike_scale_factor(path):
    """Per-sample spike-in scale factor from one flagstat report:
    1e6 / max(spike-in mapped reads, 1)."""
    _, mapped = flagstat_counts(path)
    return 1e6 / max(mapped, 1)


def _spike_scaled_bedgraph_stage(factors, sorted_path):
    """Shell snippet multiplying a sorted bedGraph value column by a spike-in
    scale factor, in place (`factors` is averaged first: a group FE track
    uses the mean over its treat samples). Returns a leading command plus a
    newline/indent tail, so prepending it to the next command keeps the
    rendered script byte-identical when the stage is off (empty string).
    Built by concatenation: the pre-3.12 f-string parser rejects a doubled
    brace directly after a format field, and the awk program needs literal
    braces (params values are never re-formatted by snakemake)."""
    factor = sum(factors) / len(factors)
    scaled = sorted_path + ".scaled"
    awk = ("awk 'BEGIN{OFS=\"\\t\"} NF>=4{print $1, $2, $3, $4*"
           + format(factor, ".6f")
           + "; next}{print}' "
           + sorted_path + " > " + scaled
           + " && mv " + scaled + " " + sorted_path + "\n        ")
    return awk


def _spike_scale_flag_arg(path):
    """' --scaleFactor F' bamCoverage argument built from one sample flagstat."""
    return f" --scaleFactor {spike_scale_factor(path):.6f}"

# Signal tracks (v0.5 Phase 4): per-sample normalized coverage bigWigs and the
# group-track measure; both optional, defaults keep today's FE-only behavior.
BIGWIG_DEFAULTS = {"per_sample": False, "normalize": "RPGC", "bin": 25}
BIGWIG = dict(BIGWIG_DEFAULTS)
BIGWIG.update((config.get("bigwig") or {}))
PEAK_MEASURE = str((config.get("peak") or {}).get("bigwig_measure", "FE"))

# HOMER motif enrichment (v0.5 Phase 5), default off.
MOTIF_DEFAULTS = {"enabled": False, "homer_genome": "", "size": "given",
                  "background": "", "extra": ""}
MOTIF = dict(MOTIF_DEFAULTS)
MOTIF.update((config.get("motif") or {}))

# DiffBind differential binding (v0.5 Phase 3), default off. Contrasts are
# pairs of sample-table groups; each arm needs >= 2 treat replicates.
DIFFBIND_DEFAULTS = {"enabled": False, "contrasts": [], "analysis": "DESeq2",
                     "summit_flank": 250, "use_controls": False,
                     "fdr": 0.05, "foldchange": 1.0, "batch_correction": True}
DIFFBIND = dict(DIFFBIND_DEFAULTS)
DIFFBIND.update((config.get("diffbind") or {}))

# Replicate-stage group enumeration. Every group calls peaks, so every treat
# sample gets a per-replicate peak file; IDR (narrow) / overlap consensus
# (broad) are only built for groups with >= 2 treats.
NARROW_REP_GROUPS = [g for g, v in GROUPS.items() if _group_calls_narrow(v)]
BROAD_REP_GROUPS = [g for g, v in GROUPS.items() if not _group_calls_narrow(v)]
IDR_GROUPS = [g for g in NARROW_REP_GROUPS if len(GROUPS[g]["treat"]) >= 2]
BROAD_CONSENSUS_GROUPS = [g for g in BROAD_REP_GROUPS if len(GROUPS[g]["treat"]) >= 2]
IDR_PAIRS = {g: _unordered_pairs(GROUPS[g]["treat"]) for g in IDR_GROUPS}
REPLICATE_TREATS = {g: list(v["treat"]) for g, v in GROUPS.items()}

# TSS enrichment covers the treat samples of the open-chromatin assays
# (controls carry no TSS signal worth plotting).
TSS_SAMPLES = []
for _g, _v in GROUPS.items():
    if _v["seqtype"] in ("atac", "faire"):
        for _s in _v["treat"]:
            if _s not in TSS_SAMPLES:
                TSS_SAMPLES.append(_s)

if REPLICATE["enabled"]:
    _single = [g for g in GROUPS if len(GROUPS[g]["treat"]) < 2]
    if _single:
        print(f"[replicate notice] single-treat groups fall back to the pooled "
              f"peak set (no IDR/consensus exists for them): {_single}")
    _few = [g for g in IDR_GROUPS + BROAD_CONSENSUS_GROUPS
            if len(GROUPS[g]["treat"]) < REPLICATE["consensus_min_replicates"]]
    if _few:
        print(f"[config warning] groups with fewer treats than "
              f"peak.replicate.consensus_min_replicates={REPLICATE['consensus_min_replicates']} "
              f"produce an empty consensus: {_few}")
if QC_TSS and not TSS_SAMPLES:
    print("[config warning] qc.tss is enabled but the table has no atac/faire "
          "treat sample; the TSS stage is skipped")
if FOOTPRINT["enabled"] and not FOOTPRINT_GROUPS:
    print("[config warning] footprint is enabled but the table has no atac/faire "
          "group; the footprinting stage is skipped")

# Gate summary covers every sample; the FRiP filenames need each sample's
# group (first group wins for a sample listed in several groups).
GATES_SAMPLE_GROUPS = {}
for _g, _v in GROUPS.items():
    for _s in _v["treat"] + _v["control"]:
        GATES_SAMPLE_GROUPS.setdefault(_s, _g)

# DiffBind contrast enumeration (slugs reuse the '__vs__' convention; group
# names cannot contain '__' so the split is unambiguous).
_contrast_errors = []
for _pair in (DIFFBIND.get("contrasts") if isinstance(DIFFBIND.get("contrasts"), list) else None) or []:
    if not (isinstance(_pair, (list, tuple)) and len(_pair) == 2
            and all(isinstance(_x, str) and _x for _x in _pair)):
        _contrast_errors.append(f"every diffbind.contrasts entry must be a pair of group names, got {_pair!r}")
        continue
    for _g in _pair:
        if _g not in GROUPS:
            _contrast_errors.append(f"diffbind.contrasts references unknown group {_g!r}")
        elif len(GROUPS[_g]["treat"]) < 2:
            _contrast_errors.append(
                f"diffbind.contrast arm {_g!r} needs >= 2 treat replicates "
                f"(has {len(GROUPS[_g]['treat'])})")
if _contrast_errors:
    raise WorkflowError("diffbind.contrasts validation failed:\n  " + "\n  ".join(_contrast_errors))
DIFFBIND_CONTRASTS = [(str(a), str(b), idr_pair_slug(str(a), str(b)))
                      for a, b in (DIFFBIND["contrasts"] or [])]
DIFFBIND_GROUPS = []
for _a, _b, _ in DIFFBIND_CONTRASTS:
    for _g in (_a, _b):
        if _g not in DIFFBIND_GROUPS:
            DIFFBIND_GROUPS.append(_g)
if DIFFBIND["enabled"] and not DIFFBIND_CONTRASTS:
    print("[config warning] diffbind.enabled is true but diffbind.contrasts "
          "is empty; the differential stage is skipped")
if DIFFBIND["enabled"] and DIFFBIND["use_controls"]:
    _multi_ctl = [g for g in DIFFBIND_GROUPS if len(GROUPS[g]["control"]) != 1]
    if _multi_ctl:
        print(f"[config warning] diffbind.use_controls: groups without exactly "
              f"one control lose their control background: {_multi_ctl}")
if MOTIF["enabled"] and MOTIF["background"] and not os.path.exists(MOTIF["background"]):
    print(f"[config warning] motif.background does not exist (verify before "
          f"running): {MOTIF['background']}")

# ---------------------------------------------------------------------
# Output redirection (aligned with rna-seq): raw inputs (1.rawdata/) stay
# at the project working-directory root; every derived artifact lives
# under the configurable results root (default "results").
# ---------------------------------------------------------------------
RD = str(config.get("results_dir", "results")).rstrip(os.sep) + os.sep


def R(path=""):
    """Return a path under the configured results directory."""
    return f"{RD}{path}"


# R interpreter channel: run.sh resolves r.rscript from software.yaml and
# exports it as CHIP_RSCRIPT (resolved once at orchestrator parse time and
# baked into the jobscript). Rules must invoke {params.rscript} instead of a
# bare `Rscript`, which on cluster nodes would resolve to whatever R happens
# to be on PATH and bypass the configured R runtime/libraries.
RSCRIPT = os.environ.get("CHIP_RSCRIPT", "Rscript")

# idr resolves like the R channel: run.sh exports CHIP_IDR from the
# software.yaml paths: section (the classic idr tool is python2-based and
# deliberately stays OUT of the main conda environment; only needed when
# peak.replicate.enabled is true).
IDR_BIN = os.environ.get("CHIP_IDR", "idr")

# HOMER findMotifsGenome.pl resolves the same way (external distribution,
# configured via software.yaml paths: -> CHIP_HOMER_FINDMOTIFS; only needed
# when motif.enabled is true).
HOMER_BIN = os.environ.get("CHIP_HOMER_FINDMOTIFS", "findMotifsGenome.pl")

# SEACR (the alternative pooled peak caller, peak.caller=seacr only) resolves
# the same way: an external bash script outside the conda template, configured
# via software.yaml paths: -> CHIP_SEACR. The rules invoke it as
# `bash {SEACR_BIN} ...` (SEACR is distributed as a shell script, not a binary).
SEACR_BIN = os.environ.get("CHIP_SEACR", "SEACR_1.3.sh")

# TOBIAS (the footprinting stage, footprint.enabled only) resolves the same
# way: a heavyweight external tool deliberately kept OUT of the conda
# template (its own environment, e.g. conda create -n chip-tobias -c bioconda
# tobias), configured via software.yaml paths: -> CHIP_TOBIAS. The entry point
# is the env's bin/TOBIAS console script (an absolute path works standalone
# via its shebang interpreter).
TOBIAS_BIN = os.environ.get("CHIP_TOBIAS", "TOBIAS")

# Python interpreter for the workflow's own scripts (run.sh exports CHIP_PYTHON
# from the runtime resolution; python3 is the sane fallback).
PYTHON_BIN = os.environ.get("CHIP_PYTHON", "python3")


# ---------------------------------------------------------------------
# Shared query helpers (used by the rules/*.smk modules)
# ---------------------------------------------------------------------

def raw_fastq_pair(sample):
    """Resolve the raw paired-end FASTQ paths for one sample from 1.rawdata/.

    Supported naming variants, checked in this priority order (the first
    complete pair wins; identical to the rna-seq workflow's resolution
    order so mixed repository deployments behave consistently):
      {id}_1.fastq.gz + {id}_2.fastq.gz
      {id}_1.fq.gz + {id}_2.fq.gz
      {id}_R1.fastq.gz + {id}_R2.fastq.gz
      {id}_R1.fq.gz + {id}_R2.fq.gz
    Raises WorkflowError when no complete pair exists (the pipeline is
    paired-end only; single-end files alone are not accepted).
    """
    for suffix1, suffix2 in (
        ("_1.fastq.gz", "_2.fastq.gz"),
        ("_1.fq.gz", "_2.fq.gz"),
        ("_R1.fastq.gz", "_R2.fastq.gz"),
        ("_R1.fq.gz", "_R2.fq.gz"),
    ):
        fq1 = f"1.rawdata/{sample}{suffix1}"
        fq2 = f"1.rawdata/{sample}{suffix2}"
        if os.path.exists(fq1) and os.path.exists(fq2):
            return fq1, fq2
    raise WorkflowError(
        f"sample {sample!r}: no paired-end raw FASTQ pair found under 1.rawdata/. "
        f"Supported naming variants are {{id}}_1.fastq.gz + {{id}}_2.fastq.gz, "
        f"{{id}}_1.fq.gz + {{id}}_2.fq.gz, {{id}}_R1.fastq.gz + {{id}}_R2.fastq.gz, "
        f"or {{id}}_R1.fq.gz + {{id}}_R2.fq.gz (both mates required; "
        "the pipeline is paired-end only)."
    )


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


# ----- Replicate stage / blacklist path helpers -----

def replicate_peak_file(group, sample):
    """Per-replicate peak file for one treat sample (replicate stage)."""
    suffix = "narrowPeak" if _group_calls_narrow(GROUPS[group]) else "broadPeak"
    return f"{RD}4.peak/replicates/{group}/{sample}_peaks.{suffix}"


def idr_pair_file(group, a, b):
    """One pairwise IDR comparison of two treat replicates."""
    return f"{RD}4.peak/idr/{group}/{idr_pair_slug(a, b)}.narrowPeak"


def group_idr_file(group):
    """Final reproducible narrow peak set (IDR) for a >=2-treat group."""
    return f"{RD}4.peak/{group}_IDR_peaks.narrowPeak"


def group_idr_support_file(group):
    """Per-peak replicate-support BED alongside the IDR peak set."""
    return f"{RD}4.peak/{group}_IDR_support.bed"


def group_consensus_file(group):
    """Final reproducible broad peak set (overlap consensus)."""
    return f"{RD}4.peak/{group}_consensus_peaks.broadPeak"


def group_consensus_support_file(group):
    """Per-peak replicate-support BED alongside the broad consensus."""
    return f"{RD}4.peak/{group}_consensus_support.bed"


def _blacklist_variant(path):
    """Blacklist-filtered copy path for a peak file (same basename, dedicated
    directory; only produced when a blacklist is configured)."""
    return f"{RD}4.peak/blacklist_filtered/{os.path.basename(path)}"


def group_final_peak_file(group):
    """The group's reproducible peak set when the replicate stage is enabled
    (IDR for narrow, overlap consensus for broad); the pooled set otherwise.
    Single-treat groups have no consensus and always fall back to pooled."""
    if REPLICATE["enabled"]:
        if group in IDR_GROUPS:
            return group_idr_file(group)
        if group in BROAD_CONSENSUS_GROUPS:
            return group_consensus_file(group)
    return group_peak_file(group)


def frip_peak_file(group):
    """Peak set the FRiP metric is computed against (frip_on switch), in its
    blacklist-filtered variant when a blacklist is configured."""
    peak = (group_final_peak_file(group) if REPLICATE["frip_on"] == "consensus"
            else group_peak_file(group))
    return _blacklist_variant(peak) if BLACKLIST else peak


def final_peak_files():
    """Final reproducible peak set per group (pooled when the replicate stage
    is off), BEFORE any blacklist filtering. The unfiltered sources the
    blacklist stage reads from."""
    return [group_final_peak_file(g) for g in GROUPS]


def annot_peak_files():
    """Peak sets the ChIPseeker annotation covers: the final reproducible set
    per group, blacklist-filtered copies when a blacklist is configured."""
    files = final_peak_files()
    if BLACKLIST:
        files = [_blacklist_variant(f) for f in files]
    return files


def blacklist_sources():
    """Ordered {basename: source path} map for the blacklist stage (basenames
    are unique: pooled/IDR/consensus files carry distinct suffixes). Sources
    are the UNFILTERED final reproducible peak files per group, plus the
    pooled set whenever it differs — FRiP reads the pooled set under
    frip_on=pooled and must get its filtered copy too."""
    files = list(final_peak_files())
    for g in GROUPS:
        pooled = group_peak_file(g)
        if pooled not in files:
            files.append(pooled)
    return {os.path.basename(f): f for f in files}


def motif_peak_file(group):
    """Peak set the motif enrichment reads for one group: the same final
    deliverable annotation uses (blacklist-filtered when configured)."""
    return dict(zip(GROUPS, annot_peak_files()))[group]


def diffbind_sample_peaks(group, sample):
    """Per-sample peak file for the differential count matrix: the replicate
    call when the replicate stage is enabled, otherwise the pooled group set
    (degraded — DiffBind then sees identical peak sets; enabling
    peak.replicate is recommended)."""
    if REPLICATE["enabled"]:
        return replicate_peak_file(group, sample)
    return group_peak_file(group)


def _contrast_groups(slug):
    """The [groupA, groupB] pair behind one diffbind contrast slug."""
    for a, b, s in DIFFBIND_CONTRASTS:
        if s == slug:
            return [a, b]
    raise WorkflowError(f"unknown diffbind contrast slug: {slug!r}")


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
    "seacr_bedgraph_treat": {"threads": 2, "mem_mb": 8192, "runtime_min": 120},
    "seacr_bedgraph_control": {"threads": 2, "mem_mb": 8192, "runtime_min": 120},
    "seacr_callpeak": {"threads": 1, "mem_mb": 8192, "runtime_min": 180},
    "seacr_bigwig": {"threads": 1, "mem_mb": 4096, "runtime_min": 60},
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
    "callpeak_narrow_replicate": {"threads": 1, "mem_mb": 8192, "runtime_min": 180},
    "callpeak_broad_replicate": {"threads": 1, "mem_mb": 8192, "runtime_min": 180},
    "idr_pair": {"threads": 1, "mem_mb": 4096, "runtime_min": 60},
    "idr_final": {"threads": 1, "mem_mb": 2048, "runtime_min": 15},
    "broad_consensus": {"threads": 1, "mem_mb": 2048, "runtime_min": 15},
    "replicate_summary": {"threads": 1, "mem_mb": 1024, "runtime_min": 10},
    "blacklist_filter": {"threads": 1, "mem_mb": 2048, "runtime_min": 15},
    "blacklist_summary": {"threads": 1, "mem_mb": 1024, "runtime_min": 10},
    "tss_bed": {"threads": 1, "mem_mb": 1024, "runtime_min": 10},
    "tss_matrix": {"threads": 2, "mem_mb": 8192, "runtime_min": 120},
    "tss_summary": {"threads": 1, "mem_mb": 1024, "runtime_min": 10},
    "organelle_idxstats": {"threads": 1, "mem_mb": 2048, "runtime_min": 30},
    "organelle_summary": {"threads": 1, "mem_mb": 1024, "runtime_min": 10},
    "gates_flagstat": {"threads": 1, "mem_mb": 2048, "runtime_min": 15},
    "qc_gates": {"threads": 1, "mem_mb": 1024, "runtime_min": 10},
    "spike_bowtie2_index": {"threads": 4, "mem_mb": 2048, "runtime_min": 30},
    "spike_align": {"threads": 8, "mem_mb": 16384, "runtime_min": 240},
    "spike_summary": {"threads": 1, "mem_mb": 1024, "runtime_min": 10},
    "bigwig_sample": {"threads": 2, "mem_mb": 8192, "runtime_min": 60},
    "motif_enrichment": {"threads": 2, "mem_mb": 8192, "runtime_min": 720},
    "diffbind_sheet": {"threads": 1, "mem_mb": 1024, "runtime_min": 10},
    "diffbind_report": {"threads": 1, "mem_mb": 16384, "runtime_min": 240},
    "footprint_ataccorrect": {"threads": 2, "mem_mb": 8192, "runtime_min": 240},
    "footprint_scorebigwig": {"threads": 1, "mem_mb": 4096, "runtime_min": 60},
    "footprint_bindetect": {"threads": 2, "mem_mb": 8192, "runtime_min": 240},
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
# Per-sample coverage tracks are produced by bigwig_sample, which lives in the
# MACS2 callpeak module (callpeak.smk); under peak.caller=seacr that module is
# not included, so the per-sample targets drop out of the DAG (the group track
# keeps coming from the seacr module's raw-depth route).
BW_TARGETS = [f"{RD}4.peak/{g}_FE.bw" for g in GROUPS]
if BIGWIG["per_sample"]:
    if PEAK_CALLER["caller"] == "macs2":
        BW_TARGETS += [f"{RD}4.peak/samples/{s}.bw" for s in SAMPLES]
    else:
        print("[config warning] bigwig.per_sample rides the MACS2 callpeak module "
              "and has no rule under peak.caller=seacr; the per-sample tracks are "
              "skipped (the group tracks still come from the seacr module)")

MOTIF_TARGETS = []
if MOTIF["enabled"]:
    MOTIF_TARGETS += [R(f"6.motif/{g}") for g in GROUPS]

DIFFBIND_TARGETS = []
if DIFFBIND["enabled"]:
    DIFFBIND_TARGETS += [R(f"6.diffbind/{slug}/samplesheet.tsv")
                         for _a, _b, slug in DIFFBIND_CONTRASTS]
    DIFFBIND_TARGETS += [R(f"6.diffbind/{slug}/DB_results.tsv")
                         for _a, _b, slug in DIFFBIND_CONTRASTS]

# TOBIAS footprinting targets (only aggregated when the stage is enabled; the
# deliverables are the corrected bigWig file and one directory per scoring /
# binding stage under the group's stage root — see workflow/rules/footprint.smk).
FOOTPRINT_TARGETS = []
if FOOTPRINT["enabled"]:
    FOOTPRINT_TARGETS += [R(f"7.footprint/{g}/ataccorrect/{g}_corrected.bw")
                          for g in FOOTPRINT_GROUPS]
    FOOTPRINT_TARGETS += [R(f"7.footprint/{g}/scorebigwig") for g in FOOTPRINT_GROUPS]
    if FOOTPRINT["bindetect"]:
        FOOTPRINT_TARGETS += [R(f"7.footprint/{g}/bindetect") for g in FOOTPRINT_GROUPS]

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

# Replicate stage targets (only aggregated when the stage is enabled).
REPLICATE_TARGETS = []
if REPLICATE["enabled"]:
    REPLICATE_TARGETS += [
        R("5.QC/replicate_peaks/Replicate_summary.tsv"),
        R("5.QC/replicate_peaks/Replicate_summary_mqc.tsv"),
    ]
    REPLICATE_TARGETS += [replicate_peak_file(g, s)
                          for g, treats in REPLICATE_TREATS.items() for s in treats]
    REPLICATE_TARGETS += [idr_pair_file(g, a, b)
                          for g in IDR_GROUPS for a, b in IDR_PAIRS[g]]
    REPLICATE_TARGETS += [group_idr_file(g) for g in IDR_GROUPS]
    REPLICATE_TARGETS += [group_idr_support_file(g) for g in IDR_GROUPS]
    REPLICATE_TARGETS += [group_consensus_file(g) for g in BROAD_CONSENSUS_GROUPS]
    REPLICATE_TARGETS += [group_consensus_support_file(g) for g in BROAD_CONSENSUS_GROUPS]

# Blacklist-filtered copies of the final peak sets (empty when no blacklist).
BLACKLIST_TARGETS = []
if BLACKLIST:
    BLACKLIST_TARGETS += [R("5.QC/blacklist/blacklist_summary.tsv")]
    BLACKLIST_TARGETS += [R(f"4.peak/blacklist_filtered/{b}")
                          for b in blacklist_sources()]

# TSS enrichment / organelle fraction targets.
TSS_TARGETS = []
if QC_TSS and TSS_SAMPLES:
    TSS_TARGETS += [R("5.QC/tss/tss.bed")]
    TSS_TARGETS += [R(f"5.QC/tss/{s}_matrix.gz") for s in TSS_SAMPLES]
    TSS_TARGETS += [R(f"5.QC/tss/{s}_tss_profile.png") for s in TSS_SAMPLES]
    TSS_TARGETS += [R(f"5.QC/tss/{s}_TSSE.txt") for s in TSS_SAMPLES]
    TSS_TARGETS += [R("5.QC/tss/TSSE_summary.tsv"), R("5.QC/tss/TSSE_summary_mqc.tsv")]

ORGANELLE_TARGETS = []
if QC_ORGANELLE:
    ORGANELLE_TARGETS += [R(f"5.QC/organelle/{s}_idxstats.tsv") for s in SAMPLES]
    ORGANELLE_TARGETS += [R("5.QC/organelle/Organelle_summary.tsv"),
                          R("5.QC/organelle/Organelle_summary_mqc.tsv")]

# QC gate summary targets (only aggregated when the stage is enabled).
GATES_TARGETS = []
if GATES["enabled"]:
    GATES_TARGETS += [R(f"5.QC/gates/{s}_flagstat.txt") for s in SAMPLES]
    GATES_TARGETS += [R("5.QC/gates/gate_summary.tsv"),
                      R("5.QC/gates/gate_summary_mqc.tsv")]

# Spike-in normalization targets (only aggregated when the stage is enabled;
# the spike-in BAMs are produced transitively by the summary's flagstat
# inputs, so only the QC deliverables are aggregated here).
SPIKEIN_TARGETS = []
if SPIKE_IN["enabled"]:
    SPIKEIN_TARGETS += [R(f"5.QC/spike_in/{s}_idxstats.txt") for s in SAMPLES]
    SPIKEIN_TARGETS += [R("5.QC/spike_in/Spikein_summary.tsv"),
                        R("5.QC/spike_in/Spikein_summary_mqc.tsv")]

# Software-version record: generated by the software_versions rule in
# meta.smk (no input dependency, scheduled freely within the DAG); always
# part of rule all.
VERSION_TARGETS = [R("5.QC/software_versions.yaml")]
