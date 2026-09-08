#!/usr/bin/env python3
"""Unit tests for the bs-seq pure-Python workflow logic (docs/TODO.md item 3).

The functions under test live in workflow/rules/common.smk, which imports
snakemake.exceptions.WorkflowError at module level (plain python cannot).
The tests therefore exec-extract the REAL function blocks from common.smk
(single source of truth; same pattern as chip_cuttag_atac_faire's
tests/run_tests.py) with a signature-compatible WorkflowError stand-in,
and sanity-check the two thin shared-layer wrappers in workflow/scripts
through their own __file__-relative shared-path resolution.

Since v0.2 the sample table accepts optional group/batch design columns
(consumed by the differential-methylation stage) and the dmr config section
plus its design validation are pinned here too.

Run: python -m pytest tests -q   (no snakemake involvement anywhere;
the suite is independent of the working directory).
"""
import csv
import importlib.util
import os
import re
import sys
from copy import deepcopy
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
COMMON_SMK = REPO / "workflow" / "rules" / "common.smk"
SAMPLES_CSV = REPO / "config" / "samples.csv"


class WorkflowError(Exception):
    """Stand-in for snakemake.exceptions.WorkflowError (single str arg)."""


# ---------------------------------------------------------------------
# Extraction of the real function blocks from common.smk.
# Module-level statements (SAMPLES = ..., apply_species_presets(...),
# validate_config(config), TARGETS aggregation) reference globals such as
# config/_SPECIES_PRESET and are excluded by the end markers below.
# ---------------------------------------------------------------------
def _segment(src, start_marker, end_markers):
    """Return src[start:] cut at the earliest end marker found."""
    start = src.index(start_marker)
    end = len(src)
    for marker in end_markers:
        pos = src.find(marker, start + len(start_marker))
        if pos != -1:
            end = min(end, pos)
    return src[start:end]


_COMMON_SRC = COMMON_SMK.read_text(encoding="utf-8")

# _NAME_RE + the accepted-header table + _resolve_sample_table +
# load_sample_table (+ the pure apply_species_presets def); the module-level
# merge call and the SAMPLES statement are excluded.
_BLOCK_TABLE = _segment(
    _COMMON_SRC,
    "_NAME_RE = re.compile",
    ["\napply_species_presets(config", "\nSAMPLES = load_sample_table"],
)
# validate_config + validate_dmr_design defs only (the module-level
# validate_config(config) call is excluded).
_BLOCK_VALIDATE = _segment(
    _COMMON_SRC,
    "def validate_config(",
    ["\nvalidate_config(config)"],
)
# apply_species_presets def only (module-level merge call excluded).
_BLOCK_PRESETS = _segment(
    _COMMON_SRC,
    "def apply_species_presets(",
    ["\napply_species_presets(config"],
)


def _exec_block(block):
    ns = {
        "csv": csv,
        "os": os,
        "re": re,
        "WorkflowError": WorkflowError,
        "BASE_DIR": str(REPO),
    }
    exec(compile(block, "common.smk(extracted)", "exec"), ns)
    return ns


_NS_TABLE = _exec_block(_BLOCK_TABLE)
load_sample_table = _NS_TABLE["load_sample_table"]
_resolve_sample_table = _NS_TABLE["_resolve_sample_table"]
apply_species_presets = _exec_block(_BLOCK_PRESETS)["apply_species_presets"]
_NS_VALIDATE = _exec_block(_BLOCK_VALIDATE)
validate_config = _NS_VALIDATE["validate_config"]
validate_dmr_design = _NS_VALIDATE["validate_dmr_design"]


def _write_table(path, rows):
    """rows[0] is the header; writes a UTF-8 CSV and returns its path."""
    with open(path, "w", newline="", encoding="utf-8") as fh:
        csv.writer(fh).writerows(rows)
    return str(path)


HEADER = ["sample_id"]


# ---------------------------------------------------------------------
# load_sample_table
# ---------------------------------------------------------------------
def test_real_samples_csv_loads_in_order():
    assert load_sample_table(str(SAMPLES_CSV)) == ["s1", "s2"]


def test_synthetic_table_keeps_order_and_accepts_legal_characters(tmp_path):
    path = _write_table(
        tmp_path / "samples.csv", [HEADER, ["z-1"], ["a.2"], ["B_3"], ["9lives"]]
    )
    assert load_sample_table(path) == ["z-1", "a.2", "B_3", "9lives"]


def test_bad_header_rejected(tmp_path):
    variants = (
        [["sample"], ["s1"]],
        [HEADER + ["extra"], ["s1", "x"]],
        [HEADER + ["batch"], ["s1", "b1"]],  # batch without group
    )
    for rows in variants:
        path = _write_table(tmp_path / "samples.csv", rows)
        with pytest.raises(WorkflowError) as excinfo:
            load_sample_table(path)
        msg = str(excinfo.value)
        assert "header must be 'sample_id', 'sample_id,group' or 'sample_id,group,batch'" in msg
        assert "'sample_id'" in msg


def test_empty_sample_id_rejected(tmp_path):
    path = _write_table(tmp_path / "samples.csv", [HEADER, ["s1"], ["   "]])
    with pytest.raises(WorkflowError) as excinfo:
        load_sample_table(path)
    msg = str(excinfo.value)
    assert "line 3" in msg and "must not be empty" in msg


@pytest.mark.parametrize("bad", ["a b", "-lead", "a,b", "s/pam"])
def test_illegal_characters_rejected(tmp_path, bad):
    path = _write_table(tmp_path / "samples.csv", [HEADER, ["s1"], [bad]])
    with pytest.raises(WorkflowError) as excinfo:
        load_sample_table(path)
    assert "illegal characters" in str(excinfo.value)


def test_double_underscore_rejected(tmp_path):
    path = _write_table(tmp_path / "samples.csv", [HEADER, ["s1"], ["a__b"]])
    with pytest.raises(WorkflowError) as excinfo:
        load_sample_table(path)
    assert "illegal characters" in str(excinfo.value)


def test_duplicate_sample_id_rejected(tmp_path):
    path = _write_table(tmp_path / "samples.csv", [HEADER, ["s1"], ["s2"], ["s1"]])
    with pytest.raises(WorkflowError) as excinfo:
        load_sample_table(path)
    msg = str(excinfo.value)
    assert "duplicate sample_id" in msg and "'s1'" in msg


def test_header_only_table_rejected(tmp_path):
    path = _write_table(tmp_path / "samples.csv", [HEADER])
    with pytest.raises(WorkflowError) as excinfo:
        load_sample_table(path)
    assert "has no data rows" in str(excinfo.value)


# ---------------------------------------------------------------------
# Sample-table design columns (group/batch, v0.2 differential methylation)
# ---------------------------------------------------------------------
def _loaded_dicts(path):
    """Call the real loader and return the module-level design dicts it
    fills (same parse as the returned id list)."""
    ids = load_sample_table(path)
    return ids, dict(_NS_TABLE["SAMPLE_GROUPS"]), dict(_NS_TABLE["SAMPLE_BATCH"])


def test_single_column_fills_design_dicts_with_none(tmp_path):
    path = _write_table(tmp_path / "samples.csv", [HEADER, ["s1"], ["s2"]])
    ids, groups, batches = _loaded_dicts(path)
    assert ids == ["s1", "s2"]
    assert groups == {"s1": None, "s2": None}
    assert batches == {"s1": None, "s2": None}


def test_group_column_parses_into_groups_dict(tmp_path):
    path = _write_table(
        tmp_path / "samples.csv",
        [HEADER + ["group"], ["s1", "control"], ["s2", "control"], ["t1", "heat"]],
    )
    ids, groups, batches = _loaded_dicts(path)
    assert ids == ["s1", "s2", "t1"]
    assert groups == {"s1": "control", "s2": "control", "t1": "heat"}
    assert batches == {"s1": None, "s2": None, "t1": None}


def test_group_and_batch_columns_parse(tmp_path):
    path = _write_table(
        tmp_path / "samples.csv",
        [
            HEADER + ["group", "batch"],
            ["s1", "control", "b1"],
            ["s2", "control", "b2"],
            ["t1", "heat", "b1"],
        ],
    )
    ids, groups, batches = _loaded_dicts(path)
    assert ids == ["s1", "s2", "t1"]
    assert groups == {"s1": "control", "s2": "control", "t1": "heat"}
    assert batches == {"s1": "b1", "s2": "b2", "t1": "b1"}


@pytest.mark.parametrize("bad", ["a b", "-lead", "a__b", "x/y"])
def test_illegal_group_value_rejected(tmp_path, bad):
    path = _write_table(tmp_path / "samples.csv", [HEADER + ["group"], ["s1", "control"], ["s2", bad]])
    with pytest.raises(WorkflowError) as excinfo:
        load_sample_table(path)
    msg = str(excinfo.value)
    assert "group" in msg and "illegal characters" in msg


def test_illegal_batch_value_rejected(tmp_path):
    path = _write_table(
        tmp_path / "samples.csv",
        [HEADER + ["group", "batch"], ["s1", "control", "b1"], ["s2", "control", "a__b"]],
    )
    with pytest.raises(WorkflowError) as excinfo:
        load_sample_table(path)
    msg = str(excinfo.value)
    assert "batch" in msg and "illegal characters" in msg


def test_empty_group_value_rejected(tmp_path):
    path = _write_table(tmp_path / "samples.csv", [HEADER + ["group"], ["s1", "control"], ["s2", ""]])
    with pytest.raises(WorkflowError) as excinfo:
        load_sample_table(path)
    msg = str(excinfo.value)
    assert "line 3" in msg and "group must not be empty" in msg and "'s2'" in msg


def test_empty_batch_value_rejected(tmp_path):
    path = _write_table(
        tmp_path / "samples.csv",
        [HEADER + ["group", "batch"], ["s1", "control", ""]],
    )
    with pytest.raises(WorkflowError) as excinfo:
        load_sample_table(path)
    msg = str(excinfo.value)
    assert "line 2" in msg and "batch must not be empty" in msg


def test_short_row_missing_batch_rejected(tmp_path):
    # A row with fewer fields than the header reads as an empty trailing cell.
    path = _write_table(tmp_path / "samples.csv", [HEADER + ["group", "batch"], ["s1", "control"]])
    with pytest.raises(WorkflowError) as excinfo:
        load_sample_table(path)
    assert "batch must not be empty" in str(excinfo.value)


# ---------------------------------------------------------------------
# _resolve_sample_table (absolute > cwd-relative > BASE_DIR-relative)
# ---------------------------------------------------------------------
def test_resolve_absolute_path_unchanged(tmp_path):
    absolute = str(tmp_path / "samples.csv")
    assert _resolve_sample_table(absolute) == absolute


def test_resolve_existing_relative_path_uses_cwd(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _write_table(tmp_path / "samples.csv", [HEADER, ["s1"]])
    resolved = _resolve_sample_table("samples.csv")
    assert os.path.isabs(resolved)
    assert Path(resolved) == tmp_path / "samples.csv"


def test_resolve_missing_relative_path_falls_back_to_base_dir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)  # ensure the table is not cwd-visible
    assert _resolve_sample_table(os.path.join("config", "samples.csv")) == str(
        SAMPLES_CSV
    )


# ---------------------------------------------------------------------
# validate_config
# ---------------------------------------------------------------------
GOOD_GENOME = "/reference/Oryza_sativa.IRGSP-1.0.dna.genome.fa"

DMR_OFF = {
    "enabled": False,
    "control_group": "control",
    "qvalue": 0.01,
    "min_diff": 25,
    "tile_len": 1000,
    "tile_step": 100,
    "min_cpg": 3,
    "batch_correction": "F",
    "min_cov": 10,
    "max_cov": 500,
}

GOOD_CFG = {
    "SampleListFile": "config/samples.csv",
    "results_dir": "results",
    "threads": 12,
    "species": "osa",
    "genome": GOOD_GENOME,
    "trim": {
        "enabled": True,
        "quality": 20,
        "min_len": 20,
        "adapter": "AGATCGGAAGAGC",
        "stringency": 3,
        "error_rate": 0.1,
        "extra": "",
    },
    "bismark": {"align_extra": "--phred33-quals"},
    "methylation_extractor": {"cx_report": False, "merge_cpg": True, "buffer_frac": 4},
    "dmr": dict(DMR_OFF),
}


def _cfg_copy():
    return deepcopy(GOOD_CFG)


def test_valid_config_passes_silently(capsys):
    validate_config(_cfg_copy())
    assert capsys.readouterr().out == ""


def test_placeholder_genome_passes_with_warning(capsys):
    cfg = _cfg_copy()
    cfg["genome"] = "/path/to/reference/osa/genome.fa"
    validate_config(cfg)  # must NOT raise
    out = capsys.readouterr().out
    assert "[config warning]" in out
    assert "still a placeholder" in out


def test_missing_required_keys_aggregated_into_one_error():
    cfg = _cfg_copy()
    for key in ("SampleListFile", "results_dir", "genome", "trim", "bismark"):
        cfg.pop(key)
    with pytest.raises(WorkflowError) as excinfo:
        validate_config(cfg)
    msg = str(excinfo.value)
    # 5 missing-key issues + the follow-up type issues for trim/bismark/genome
    assert "config validation failed (8 issues)" in msg
    for key in ("SampleListFile", "results_dir", "genome", "trim", "bismark"):
        assert f"missing required config key: {key}" in msg


def test_multiple_error_kinds_aggregated_in_one_report():
    cfg = _cfg_copy()
    cfg.pop("methylation_extractor")
    cfg["threads"] = 0
    cfg["trim"]["adapter"] = "XYZ"
    with pytest.raises(WorkflowError) as excinfo:
        validate_config(cfg)
    msg = str(excinfo.value)
    # missing methylation_extractor + threads + adapter + the follow-up
    # "methylation_extractor must be a mapping" issue
    assert "config validation failed (4 issues)" in msg
    assert "missing required config key: methylation_extractor" in msg
    assert "threads must be an integer >= 1" in msg
    assert "trim.adapter must be a non-empty ACGTN string" in msg


@pytest.mark.parametrize("bad", ["12", True, 0, -2, 2.5])
def test_bad_threads_rejected(bad):
    cfg = _cfg_copy()
    cfg["threads"] = bad
    with pytest.raises(WorkflowError) as excinfo:
        validate_config(cfg)
    assert "threads must be an integer >= 1" in str(excinfo.value)


@pytest.mark.parametrize("bad", ["XYZ", "acgt", "AG AT", "", 123])
def test_bad_trim_adapter_rejected(bad):
    cfg = _cfg_copy()
    cfg["trim"]["adapter"] = bad
    with pytest.raises(WorkflowError) as excinfo:
        validate_config(cfg)
    assert "trim.adapter must be a non-empty ACGTN string" in str(excinfo.value)


@pytest.mark.parametrize(
    "key,bad", [("quality", -1), ("min_len", 0), ("stringency", 0), ("quality", True)]
)
def test_bad_trim_integers_rejected(key, bad):
    cfg = _cfg_copy()
    cfg["trim"][key] = bad
    with pytest.raises(WorkflowError) as excinfo:
        validate_config(cfg)
    assert f"trim.{key} must be an integer" in str(excinfo.value)


@pytest.mark.parametrize("bad", [0, 1.5, "abc", None])
def test_bad_trim_error_rate_rejected(bad):
    cfg = _cfg_copy()
    cfg["trim"]["error_rate"] = bad
    with pytest.raises(WorkflowError) as excinfo:
        validate_config(cfg)
    assert "trim.error_rate must be" in str(excinfo.value)


@pytest.mark.parametrize(
    "key,bad", [("cx_report", "false"), ("merge_cpg", 1), ("cx_report", 0)]
)
def test_bad_methylation_extractor_switches_rejected(key, bad):
    cfg = _cfg_copy()
    cfg["methylation_extractor"][key] = bad
    with pytest.raises(WorkflowError) as excinfo:
        validate_config(cfg)
    assert f"methylation_extractor.{key} must be true/false" in str(excinfo.value)


@pytest.mark.parametrize("bad", [0, True, 1.5, "4"])
def test_bad_methylation_extractor_buffer_frac_rejected(bad):
    cfg = _cfg_copy()
    cfg["methylation_extractor"]["buffer_frac"] = bad
    with pytest.raises(WorkflowError) as excinfo:
        validate_config(cfg)
    assert "methylation_extractor.buffer_frac must be an integer >= 1" in str(
        excinfo.value
    )


def test_bad_bismark_align_extra_rejected():
    cfg = _cfg_copy()
    cfg["bismark"]["align_extra"] = 33
    with pytest.raises(WorkflowError) as excinfo:
        validate_config(cfg)
    assert "bismark must be a mapping and bismark.align_extra a string" in str(
        excinfo.value
    )


def test_resources_unknown_field_rejected():
    cfg = _cfg_copy()
    cfg["resources"] = {"trim": {"cpus": 2}}
    with pytest.raises(WorkflowError) as excinfo:
        validate_config(cfg)
    assert "resources.trim.cpus: unknown field (supported: threads/mem_mb/runtime_min)" in str(
        excinfo.value
    )


def test_resources_non_integer_value_rejected():
    cfg = _cfg_copy()
    cfg["resources"] = {"bismark_align": {"mem_mb": "32000"}}
    with pytest.raises(WorkflowError) as excinfo:
        validate_config(cfg)
    assert "resources.bismark_align.mem_mb must be an integer >= 1" in str(
        excinfo.value
    )


def test_resources_valid_override_accepted(capsys):
    cfg = _cfg_copy()
    cfg["resources"] = {"trim": {"threads": 8, "mem_mb": 8192, "runtime_min": 60}}
    validate_config(cfg)
    assert capsys.readouterr().out == ""


# ---------------------------------------------------------------------
# dmr config section (v0.2, default off)
# ---------------------------------------------------------------------
def test_dmr_defaults_pass_silently(capsys):
    validate_config(_cfg_copy())  # GOOD_CFG ships the full dmr: defaults
    assert capsys.readouterr().out == ""


def test_dmr_absent_section_still_validates(capsys):
    """Pre-v0.2 project configs carry no dmr key: absent stays legal (and
    means disabled)."""
    cfg = _cfg_copy()
    del cfg["dmr"]
    validate_config(cfg)
    assert capsys.readouterr().out == ""


def test_dmr_non_mapping_rejected():
    cfg = _cfg_copy()
    cfg["dmr"] = True
    with pytest.raises(WorkflowError) as excinfo:
        validate_config(cfg)
    assert "dmr must be a mapping with sub-keys" in str(excinfo.value)


@pytest.mark.parametrize(
    "key,bad,needle",
    [
        ("enabled", "yes", "dmr.enabled must be true/false"),
        ("control_group", "", "dmr.control_group must be a non-empty string"),
        ("qvalue", 0, "dmr.qvalue must be"),
        ("qvalue", 1.5, "dmr.qvalue must be"),
        ("min_diff", -1, "dmr.min_diff must be a number in [0, 100]"),
        ("min_diff", 101, "dmr.min_diff must be a number in [0, 100]"),
        ("tile_len", 0, "dmr.tile_len must be an integer >= 1"),
        ("tile_step", True, "dmr.tile_step must be an integer >= 1"),
        ("min_cpg", "3", "dmr.min_cpg must be an integer >= 1"),
        ("min_cov", 0, "dmr.min_cov must be an integer >= 1"),
        ("max_cov", 1.5, "dmr.max_cov must be an integer >= 1"),
        ("batch_correction", "X", "dmr.batch_correction must be 'T' or 'F'"),
    ],
)
def test_dmr_bad_value_rejected(key, bad, needle):
    cfg = _cfg_copy()
    cfg["dmr"][key] = bad
    with pytest.raises(WorkflowError) as excinfo:
        validate_config(cfg)
    assert needle in str(excinfo.value)


def test_dmr_tile_step_above_tile_len_rejected():
    cfg = _cfg_copy()
    cfg["dmr"]["tile_step"] = 2000
    with pytest.raises(WorkflowError) as excinfo:
        validate_config(cfg)
    assert "dmr.tile_step must be <= dmr.tile_len" in str(excinfo.value)


def test_dmr_partial_section_aggregates_all_missing_keys():
    cfg = _cfg_copy()
    cfg["dmr"] = {"enabled": False}
    with pytest.raises(WorkflowError) as excinfo:
        validate_config(cfg)
    msg = str(excinfo.value)
    # every sub-key except enabled is reported in the single aggregated error
    assert "config validation failed (9 issues)" in msg
    for key in ("control_group", "qvalue", "min_diff", "tile_len", "tile_step",
                "min_cpg", "min_cov", "max_cov", "batch_correction"):
        assert f"dmr.{key} " in msg


# ---------------------------------------------------------------------
# dmr design validation (sample-table group/batch vs the dmr section)
# ---------------------------------------------------------------------
GROUPS_2X2 = {"s1": "control", "s2": "control", "t1": "heat", "t2": "heat"}
BATCH_2X2 = {"s1": "b1", "s2": "b2", "t1": "b1", "t2": "b2"}


def _design_cfg(**over):
    cfg = _cfg_copy()
    cfg["dmr"] = {**DMR_OFF, "enabled": True, **over}
    return cfg


def test_dmr_design_valid_2x2_passes():
    validate_dmr_design(_design_cfg(), GROUPS_2X2, BATCH_2X2)  # must NOT raise


def test_dmr_design_requires_group_column():
    single = {sid: None for sid in ("s1", "s2", "t1", "t2")}
    with pytest.raises(WorkflowError) as excinfo:
        validate_dmr_design(_design_cfg(), single, single)
    msg = str(excinfo.value)
    assert "requires a 'group' column" in msg
    assert "dmr design validation failed (1 issues)" in msg


def test_dmr_design_missing_control_group():
    groups = {"s1": "wt", "s2": "wt", "t1": "heat", "t2": "heat"}
    with pytest.raises(WorkflowError) as excinfo:
        validate_dmr_design(_design_cfg(), groups, {k: None for k in groups})
    assert "dmr.control_group='control' is not a group in the sample table" in str(excinfo.value)


def test_dmr_design_replicates_are_a_hard_error():
    groups = {"s1": "control", "t1": "heat", "t2": "heat"}
    with pytest.raises(WorkflowError) as excinfo:
        validate_dmr_design(_design_cfg(), groups, {k: None for k in groups})
    msg = str(excinfo.value)
    assert "every group needs >= 2 replicates for methylKit" in msg
    assert "'control'" in msg  # the lonely group is named


def test_dmr_design_no_non_control_group():
    groups = {"s1": "control", "s2": "control"}
    with pytest.raises(WorkflowError) as excinfo:
        validate_dmr_design(_design_cfg(), groups, {k: None for k in groups})
    assert "no non-control group found" in str(excinfo.value)


def test_dmr_design_batch_correction_requires_batch_column():
    no_batch = {k: None for k in GROUPS_2X2}
    with pytest.raises(WorkflowError) as excinfo:
        validate_dmr_design(_design_cfg(batch_correction="T"), GROUPS_2X2, no_batch)
    assert "dmr.batch_correction='T' requires a 'batch' column" in str(excinfo.value)


def test_dmr_design_batch_correction_with_batch_column_passes():
    validate_dmr_design(
        _design_cfg(batch_correction="T"), GROUPS_2X2, BATCH_2X2
    )  # must NOT raise


def test_dmr_design_requires_merged_cpg_tables():
    cfg = _design_cfg()
    cfg["methylation_extractor"]["merge_cpg"] = False
    with pytest.raises(WorkflowError) as excinfo:
        validate_dmr_design(cfg, GROUPS_2X2, BATCH_2X2)
    assert "set methylation_extractor.merge_cpg: true" in str(excinfo.value)


def test_dmr_design_aggregates_multiple_issues():
    cfg = _design_cfg()
    cfg["methylation_extractor"]["merge_cpg"] = False
    single = {sid: None for sid in ("s1", "s2")}
    with pytest.raises(WorkflowError) as excinfo:
        validate_dmr_design(cfg, single, single)
    msg = str(excinfo.value)
    assert "dmr design validation failed (2 issues)" in msg
    assert "requires a 'group' column" in msg
    assert "set methylation_extractor.merge_cpg: true" in msg


# ---------------------------------------------------------------------
# apply_species_presets
# ---------------------------------------------------------------------
def test_preset_fills_only_unset_keys():
    cfg = {"genome": "/custom/genome.fa", "threads": 8}
    ret = apply_species_presets(
        cfg, {"genome": "/preset/genome.fa", "label": "Rice", "threads": 4}
    )
    assert ret is None  # mutates in place, returns None
    assert cfg == {"genome": "/custom/genome.fa", "threads": 8, "label": "Rice"}


def test_preset_empty_dict_is_noop():
    cfg = {"genome": "/x.fa"}
    assert apply_species_presets(cfg, {}) is None
    assert cfg == {"genome": "/x.fa"}


def test_preset_none_is_noop():
    cfg = {"genome": "/x.fa"}
    assert apply_species_presets(cfg, None) is None
    assert cfg == {"genome": "/x.fa"}


def test_preset_falsy_values_still_fill_unset_keys():
    cfg = {}
    apply_species_presets(cfg, {"genome": "", "flag": 0})
    assert cfg == {"genome": "", "flag": 0}


def test_species_merge_wiring_order_preserved():
    """The module-level merge call must sit at the historical position:
    before the sample table is read and before validate_config (a preset
    genome is what satisfies that validation)."""
    merge = _COMMON_SRC.index("apply_species_presets(config, _SPECIES_PRESET)")
    samples = _COMMON_SRC.index("SAMPLES = load_sample_table(")
    validation = _COMMON_SRC.index("validate_config(config)")
    assert merge < samples < validation
    snakefile = (REPO / "workflow" / "Snakefile").read_text(encoding="utf-8")
    assert '_SPECIES_PRESET = config.get(SPECIES) or {}' in snakefile
    assert '"osa", "hsa", "none"' in snakefile  # allowlist preserved


def test_dmr_wiring_preserved():
    """v0.2 DMR wiring: the RSCRIPT env fallback, the resource defaults
    entry, the TARGETS gating on DMR_ENABLED, and the dmr.smk include placed
    after methylation.smk."""
    assert 'RSCRIPT = os.environ.get("BSSEQ_RSCRIPT", "Rscript")' in _COMMON_SRC
    assert '"dmr": {"threads": 4, "mem_mb": 16000, "runtime_min": 240}' in _COMMON_SRC
    targets = _COMMON_SRC[_COMMON_SRC.index("TARGETS = ["):]
    gate = targets.index("if DMR_ENABLED:")
    assert targets.index('R("6.DMR/flag.log")') > gate  # only requested when enabled
    snakefile = (REPO / "workflow" / "Snakefile").read_text(encoding="utf-8")
    methylation_include = snakefile.index('"rules", "methylation.smk"')
    dmr_include = snakefile.index('"rules", "dmr.smk"')
    assert methylation_include < dmr_include


def test_dmr_rule_gated_on_dmr_enabled():
    """The dmr module must define its rule inside a parse-time guard so the
    default (dmr disabled) DAG is unchanged."""
    dmr_smk = (REPO / "workflow" / "rules" / "dmr.smk").read_text(encoding="utf-8")
    guard = dmr_smk.index("if DMR_ENABLED:")
    rule = dmr_smk.index("rule dmr_methylkit:")
    assert 0 < guard < rule
    assert dmr_smk.count("rule dmr_methylkit:") == 1
    assert "configfile" not in dmr_smk  # no config side effects in the module


# ---------------------------------------------------------------------
# Shared-layer wrapper sanity (imported through their own __file__-relative
# shared-path resolution; no snakemake involvement)
# ---------------------------------------------------------------------
def _import_script(script_name):
    path = REPO / "workflow" / "scripts" / script_name
    spec = importlib.util.spec_from_file_location("bsseq_tests_" + path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_runtime_config_spec_invariants():
    rc = _import_script("runtime_config.py")
    assert rc.SPEC.env_prefix == "BSSEQ"
    assert "rscript" in rc.DEFAULT_TOOLS  # required by the shared framework
    assert "rscript" not in rc.SPEC.pipeline_tools["default"]  # no R preflight
    assert rc.SPEC.r_packages == {}
    assert rc.SPEC.tool_defaults == {"python": "python3"}
    # The shared framework must resolve through the wrapper's own
    # __file__-relative path into THIS repository's shared/ layer
    # (shared/ sits at the repository root, next to the project dirs).
    shared = sys.modules.get("bioworkflows_runtime")
    assert shared is not None
    assert Path(shared.__file__).resolve().parents[1] == REPO.parent / "shared"


def test_collect_versions_invariants():
    cv = _import_script("collect_versions.py")
    assert cv.TOOLS and all(isinstance(tool, str) for tool in cv.TOOLS)
    assert "bismark" in cv.TOOLS and "python" in cv.TOOLS
    assert cv.TOOL_DEFAULTS == {"python": "python3"}
    assert cv.DATABASES == {}
