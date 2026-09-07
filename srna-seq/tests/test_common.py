"""Unit tests for the pure-Python logic of srna-seq (docs/TODO.md item 3).

The functions under test are extracted from the real workflow source and
executed -- never copied -- so the tests stay in sync with the
implementation (same technique as chip_cuttag_atac_faire/tests/run_tests.py):

- workflow/rules/common.smk imports snakemake.exceptions.WorkflowError at
  module level, which plain python cannot import, so only the pure def
  blocks are extracted and exec'd against a WorkflowError stand-in. The
  module-level statements (SAMPLES = ..., validate_config(config)) are
  never executed here.
- The species preset merge exists twice by necessity:
  apply_species_presets() in workflow/rules/common.smk (the unit-tested
  def) and an inline block in workflow/Snakefile -- which must run before
  common.smk is included, because common.smk's parse-time validation
  consumes the merged references (see
  test_validate_config_fails_without_species_merge), and the main
  Snakefile cannot carry top-level function definitions (snakemake 7 lint
  rejects mixed rules and functions in the same snakefile). The parity
  test at the bottom pins the inline block's behavior to the def's.

Run: python -m pytest tests -q  (paths derive from __file__, so the
working directory does not matter).
"""
import copy
import csv
import os
import re

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class WorkflowError(Exception):
    """Stand-in for snakemake.exceptions.WorkflowError (single-str arg)."""


def _segment(src, start_marker, end_markers):
    """Cut one contiguous source block out of a workflow file."""
    start = src.index(start_marker)
    end = len(src)
    for marker in end_markers:
        pos = src.find(marker, start + len(start_marker))
        if pos != -1:
            end = min(end, pos)
    return src[start:end]


def _read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def _extract_workflow_source():
    """Exec the pure def blocks of common.smk into one namespace with a
    WorkflowError stand-in and the globals they need."""
    common_src = _read(os.path.join(REPO, "workflow", "rules", "common.smk"))
    blocks = [
        # Module-level regex compile; referenced by the extracted defs.
        _segment(common_src, "_NAME_RE = ", ["\nSCRIPTS = "]),
        _segment(common_src, "def _resolve_sample_table(", ["\n_SAMPLE_TABLE = "]),
        _segment(common_src, "def _validate_deg_design(", ["\ndef validate_config("]),
        _segment(common_src, "def validate_config(",
                 ["\nvalidate_config(config)", "\ndef "]),
        _segment(common_src, "def apply_species_presets(",
                 ["\ndef cascade_classes("]),
    ]
    ns = {"csv": csv, "os": os, "re": re, "WorkflowError": WorkflowError,
          "BASE_DIR": REPO,
          # Module-level sample annotations; load_sample_table populates them
          # in the same single parse that returns the sample-id list.
          "SAMPLE_GROUPS": {}, "SAMPLE_BATCH": {}}
    for block in blocks:
        exec(compile(block, "common.smk(extracted)", "exec"), ns)
    return ns


WF = _extract_workflow_source()
_NAME_RE = WF["_NAME_RE"]
_resolve_sample_table = WF["_resolve_sample_table"]
load_sample_table = WF["load_sample_table"]
validate_config = WF["validate_config"]
_validate_deg_design = WF["_validate_deg_design"]
apply_species_presets = WF["apply_species_presets"]
SAMPLE_GROUPS = WF["SAMPLE_GROUPS"]
SAMPLE_BATCH = WF["SAMPLE_BATCH"]


# ---------------------------------------------------------------------
# load_sample_table
# ---------------------------------------------------------------------
def _write_table(tmp_path, rows):
    """rows[0] is the header; returns the file path as str."""
    path = tmp_path / "samples.csv"
    with open(path, "w", newline="", encoding="utf-8") as fh:
        csv.writer(fh).writerows(rows)
    return str(path)


def _expect_table_error(tmp_path, rows, needle):
    path = _write_table(tmp_path, rows)
    with pytest.raises(WorkflowError) as excinfo:
        load_sample_table(path)
    assert needle in str(excinfo.value)


def test_load_sample_table_real_repository_file():
    assert load_sample_table(os.path.join(REPO, "config", "samples.csv")) == ["s1", "s2"]


def test_load_sample_table_keeps_order_and_accepts_legal_names(tmp_path):
    path = _write_table(tmp_path, [["sample_id"], ["s1"], ["A.b-2"], ["x_3"]])
    assert load_sample_table(path) == ["s1", "A.b-2", "x_3"]


def test_load_sample_table_rejects_bad_header(tmp_path):
    _expect_table_error(tmp_path, [["sample"], ["s1"]],
                        "must have exactly one of these headers")


def test_load_sample_table_single_column_builds_none_annotations(tmp_path):
    """The v0.1 single-column contract is unchanged; the annotation dicts are
    populated with None values in the same single parse."""
    path = _write_table(tmp_path, [["sample_id"], ["s1"], ["s2"]])
    assert load_sample_table(path) == ["s1", "s2"]
    assert SAMPLE_GROUPS == {"s1": None, "s2": None}
    assert SAMPLE_BATCH == {"s1": None, "s2": None}


def test_load_sample_table_group_column_parses(tmp_path):
    path = _write_table(tmp_path, [["sample_id", "group"],
                                   ["s1", "control"], ["s2", "treat"]])
    assert load_sample_table(path) == ["s1", "s2"]
    assert SAMPLE_GROUPS == {"s1": "control", "s2": "treat"}
    assert SAMPLE_BATCH == {"s1": None, "s2": None}


def test_load_sample_table_group_and_batch_columns_parse(tmp_path):
    path = _write_table(tmp_path, [["sample_id", "group", "batch"],
                                   ["s1", "control", "b1"],
                                   ["s2", "treat", "b2"]])
    assert load_sample_table(path) == ["s1", "s2"]
    assert SAMPLE_GROUPS == {"s1": "control", "s2": "treat"}
    assert SAMPLE_BATCH == {"s1": "b1", "s2": "b2"}


def test_load_sample_table_duplicate_group_values_are_legal(tmp_path):
    """Duplicates are only rejected on sample_id, never on group/batch."""
    path = _write_table(tmp_path, [["sample_id", "group"],
                                   ["s1", "control"], ["s2", "control"]])
    assert load_sample_table(path) == ["s1", "s2"]
    assert SAMPLE_GROUPS == {"s1": "control", "s2": "control"}


def test_load_sample_table_rejects_batch_without_group(tmp_path):
    _expect_table_error(tmp_path, [["sample_id", "batch"], ["s1", "b1"]],
                        "must have exactly one of these headers")


def test_load_sample_table_rejects_wrong_column_order(tmp_path):
    _expect_table_error(tmp_path, [["group", "sample_id"], ["control", "s1"]],
                        "must have exactly one of these headers")


def test_load_sample_table_rejects_unknown_extra_column(tmp_path):
    _expect_table_error(tmp_path, [["sample_id", "group", "strain"],
                                   ["s1", "control", "x"]],
                        "must have exactly one of these headers")


def test_load_sample_table_rejects_empty_group(tmp_path):
    _expect_table_error(tmp_path, [["sample_id", "group"], ["s1", ""]],
                        "group must not be empty")


def test_load_sample_table_rejects_illegal_group(tmp_path):
    _expect_table_error(tmp_path, [["sample_id", "group"], ["s1", "a b"]],
                        "group='a b' contains illegal characters")


def test_load_sample_table_rejects_double_underscore_group(tmp_path):
    _expect_table_error(tmp_path, [["sample_id", "group"], ["s1", "a__b"]],
                        "group='a__b' contains illegal characters")


def test_load_sample_table_rejects_illegal_batch(tmp_path):
    _expect_table_error(tmp_path,
                        [["sample_id", "group", "batch"], ["s1", "control", "-b1"]],
                        "batch='-b1' contains illegal characters")


def test_load_sample_table_duplicates_still_checked_on_sample_id_only(tmp_path):
    _expect_table_error(tmp_path,
                        [["sample_id", "group"], ["s1", "control"], ["s1", "treat"]],
                        "duplicate sample_id")


def test_load_sample_table_rejects_empty_id(tmp_path):
    _expect_table_error(tmp_path, [["sample_id"], ["   "]], "must not be empty")


def test_load_sample_table_rejects_illegal_characters(tmp_path):
    _expect_table_error(tmp_path, [["sample_id"], ["s 1"]], "illegal characters")


def test_load_sample_table_rejects_leading_dash(tmp_path):
    _expect_table_error(tmp_path, [["sample_id"], ["-s1"]], "illegal characters")


def test_load_sample_table_rejects_double_underscore(tmp_path):
    _expect_table_error(tmp_path, [["sample_id"], ["s__1"]], "illegal characters")


def test_load_sample_table_rejects_duplicates(tmp_path):
    _expect_table_error(tmp_path, [["sample_id"], ["s1"], ["s1"]],
                        "duplicate sample_id")


def test_load_sample_table_rejects_table_without_data_rows(tmp_path):
    _expect_table_error(tmp_path, [["sample_id"]], "no data rows")


# ---------------------------------------------------------------------
# _resolve_sample_table
# ---------------------------------------------------------------------
def test_resolve_sample_table_returns_absolute_path_unchanged(tmp_path):
    target = str(tmp_path / "samples.csv")
    assert _resolve_sample_table(target) == target


def test_resolve_sample_table_prefers_existing_cwd_relative_file(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "samples.csv").write_text("sample_id\ns1\n", encoding="utf-8")
    assert _resolve_sample_table("samples.csv") == str(tmp_path / "samples.csv")


def test_resolve_sample_table_falls_back_to_base_dir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)  # empty directory: the name exists nowhere
    missing = "no_such_samples_table_9f3a.csv"
    assert _resolve_sample_table(missing) == os.path.join(REPO, missing)


# ---------------------------------------------------------------------
# validate_config
# ---------------------------------------------------------------------
GOOD_CFG = {
    "SampleListFile": "config/samples.csv",
    "results_dir": "results",
    "threads": 12,
    "trim": {"quality": 25, "min_len": 15, "stringency": 3, "error_rate": 0.1,
             "adapter": "AGATCGGAAGAGC", "extra": ""},
    "cascade": [
        {"name": "rRNA", "fasta": "/ref/rRNA.fa"},
        {"name": "miRNA", "fasta": "/ref/miRNA.fa"},
    ],
    "genome": {"fasta": "/ref/genome.fa"},
    "bowtie": {"extra": ""},
}


def _validated_errors(cfg):
    """Run the extracted validate_config; return the aggregated message
    or None when the config is accepted."""
    try:
        validate_config(cfg)
    except WorkflowError as exc:
        return str(exc)
    return None


def test_validate_config_accepts_a_valid_config():
    assert _validated_errors(copy.deepcopy(GOOD_CFG)) is None


def test_validate_config_aggregates_all_errors_into_one_report():
    cfg = copy.deepcopy(GOOD_CFG)
    del cfg["SampleListFile"]
    del cfg["results_dir"]
    cfg["threads"] = 0
    message = _validated_errors(cfg)
    assert message is not None
    assert message.startswith("config validation failed (3 issues):")
    for needle in ("missing required config key: SampleListFile",
                   "missing required config key: results_dir",
                   "threads must be an integer >= 1"):
        assert needle in message


@pytest.mark.parametrize("bad_threads", [0, -2, "12", True, 1.5])
def test_validate_config_rejects_bad_threads(bad_threads):
    cfg = copy.deepcopy(GOOD_CFG)
    cfg["threads"] = bad_threads
    assert "threads must be an integer >= 1" in _validated_errors(cfg)


def test_validate_config_rejects_non_string_results_dir():
    cfg = copy.deepcopy(GOOD_CFG)
    cfg["results_dir"] = 5
    assert "results_dir must be a string" in _validated_errors(cfg)


def test_validate_config_missing_cascade_key_is_reported():
    cfg = copy.deepcopy(GOOD_CFG)
    del cfg["cascade"]
    message = _validated_errors(cfg)
    # Missing key AND shape are both reported (cfg.get("cascade") is None,
    # so the non-empty-list check fires as well) -- pin this behavior.
    assert "missing required config key: cascade" in message
    assert "cascade must be a non-empty ordered list" in message


def test_validate_config_rejects_empty_cascade():
    cfg = copy.deepcopy(GOOD_CFG)
    cfg["cascade"] = []
    assert "cascade must be a non-empty ordered list" in _validated_errors(cfg)


def test_validate_config_rejects_non_list_cascade():
    cfg = copy.deepcopy(GOOD_CFG)
    cfg["cascade"] = "rRNA"
    assert "cascade must be a non-empty ordered list" in _validated_errors(cfg)


def test_validate_config_rejects_non_mapping_cascade_entry():
    cfg = copy.deepcopy(GOOD_CFG)
    cfg["cascade"] = ["rRNA"]
    assert "cascade[0] must be a mapping" in _validated_errors(cfg)


def test_validate_config_rejects_illegal_cascade_name():
    cfg = copy.deepcopy(GOOD_CFG)
    cfg["cascade"] = [{"name": "a__b", "fasta": "/ref/a.fa"}]
    assert "cascade[0].name='a__b' is not a legal class name" in _validated_errors(cfg)


def test_validate_config_rejects_duplicate_cascade_names():
    cfg = copy.deepcopy(GOOD_CFG)
    cfg["cascade"] = [{"name": "rRNA", "fasta": "/ref/1.fa"},
                      {"name": "rRNA", "fasta": "/ref/2.fa"}]
    assert "duplicates an earlier entry" in _validated_errors(cfg)


def test_validate_config_rejects_non_string_cascade_fasta():
    cfg = copy.deepcopy(GOOD_CFG)
    cfg["cascade"] = [{"name": "rRNA", "fasta": 5}]
    assert "cascade[0].fasta must be a string" in _validated_errors(cfg)


def test_validate_config_requires_at_least_one_configured_fasta():
    cfg = copy.deepcopy(GOOD_CFG)
    cfg["cascade"] = [{"name": "rRNA", "fasta": ""}]
    cfg["genome"] = {"fasta": ""}
    assert "nothing to align" in _validated_errors(cfg)


@pytest.mark.parametrize("bad_adapter", ["", "AGATCGX", "agatcg", 123, None])
def test_validate_config_rejects_bad_adapter(bad_adapter):
    cfg = copy.deepcopy(GOOD_CFG)
    cfg["trim"]["adapter"] = bad_adapter
    assert "trim.adapter must be a non-empty ACGTN string" in _validated_errors(cfg)


@pytest.mark.parametrize("field,lo", [("quality", 0), ("min_len", 1), ("stringency", 1)])
def test_validate_config_rejects_bad_trim_integers(field, lo):
    cfg = copy.deepcopy(GOOD_CFG)
    cfg["trim"][field] = lo - 1
    assert f"trim.{field} must be an integer >= {lo}" in _validated_errors(cfg)


@pytest.mark.parametrize("bad", [0, 1.5, -0.1, "abc", None])
def test_validate_config_rejects_bad_error_rate(bad):
    cfg = copy.deepcopy(GOOD_CFG)
    cfg["trim"]["error_rate"] = bad
    assert "trim.error_rate" in _validated_errors(cfg)


def test_validate_config_rejects_non_string_trim_extra():
    cfg = copy.deepcopy(GOOD_CFG)
    cfg["trim"]["extra"] = ["-q 5"]
    assert "trim.extra must be a string" in _validated_errors(cfg)


def test_validate_config_rejects_non_mapping_trim():
    cfg = copy.deepcopy(GOOD_CFG)
    cfg["trim"] = "fast"
    assert "trim must be a mapping with sub-keys" in _validated_errors(cfg)


def test_validate_config_rejects_non_string_bowtie_extra():
    cfg = copy.deepcopy(GOOD_CFG)
    cfg["bowtie"] = {"extra": 5}
    assert "bowtie must be a mapping and bowtie.extra a string" in _validated_errors(cfg)


def test_validate_config_rejects_non_mapping_resources():
    cfg = copy.deepcopy(GOOD_CFG)
    cfg["resources"] = ["trim"]
    assert "resources must be a mapping" in _validated_errors(cfg)


def test_validate_config_rejects_non_mapping_resource_entry():
    cfg = copy.deepcopy(GOOD_CFG)
    cfg["resources"] = {"trim": 4}
    assert "resources.trim must be a mapping" in _validated_errors(cfg)


def test_validate_config_rejects_unknown_resource_field():
    cfg = copy.deepcopy(GOOD_CFG)
    cfg["resources"] = {"trim": {"gpu": 1}}
    assert ("resources.trim.gpu: unknown field "
            "(supported: threads/mem_mb/runtime_min)") in _validated_errors(cfg)


def test_validate_config_accepts_known_resource_fields():
    cfg = copy.deepcopy(GOOD_CFG)
    cfg["resources"] = {"trim": {"threads": 4, "mem_mb": 4096, "runtime_min": 60}}
    assert _validated_errors(cfg) is None


def test_validate_config_warns_for_placeholders_and_skipped_classes(capsys):
    cfg = copy.deepcopy(GOOD_CFG)
    cfg["cascade"] = [
        {"name": "rRNA", "fasta": "/path/to/reference/rRNA.fa"},
        {"name": "miRNA", "fasta": ""},
    ]
    cfg["genome"] = {"fasta": "/path/to/reference/genome.fa"}
    assert _validated_errors(cfg) is None
    out = capsys.readouterr().out
    assert ("[config warning] cascade.rRNA.fasta is still a placeholder: "
            "/path/to/reference/rRNA.fa") in out
    assert "[config warning] cascade class 'miRNA' has no fasta and is skipped" in out
    assert ("[config warning] genome.fasta is still a placeholder: "
            "/path/to/reference/genome.fa") in out


def test_validate_config_fails_without_species_merge():
    """Documents why apply_species_presets must run before common.smk is
    included: the repository default config leaves every cascade and
    genome fasta empty, so parse-time validation of the un-merged config
    raises 'nothing to align'."""
    cfg = {
        "SampleListFile": "config/samples.csv",
        "results_dir": "results",
        "threads": 12,
        "trim": dict(GOOD_CFG["trim"]),
        "cascade": [{"name": name, "fasta": ""} for name in
                    ("rRNA", "snoRNA", "snRNA", "tRNA", "miRNA", "mRNA", "rhizo")],
        "genome": {"fasta": ""},
        "bowtie": {"extra": ""},
    }
    assert "nothing to align" in _validated_errors(cfg)


# ---------------------------------------------------------------------
# deg section: optional differential-expression stage (default off)
# ---------------------------------------------------------------------
DEG_ENABLED = {"enabled": True, "classes": ["miRNA"], "control_group": "control",
               "foldchange": 2, "padj": 0.05, "batch_correction": "F",
               "pca_ntop": 2000}


def _set_annotations(groups, batches=None):
    """Seed the module-level sample-annotation dicts as load_sample_table
    would (used to test the design checks without writing table files)."""
    SAMPLE_GROUPS.clear()
    SAMPLE_GROUPS.update(groups)
    SAMPLE_BATCH.clear()
    SAMPLE_BATCH.update(batches if batches is not None else {k: None for k in groups})


def _deg_cfg(**overrides):
    """GOOD_CFG with the full deg contract section, deg fields overridable."""
    cfg = copy.deepcopy(GOOD_CFG)
    cfg["deg"] = {**DEG_ENABLED, **overrides}
    return cfg


def test_validate_config_accepts_full_deg_defaults_when_disabled():
    _set_annotations({"s1": None, "s2": None})  # no group column needed
    assert _validated_errors(_deg_cfg(enabled=False)) is None


def test_validate_config_rejects_non_mapping_deg():
    cfg = copy.deepcopy(GOOD_CFG)
    cfg["deg"] = ["enabled"]
    assert "deg must be a mapping" in _validated_errors(cfg)


@pytest.mark.parametrize("bad", ["yes", 1, None])
def test_validate_config_rejects_bad_deg_enabled(bad):
    assert "deg.enabled must be a boolean" in _validated_errors(_deg_cfg(enabled=bad))


@pytest.mark.parametrize("bad", [[], "miRNA", ["miRNA", "a__b"], ["x y"], [5]])
def test_validate_config_rejects_bad_deg_classes(bad):
    assert "deg.classes must be a non-empty list of legal class names" in _validated_errors(_deg_cfg(classes=bad))


@pytest.mark.parametrize("bad", [1, 0.5, "2", True, None])
def test_validate_config_rejects_bad_deg_foldchange(bad):
    assert "deg.foldchange must be a number > 1" in _validated_errors(_deg_cfg(foldchange=bad))


@pytest.mark.parametrize("bad", [0, 1, -0.1, "0.05", True])
def test_validate_config_rejects_bad_deg_padj(bad):
    assert "deg.padj must be in (0, 1)" in _validated_errors(_deg_cfg(padj=bad))


@pytest.mark.parametrize("bad", ["true", False, 1, None])
def test_validate_config_rejects_bad_deg_batch_correction(bad):
    assert "deg.batch_correction must be 'T' or 'F'" in _validated_errors(_deg_cfg(batch_correction=bad))


@pytest.mark.parametrize("bad", [0, 1.5, "2000", True, None])
def test_validate_config_rejects_bad_deg_pca_ntop(bad):
    assert "deg.pca_ntop must be an integer >= 1" in _validated_errors(_deg_cfg(pca_ntop=bad))


def test_validate_config_warns_deg_class_not_in_cascade(capsys):
    _set_annotations({"s1": "control", "s2": "treat"})
    assert _validated_errors(_deg_cfg(classes=["piRNA"])) is None
    out = capsys.readouterr().out
    assert ("[config warning] deg class 'piRNA' is not a configured cascade "
            "class and will be skipped") in out


def test_validate_config_deg_design_requires_group_column():
    _set_annotations({"s1": None, "s2": None})
    message = _validated_errors(_deg_cfg())
    assert "deg.enabled is true but the sample table has no 'group' column" in message


def test_validate_config_deg_design_requires_control_among_groups():
    _set_annotations({"s1": "treat", "s2": "treat2"})
    message = _validated_errors(_deg_cfg())
    assert ("deg.control_group='control' is not among the sample-table "
            "group values ['treat', 'treat2']") in message


def test_validate_config_deg_design_requires_two_distinct_groups():
    _set_annotations({"s1": "control", "s2": "control"})
    message = _validated_errors(_deg_cfg())
    assert "deg design needs at least 2 distinct groups" in message


def test_validate_config_deg_design_warns_on_single_replicate_group(capsys):
    _set_annotations({"s1": "control", "s2": "treat"})
    assert _validated_errors(_deg_cfg()) is None
    out = capsys.readouterr().out
    assert "[config warning] deg design: group 'control' has 1 replicate(s)" in out
    assert "[config warning] deg design: group 'treat' has 1 replicate(s)" in out


def test_validate_config_deg_design_accepts_replicated_groups():
    _set_annotations({"s1": "control", "s2": "control", "s3": "treat", "s4": "treat"})
    assert _validated_errors(_deg_cfg()) is None


def test_validate_config_deg_batch_correction_requires_batch_column():
    _set_annotations({"s1": "control", "s2": "treat"})
    message = _validated_errors(_deg_cfg(batch_correction="T"))
    assert "deg.batch_correction is 'T' but the sample table has no 'batch' column" in message


def test_validate_config_deg_batch_correction_accepts_batch_column():
    _set_annotations({"s1": "control", "s2": "treat"}, {"s1": "b1", "s2": "b2"})
    assert _validated_errors(_deg_cfg(batch_correction="T")) is None


def test_validate_config_deg_design_aggregates_with_other_errors():
    """Design errors land in the same aggregated report as shape errors."""
    _set_annotations({"s1": None, "s2": None})
    cfg = _deg_cfg(foldchange=0)
    message = _validated_errors(cfg)
    assert message.startswith("config validation failed (2 issues):")
    assert "deg.foldchange must be a number > 1" in message
    assert "deg.enabled is true but the sample table has no 'group' column" in message


# ---------------------------------------------------------------------
# apply_species_presets (real source: workflow/rules/common.smk)
# ---------------------------------------------------------------------
def test_apply_species_presets_fills_empty_cascade_fasta_only():
    cfg = {"cascade": [{"name": "rRNA", "fasta": ""},
                       {"name": "miRNA", "fasta": "/own/miRNA.fa"},
                       {"name": "rhizo", "fasta": ""}],
           "genome": {"fasta": "/own/genome.fa"}}
    preset = {"cascade": {"rRNA": "/osa/rRNA.fa", "snoRNA": "/osa/snoRNA.fa"}}
    assert apply_species_presets(cfg, preset) is None
    assert cfg["cascade"][0]["fasta"] == "/osa/rRNA.fa"   # empty -> preset path
    assert cfg["cascade"][1]["fasta"] == "/own/miRNA.fa"  # explicit value wins
    assert cfg["cascade"][2]["fasta"] == ""               # class absent from preset
    assert cfg["genome"]["fasta"] == "/own/genome.fa"     # non-empty genome untouched


def test_apply_species_presets_fills_empty_genome_fasta():
    cfg = {"cascade": [{"name": "rRNA", "fasta": "/ref/rRNA.fa"}],
           "genome": {"fasta": ""}}
    apply_species_presets(cfg, {"genome": "/osa/genome.fa"})
    assert cfg["genome"]["fasta"] == "/osa/genome.fa"


def test_apply_species_presets_creates_genome_key_when_missing():
    cfg = {"cascade": [{"name": "rRNA", "fasta": ""}]}
    apply_species_presets(cfg, {"cascade": {"rRNA": "/osa/rRNA.fa"},
                                "genome": "/osa/genome.fa"})
    assert cfg["genome"] == {"fasta": "/osa/genome.fa"}


def test_apply_species_presets_empty_preset_is_a_no_op():
    cfg = {"cascade": [{"name": "rRNA", "fasta": ""}], "genome": {"fasta": ""}}
    before = copy.deepcopy(cfg)
    apply_species_presets(cfg, {})
    assert cfg == before


def test_apply_species_presets_mirrors_real_species_yaml_shape():
    """Pins the exact merge semantics with the real species.yaml preset
    shape: preset['genome'] is a {'fasta': path} mapping and is stored
    verbatim (byte-identical to the historical inline merge)."""
    cfg = {"cascade": [{"name": "rRNA", "fasta": ""}], "genome": {"fasta": ""}}
    preset = {"cascade": {"rRNA": "/path/to/reference/osa/rRNA.fa"},
              "genome": {"fasta": "/path/to/reference/osa/genome.fa"}}
    apply_species_presets(cfg, preset)
    assert cfg["cascade"][0]["fasta"] == "/path/to/reference/osa/rRNA.fa"
    assert cfg["genome"]["fasta"] == {"fasta": "/path/to/reference/osa/genome.fa"}


def test_species_merge_over_real_repository_configs():
    """Integration pin: the real config/species.yaml preset fills every
    empty cascade fasta of the real config/config.yaml, and the merged
    config passes parse-time validation (placeholder warnings only)."""
    yaml = pytest.importorskip("yaml")
    with open(os.path.join(REPO, "config", "species.yaml"), encoding="utf-8") as fh:
        presets = yaml.safe_load(fh) or {}
    with open(os.path.join(REPO, "config", "config.yaml"), encoding="utf-8") as fh:
        project = yaml.safe_load(fh) or {}
    merged = {**presets, **project}
    species = merged.get("species", "osa")
    assert species in ("osa", "none")
    preset = (merged.get(species) or {}) if species != "none" else {}
    apply_species_presets(merged, preset)
    fastas = {entry["name"]: str(entry.get("fasta") or "").strip()
              for entry in merged["cascade"]}
    assert all(fastas.values()), fastas
    try:
        validate_config(merged)
    except WorkflowError as exc:
        pytest.fail(f"merged repository config must validate: {exc}")


# ---------------------------------------------------------------------
# Parity: workflow/Snakefile's inline species merge == apply_species_presets
# ---------------------------------------------------------------------
_OSA_PRESET = {
    "label": "Rice Oryza sativa",
    "cascade": {"rRNA": "/osa/rRNA.fa", "miRNA": "/osa/miRNA.fa"},
    "genome": {"fasta": "/osa/genome.fa"},
}


def _inline_merge_fixtures():
    """Fixture configs for the inline merge in workflow/Snakefile."""
    empty_filled = {
        "species": "osa", "osa": _OSA_PRESET,
        "cascade": [{"name": "rRNA", "fasta": ""},
                    {"name": "miRNA", "fasta": "/own/miRNA.fa"},
                    {"name": "rhizo", "fasta": ""}],
        "genome": {"fasta": ""},
    }
    none_species = copy.deepcopy(empty_filled)
    none_species["species"] = "none"          # "none" must skip the merge
    empty_preset = copy.deepcopy(empty_filled)
    empty_preset["osa"] = {}                  # preset entry present but empty
    return [empty_filled, none_species, empty_preset]


def _snakefile_merge_source():
    src = _read(os.path.join(REPO, "workflow", "Snakefile"))
    return _segment(src, 'SPECIES = config.get("species"',
                    ["\n# Shared definitions", "\nwildcard_constraints:"])


def _run_snakefile_inline_merge(cfg):
    ns = {"config": copy.deepcopy(cfg)}
    exec(compile(_snakefile_merge_source(), "Snakefile(species merge)", "exec"), ns)
    return ns["config"]


@pytest.mark.parametrize("cfg", _inline_merge_fixtures(),
                         ids=["osa-fills-empties", "none-skips-merge",
                              "empty-preset-no-op"])
def test_snakefile_inline_merge_matches_apply_species_presets(cfg):
    expected = copy.deepcopy(cfg)
    species = expected.get("species", "osa")
    preset = (expected.get(species) or {}) if species != "none" else {}
    apply_species_presets(expected, preset)
    assert _run_snakefile_inline_merge(cfg) == expected


def test_snakefile_inline_merge_rejects_unknown_species():
    with pytest.raises(ValueError, match="Choose osa or none"):
        _run_snakefile_inline_merge({"species": "mouse", "cascade": [],
                                     "genome": {"fasta": ""}})


def test_snakefile_inline_merge_fills_real_repository_config():
    """The real config/config.yaml (all fastas empty) merged through the
    inline Snakefile block must match the def-merged result."""
    yaml = pytest.importorskip("yaml")
    with open(os.path.join(REPO, "config", "species.yaml"), encoding="utf-8") as fh:
        presets = yaml.safe_load(fh) or {}
    with open(os.path.join(REPO, "config", "config.yaml"), encoding="utf-8") as fh:
        project = yaml.safe_load(fh) or {}
    merged = {**presets, **project}
    species = merged.get("species", "osa")
    assert species in ("osa", "none")
    expected = copy.deepcopy(merged)
    preset = (expected.get(species) or {}) if species != "none" else {}
    apply_species_presets(expected, preset)
    assert _run_snakefile_inline_merge(merged) == expected
