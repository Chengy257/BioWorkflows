#!/usr/bin/env python3
"""Dependency-free tests: sample-table parsing, routing helpers, config integrity.

Run: python tests/run_tests.py (from the repository root)
- Extracts the [real source] from workflow/rules/common.smk and executes it
  (not a copy), keeping tests in sync with the implementation
- pyyaml is optional: YAML checks are skipped when it is missing (installed in CI)
"""
import csv
import os
import re
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAILED = []


def check(name, cond, detail=""):
    if cond:
        print(f"  PASS  {name}")
    else:
        print(f"  FAIL  {name}  {detail}")
        FAILED.append(name)


# ---------------------------------------------------------------------
# Extract the functions under test from workflow/rules/common.smk (real source)
# ---------------------------------------------------------------------
class WorkflowError(Exception):
    """Stand-in for snakemake.exceptions.WorkflowError (signature-compatible: single str arg)"""


def extract_workflow_functions():
    with open(os.path.join(REPO, "workflow", "rules", "common.smk"), encoding="utf-8") as fh:
        src = fh.read()

    def segment(start_marker, end_markers):
        start = src.index(start_marker)
        end = len(src)
        for m in end_markers:
            pos = src.find(m, start + len(start_marker))
            if pos != -1:
                end = min(end, pos)
        return src[start:end]

    block_a = segment("ASSAYS = (", ["SAMPLES, GROUPS, SEQTYPE_OF"])
    block_b = segment("def _group_regex(", ["def _group_regex_dummy_never",
                                            "\n# -----",
                                            "\nBAM_TARGETS"])
    # Cut at the end of the _group_regex body (next top-level def or comment divider)
    m = re.search(r"\n(?=def |# -|BAM_TARGETS)", block_b[len("def _group_regex("):])
    if m:
        block_b = block_b[: len("def _group_regex(") + m.start() + 1]

    ns = {"csv": csv, "os": os, "re": re, "WorkflowError": WorkflowError,
          "BASE_DIR": REPO}
    exec(compile(block_a + "\n\n" + block_b, "common.smk(extracted)", "exec"), ns)
    return ns


WF = extract_workflow_functions()
load_sample_table = WF["load_sample_table"]
_resolve_sample_table = WF["_resolve_sample_table"]
_group_regex = WF["_group_regex"]


def write_csv(rows):
    """rows[0] is the header; returns the temp file path"""
    fd, path = tempfile.mkstemp(suffix=".csv")
    with os.fdopen(fd, "w", newline="", encoding="utf-8") as fh:
        csv.writer(fh).writerows(rows)
    return path


def expect_error(name, rows, needle):
    path = write_csv(rows)
    try:
        load_sample_table(path)
        check(name, False, "WorkflowError was not raised")
    except WorkflowError as e:
        check(name, needle in str(e), f"error message does not contain {needle!r}: {e}")
    finally:
        os.unlink(path)


HEADER = ["sample_id", "role", "group", "seqtype", "layout", "peak_type"]

# ---------------------------------------------------------------------
print("== 1. Sample-table parsing (real example file) ==")
samples, groups, seqtype_of = load_sample_table(
    os.path.join(REPO, "config", "samples.csv"))
check("example: 8 samples, deduplicated, order preserved",
      samples == ["myc", "IgG", "H3K27me3_rep1", "IgG_cuta",
                  "atac_leaf_1", "atac_leaf_2", "faire_root", "Input_faire"],
      str(samples))
check("example: myc_vs_IgG structure",
      groups["myc_vs_IgG"]["treat"] == ["myc"]
      and groups["myc_vs_IgG"]["control"] == ["IgG"]
      and groups["myc_vs_IgG"]["seqtype"] == "chip"
      and groups["myc_vs_IgG"]["peak_type"] == "narrow")
check("example: atac_leaf without control is valid",
      groups["atac_leaf"]["control"] == [] and groups["atac_leaf"]["seqtype"] == "atac")
check("example: faire group seqtype mapping", seqtype_of["faire_root"] == "faire")

s2, g2, st2 = load_sample_table(os.path.join(REPO, "example", "samples.csv"))
check("example/samples.csv: real 2-sample table parses",
      s2 == ["myc", "IgG"] and g2["myc_vs_IgG"]["peak_type"] == "narrow")

print("== 2. Sample-table validation error paths ==")
expect_error("missing columns reported", HEADER[:4] + [["a", "treat", "g", "chip"]], "missing columns")
expect_error("invalid role", [HEADER, ["a", "TREATMENT", "g", "chip", "PE", "narrow"]], "role")
expect_error("invalid seqtype", [HEADER, ["a", "treat", "g", "chipseq", "PE", "narrow"]], "seqtype")
expect_error("SE layout rejected", [HEADER, ["a", "treat", "g", "chip", "SE", "narrow"]], "PE")
expect_error("chip without peak type rejected", [HEADER, ["a", "treat", "g", "chip", "PE", "none"]], "narrow or broad")
expect_error("atac with explicit peak type rejected", [HEADER, ["a", "treat", "g", "atac", "PE", "narrow"]], "peak_type")
expect_error("mixed peak types within group rejected",
             [HEADER, ["a", "treat", "g", "chip", "PE", "narrow"],
              ["b", "control", "g", "chip", "PE", "broad"]], "same seqtype/peak_type/layout")
expect_error("group without treat rejected", [HEADER, ["b", "control", "g", "chip", "PE", "narrow"]], "treat")
expect_error("same sample across seqtypes rejected",
             [HEADER, ["a", "treat", "g1", "chip", "PE", "narrow"],
              ["a", "treat", "g2", "atac", "PE", "none"]], "different seqtypes")
expect_error("empty table rejected", [HEADER], "no data rows")
expect_error("comma in sample name rejected", [HEADER, ["a,b", "treat", "g", "chip", "PE", "narrow"]], "illegal characters")
expect_error("double underscore in sample name rejected", [HEADER, ["a__b", "treat", "g", "chip", "PE", "narrow"]], "illegal characters")
expect_error("space in group name rejected", [HEADER, ["a", "treat", "g 1", "chip", "PE", "narrow"]], "illegal characters")

print("== 3. ATAC peak-calling mode whitelist ==")
with open(os.path.join(REPO, "workflow", "rules", "common.smk"), encoding="utf-8") as fh:
    wf_src = fh.read()
mode_block = wf_src[wf_src.index('if config["peak"]["atac"]["mode"]'):]
_end = mode_block.find("\ndef validate_config(")
if _end == -1:
    _end = mode_block.find("\n# ---")
mode_block = mode_block[:_end]
for bad in ["BAMPE", "bam", ""]:
    ns = {"config": {"peak": {"atac": {"mode": bad}}}, "WorkflowError": WorkflowError}
    try:
        exec(compile(mode_block, "mode_check", "exec"), ns)
        check(f"mode={bad!r} rejected", False, "WorkflowError was not raised")
    except WorkflowError:
        check(f"mode={bad!r} rejected", True)
for good in ["bampe", "shifted"]:
    ns = {"config": {"peak": {"atac": {"mode": good}}}, "WorkflowError": WorkflowError}
    try:
        exec(compile(mode_block, "mode_check", "exec"), ns)
        check(f"mode={good!r} accepted", True)
    except WorkflowError as e:
        check(f"mode={good!r} accepted", False, str(e))

print("== 4. Centralized config validation (validate_config, extracted from real source) ==")
vc_block = wf_src[wf_src.index("def validate_config("):]
_end = vc_block.find("\nvalidate_config(config)")
if _end == -1:
    _end = vc_block.find("\n# ---")
vc_block = vc_block[:_end]
vc_ns = {"WorkflowError": WorkflowError, "os": os, "BASE_DIR": REPO,
         "ASSAYS": ("chip", "cuttag", "atac", "faire")}
exec(compile(vc_block, "validate_config(extracted)", "exec"), vc_ns)
validate_config = vc_ns["validate_config"]

import copy  # noqa: E402

GOOD_CFG = {
    "genome_fa": "/nonexistent/genome.fa", "gtf": "/nonexistent/genes.gtf",
    "bed": "/nonexistent/genes.bed", "chromsize": "/nonexistent/chrom.sizes",
    "genome_size": "3.7e8", "grouplist": "sample_info.csv",
    "threads": 12, "bowtie2_extra": "--very-sensitive", "min_mapq": 30,
    "region_flank": 3000,
    "dedup": {"chip": True, "cuttag": False, "atac": True, "faire": True},
    "peak": {"keepdup": "all", "qvalue": 0.05, "broad_cutoff": 0.05,
             "atac": {"mode": "bampe", "shift": -100, "extsize": 200}},
    "qc": {"nsc_rsc": False, "frip": True, "deeptools": True},
    "trim": {"quality": 25, "stringency": 3, "error_rate": 0.1, "extra": ""},
}


def vc_case(name, mutate, needles, expect_error=True):
    cfg = copy.deepcopy(GOOD_CFG)
    mutate(cfg)
    try:
        validate_config(cfg)
        check(name, not expect_error, "WorkflowError was not raised as expected")
    except WorkflowError as e:
        if expect_error:
            check(name, all(n in str(e) for n in needles),
                  f"error is missing {needles}: {e}")
        else:
            check(name, False, f"should not fail: {e}")


vc_case("valid complete config accepted (missing reference files warn only)",
        lambda c: None, [], expect_error=False)
vc_case("missing trim key reported", lambda c: c.pop("trim"), ["missing required config key: trim"])
vc_case("non-boolean dedup reported",
        lambda c: c["dedup"].__setitem__("chip", 1), ["dedup.chip"])
vc_case("non-boolean qc switch reported",
        lambda c: c["qc"].__setitem__("frip", "yes"), ["qc.frip"])
vc_case("negative min_mapq reported",
        lambda c: c.__setitem__("min_mapq", -1), ["min_mapq"])
vc_case("missing region_flank key reported",
        lambda c: c.pop("region_flank"), ["missing required config key: region_flank"])
vc_case("non-integer region_flank reported",
        lambda c: c.__setitem__("region_flank", "3000"), ["region_flank"])
vc_case("peak.qvalue out of range reported",
        lambda c: c["peak"].__setitem__("qvalue", 5), ["peak.qvalue"])
vc_case("negative trim.quality reported",
        lambda c: c["trim"].__setitem__("quality", -1), ["trim.quality"])
vc_case("trim.error_rate out of range reported",
        lambda c: c["trim"].__setitem__("error_rate", 2), ["trim.error_rate"])


def _two_errors(c):
    c.pop("gtf")
    c["dedup"].pop("atac")


vc_case("multiple errors aggregated in one report",
        _two_errors, ["2 issues", "gtf", "dedup.atac"])

print("== 5. Wildcard constraint regex ==")
rx = _group_regex(["myc_vs_IgG", "atac.leaf"])
check("regex: exact match (with escaping)",
      re.fullmatch(rx, "myc_vs_IgG") and re.fullmatch(rx, "atac.leaf"))
check("regex: unlisted groups and variants do not match",
      not re.fullmatch(rx, "other") and not re.fullmatch(rx, "atacXleaf"))
check("regex: empty list never matches", re.fullmatch(_group_regex([]), "anything") is None)

print("== 6. mqc shell rule bodies executed for real (snakemake-style format rendering + bash) ==")
import subprocess  # noqa: E402


def render_rule_body(smk_path, rule_name, fmt, literals=None):
    """Extract a rule's shell body and render it, mimicking snakemake's parsing chain:
    (1) the shell body is decoded with Python literal semantics (ast.literal_eval,
    equivalent to snakemake's escape/line-continuation handling of source strings,
    and immune to CRLF checkouts in the working tree);
    (2) named outputs / dotted tokens are replaced literally first (snakemake-specific
    syntax);
    (3) the rest is rendered with str.format semantics ({{}} collapses to {})."""
    with open(smk_path, encoding="utf-8") as fh:
        src = fh.read()
    m = re.search(rf'^rule {rule_name}:.*?shell:\n\s*"""\n(.*?)"""',
                  src, re.S | re.M)
    if not m:
        return None
    raw = m.group(1)
    try:
        import ast
        body = ast.literal_eval('"""' + raw + '"""')
    except (SyntaxError, ValueError):
        body = raw.replace("\r\n", "\n").replace("\r", "\n")
    for token, value in (literals or {}).items():
        body = body.replace("{" + token + "}", value)
    return body.format(**fmt)


def run_bash(body, cwd):
    # bash -c executes the body directly, avoiding Windows path mangling under MSYS bash
    return subprocess.run(["bash", "-c", "set -eo pipefail\n" + body], cwd=cwd,
                          capture_output=True, text=True, timeout=60)


tmp = tempfile.mkdtemp(prefix="mqc_")
os.makedirs(os.path.join(tmp, "5.QC", "frip"), exist_ok=True)
os.makedirs(os.path.join(tmp, "5.QC", "spp"), exist_ok=True)

# --- frip_summary ---
for name, sample, group in [("a", "myc", "myc_vs_IgG"), ("b", "IgG", "myc_vs_IgG")]:
    with open(os.path.join(tmp, "5.QC", "frip", f"{name}.tsv"), "w",
              newline="\n") as fh:
        fh.write("sample\tgroup\ttotal_reads\treads_in_peaks\tFRiP\n")
        fh.write(f"{sample}\t{group}\t1000\t50\t0.0500\n")
body = render_rule_body(
    os.path.join(REPO, "workflow", "rules", "frip.smk"), "frip_summary",
    {"input": "5.QC/frip/a.tsv 5.QC/frip/b.tsv",
     "log": "frip_summary.log"},
    literals={"output.tsv": "5.QC/frip/FRiP_summary.tsv",
              "output.mqc": "5.QC/frip/FRiP_mqc.tsv"})
r = run_bash(body, tmp)
mqc = os.path.join(tmp, "5.QC", "frip", "FRiP_mqc.tsv")
ok = (r.returncode == 0 and os.path.exists(mqc)
      and "# id: 'frip_table'" in open(mqc, encoding="utf-8").read()
      and open(mqc, encoding="utf-8").read().count("sample\tgroup") == 1
      and "myc\t" in open(mqc, encoding="utf-8").read())
check("frip_summary mqc: renders + executes + format correct", ok,
      f"rc={r.returncode} stderr={r.stderr[:200]}")

# --- spp_summary ---
for s in ["s1", "s2"]:
    with open(os.path.join(tmp, "5.QC", "spp", f"{s}_fragment_len.txt"), "w",
              newline="\n") as fh:
        fh.write("150\n")
    with open(os.path.join(tmp, "5.QC", "spp", f"{s}_NSC.txt"), "w",
              newline="\n") as fh:
        fh.write("1.15\n")
    with open(os.path.join(tmp, "5.QC", "spp", f"{s}_RSC.txt"), "w",
              newline="\n") as fh:
        fh.write("0.95\n")
body = render_rule_body(
    os.path.join(REPO, "workflow", "rules", "spp_qc.smk"), "spp_summary",
    {"output": "5.QC/spp/NSC_RSC_mqc.tsv", "log": "spp_summary.log"},
    literals={"params.samples": "s1 s2", "params.spp_dir": "5.QC/spp"})
r = run_bash(body, tmp)
mqc = os.path.join(tmp, "5.QC", "spp", "NSC_RSC_mqc.tsv")
content = open(mqc, encoding="utf-8").read() if os.path.exists(mqc) else ""
ok = (r.returncode == 0
      and "# id: 'nsc_rsc_table'" in content
      and "sample\tfragment_length\tNSC\tRSC" in content
      and "s1\t150\t1.15\t0.95" in content and "s2\t" in content)
check("spp_summary mqc: renders + executes + format correct", ok,
      f"rc={r.returncode} stderr={r.stderr[:200]}")

import shutil  # noqa: E402
shutil.rmtree(tmp, ignore_errors=True)

print("== 7. Config integrity ==")
try:
    import yaml  # noqa: F401
    HAS_YAML = True
except ImportError:
    HAS_YAML = False
    print("  SKIP  pyyaml not installed; skipping YAML parsing checks")

if HAS_YAML:
    with open(os.path.join(REPO, "config", "config.yaml"), encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    # Mimic the Snakefile: unset reference keys fall back to species presets.
    with open(os.path.join(REPO, "config", "species.yaml"), encoding="utf-8") as fh:
        presets = yaml.safe_load(fh) or {}
    _species = cfg.get("species", "osa")
    for _k, _v in (presets.get(_species) or {}).items():
        cfg.setdefault(_k, _v)
    required_top = ["genome_fa", "gtf", "bed", "chromsize", "genome_size",
                    "grouplist", "threads", "bowtie2_extra", "min_mapq",
                    "region_flank", "dedup", "peak", "qc", "trim"]
    check("config: all top-level keys present (incl. species preset merge)",
          all(k in cfg for k in required_top),
          str([k for k in required_top if k not in cfg]))
    try:
        validate_config(cfg)
        check("config: repository default config passes validate_config", True)
    except WorkflowError as e:
        check("config: repository default config passes validate_config", False, str(e))
    check("config: species key present and known",
          cfg.get("species") in ("osa", "hsa"), str(cfg.get("species")))
    check("config: dedup covers all four assays",
          set(cfg["dedup"]) == {"chip", "cuttag", "atac", "faire"})
    check("config: peak sub-keys present",
          {"keepdup", "qvalue", "broad_cutoff", "atac"} <= set(cfg["peak"]))
    check("config: qc switches present",
          {"nsc_rsc", "frip", "deeptools"} <= set(cfg["qc"]))
    # resources.yaml: every rule entry must be a known rule with valid fields.
    with open(os.path.join(REPO, "config", "resources.yaml"), encoding="utf-8") as fh:
        res_cfg = yaml.safe_load(fh) or {}
    check("resources.yaml: top-level resources key", "resources" in res_cfg)
    _bad_res = [k for k, v in (res_cfg.get("resources") or {}).items()
                if not isinstance(v, dict) or not {"threads", "mem_mb", "runtime_min"} <= set(v)]
    check("resources.yaml: every rule entry has threads/mem_mb/runtime_min", not _bad_res, str(_bad_res))

print("== 8. Per-rule resource declarations ==")
# --- resource helpers: extract the real source from common.smk, inject config ---
with open(os.path.join(REPO, "workflow", "rules", "common.smk"), encoding="utf-8") as fh:
    _src = fh.read()
_start = _src.index("RESOURCE_DEFAULTS = {")
_end = _src.find("\n# -----", _start)
res_block = _src[_start:_end if _end != -1 else len(_src)]


def _make_helpers(cfg):
    """Build (rthreads, rmem, rruntime) with an injected config (real source, exec'd)."""
    ns = {"config": cfg, "WorkflowError": WorkflowError}
    exec(compile(res_block, "common.smk(resources extracted)", "exec"), ns)
    return ns["rthreads"], ns["rmem"], ns["rruntime"]


_t, _m, _r = _make_helpers({})
check("rmem falls back to defaults without an override block", _m("frip") == 4096)
check("rruntime falls back to defaults without an override block", _r("frip") == 60)
_t_cov, _m_cov, _r_cov = _make_helpers(
    {"resources": {"bowtie2_mapping": {"mem_mb": 32768}}})
check("rmem override of mem_mb takes effect", _m_cov("bowtie2_mapping") == 32768)
check("rmem override block missing key falls back", _r_cov("bowtie2_mapping") == 240)
check("rmem leaves other rules untouched", _m_cov("frip") == 4096)
check("rthreads honors the legacy global thread cap", _make_helpers({"threads": 4})[0]("bowtie2_index") == 4)
check("rthreads uses the rule default when no cap is set", _t("bowtie2_index") == 8)
try:
    _m("no_such_rule")
    check("unknown rule name raises WorkflowError", False)
except WorkflowError:
    check("unknown rule name raises WorkflowError", True)

# --- Static scan: in every .smk containing rules, the rule count must match the
# --- runtime_sec declaration count (guards against missing declarations)
_mismatch = []
_rules_dir = os.path.join(REPO, "workflow", "rules")
for _fname in sorted(os.listdir(_rules_dir)):
    if not _fname.endswith(".smk"):
        continue
    with open(os.path.join(_rules_dir, _fname), encoding="utf-8") as fh:
        _text = fh.read()
    _n_rules = len(re.findall(r"(?m)^rule \w+:", _text))
    if _n_rules == 0:
        continue  # pure helper files such as common.smk have no rule blocks
    _n_rt = len(re.findall(r"(?m)^\s+runtime_sec=", _text))
    if _n_rules != _n_rt:
        _mismatch.append(f"{_fname}: rule={_n_rules} runtime_sec={_n_rt}")
check("every rule block declares runtime_sec (incl. meta.smk)", not _mismatch, "; ".join(_mismatch))

print("== 9. Synthetic test-data generator ==")
import gzip  # noqa: E402
import hashlib  # noqa: E402

_gen_py = os.path.join(REPO, "tests", "make_testdata.py")
_gen_out1 = tempfile.mkdtemp(prefix="testdata_a_")
_gen_out2 = tempfile.mkdtemp(prefix="testdata_b_")

try:
    # The generator runs as a subprocess (current interpreter); --reads 2000 keeps it
    # fast; outputs go to a tempdir, cleaned up in finally
    _gen = subprocess.run([sys.executable, _gen_py, "--outdir", _gen_out1,
                           "--reads", "2000"],
                          capture_output=True, text=True, timeout=300)
    _fq_count = 0
    _raw_dir = os.path.join(_gen_out1, "1.rawdata")
    if os.path.isdir(_raw_dir):
        _fq_count = len([f for f in os.listdir(_raw_dir) if f.endswith(".fq.gz")])
    check("generator exits 0 with 5 samples x PE = 10 fq.gz files",
          _gen.returncode == 0 and _fq_count == 10,
          f"rc={_gen.returncode} fq.gz={_fq_count} stderr={_gen.stderr[-300:]}")

    _fq1 = os.path.join(_raw_dir, "chip_treat_rep1_1.fq.gz")
    if os.path.exists(_fq1):
        with gzip.open(_fq1, "rt") as fh:
            _lines = [fh.readline().rstrip("\n") for _ in range(4)]
        check("FASTQ four-line structure with 50bp reads",
              _lines[0].startswith("@") and _lines[2].startswith("+")
              and len(_lines[1]) == 50 and len(_lines[3]) == 50,
              str(_lines)[:120])
    else:
        check("FASTQ four-line structure with 50bp reads", False, "chip_treat_rep1_1.fq.gz missing")

    _csv_path = os.path.join(_gen_out1, "samples.csv")
    _rows = open(_csv_path, encoding="utf-8").read().strip().splitlines()
    check("sample table: 1 header + 5 data rows",
          len(_rows) == 6 and _rows[0].startswith("sample_id,role,group"),
          f"rows={len(_rows)} header={_rows[0][:40]!r}")

    check("reference file set complete",
          all(os.path.exists(os.path.join(_gen_out1, "ref", f))
              for f in ("genome.fa", "genes.gtf", "genes.bed", "chrom.sizes")),
          str(os.listdir(os.path.join(_gen_out1, "ref")))
          if os.path.isdir(os.path.join(_gen_out1, "ref")) else "ref/ missing")

    _hdrs, _seqs = [], []
    with open(os.path.join(_gen_out1, "ref", "genome.fa"), encoding="ascii") as fh:
        for line in fh:
            if line.startswith(">"):
                _hdrs.append(line[1:].strip())
                _seqs.append([])
            else:
                _seqs[-1].append(line.strip())
    _seqs = ["".join(s) for s in _seqs]
    check("genome has two 100000bp chromosomes",
          _hdrs == ["chr1", "chr2"] and all(len(s) == 100000 for s in _seqs),
          f"headers={_hdrs} lens={[len(s) for s in _seqs]}")

    # Determinism: generate a second time with the same arguments; genome.fa must be
    # byte-identical (same md5)
    if _gen.returncode == 0:
        _gen2 = subprocess.run([sys.executable, _gen_py, "--outdir", _gen_out2,
                                "--reads", "2000"],
                               capture_output=True, text=True, timeout=300)

        def _md5(path):
            with open(path, "rb") as fh:
                return hashlib.md5(fh.read()).hexdigest()

        _fa1 = os.path.join(_gen_out1, "ref", "genome.fa")
        _fa2 = os.path.join(_gen_out2, "ref", "genome.fa")
        _det_ok = (_gen2.returncode == 0 and os.path.exists(_fa2)
                   and _md5(_fa1) == _md5(_fa2))
        check("determinism: genome.fa md5 identical across two runs with the same args", _det_ok,
              f"rc2={_gen2.returncode} stderr2={_gen2.stderr[-200:]}")
    else:
        check("determinism: genome.fa md5 identical across two runs with the same args", False, "first generation failed; skipped")
except Exception as exc:  # generator-group exceptions must not abort the remaining summary
    check("generator assertion group aborted abnormally", False, repr(exc))
finally:
    shutil.rmtree(_gen_out1, ignore_errors=True)
    shutil.rmtree(_gen_out2, ignore_errors=True)

print()
if FAILED:
    print(f"Result: {len(FAILED)} check(s) failed -> {FAILED}")
    sys.exit(1)
print("Result: all passed")
