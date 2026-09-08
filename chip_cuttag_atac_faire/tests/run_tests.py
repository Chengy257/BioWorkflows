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
    # Replicate/QC pure helpers (v0.5): from _unordered_pairs to the next divider
    block_c = segment("def _unordered_pairs(", ["\n# -----"])

    ns = {"csv": csv, "os": os, "re": re, "WorkflowError": WorkflowError,
          "BASE_DIR": REPO}
    exec(compile(block_a + "\n\n" + block_b + "\n\n" + block_c,
                 "common.smk(extracted)", "exec"), ns)
    return ns


WF = extract_workflow_functions()
load_sample_table = WF["load_sample_table"]
_resolve_sample_table = WF["_resolve_sample_table"]
_group_regex = WF["_group_regex"]
_unordered_pairs = WF["_unordered_pairs"]
idr_pair_slug = WF["idr_pair_slug"]
parse_idr_pair_slug = WF["parse_idr_pair_slug"]
_group_calls_narrow = WF["_group_calls_narrow"]
_is_organelle_contig = WF["_is_organelle_contig"]


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
samples, groups, seqtype_of, cond_of, batch_of = load_sample_table(
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

s2, g2, st2, _c2, _b2 = load_sample_table(os.path.join(REPO, "example", "samples.csv"))
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

print("== 2b. Optional condition/batch columns (differential binding) ==")
EXT = HEADER + ["condition", "batch"]
_ext_path = write_csv([EXT,
                       ["a", "treat", "g1", "chip", "PE", "narrow", "WT", "b1"],
                       ["b", "treat", "g1", "chip", "PE", "narrow", "WT", "b2"],
                       ["ctl", "control", "g1", "chip", "PE", "narrow", "", ""],
                       ["c", "treat", "g2", "chip", "PE", "narrow", "mut", "b1"],
                       ["d", "treat", "g2", "chip", "PE", "narrow", "mut", "b2"]])
try:
    _s, _g, _st, _cond, _batch = load_sample_table(_ext_path)
    check("extended table parses; condition/batch captured for treats",
          _cond == {"a": "WT", "b": "WT", "c": "mut", "d": "mut"}
          and _batch["a"] == "b1" and _batch["d"] == "b2"
          and _g["g1"]["condition"] == "WT")
    check("extended table: legacy 6-column tables unaffected (empty maps)",
          all(v == "" for v in load_sample_table(
              os.path.join(REPO, "config", "samples.csv"))[3].values()))
finally:
    os.unlink(_ext_path)
expect_error("inconsistent condition within group rejected",
             [EXT, ["a", "treat", "g", "chip", "PE", "narrow", "WT", ""],
              ["b", "treat", "g", "chip", "PE", "narrow", "mut", ""]],
             "different condition values")
expect_error("illegal condition value rejected",
             [EXT, ["a", "treat", "g", "chip", "PE", "narrow", "W T", ""]],
             "illegal characters")
expect_error("illegal batch value rejected",
             [EXT, ["a", "treat", "g", "chip", "PE", "narrow", "WT", "b;1"]],
             "illegal characters")

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
             "bigwig_measure": "FE",
             "atac": {"mode": "bampe", "shift": -100, "extsize": 200},
             "replicate": {"enabled": False, "qvalue": 0.01,
                           "idr_threshold": 0.05, "idr_rank": "p.value",
                           "consensus_min_replicates": 2, "frip_on": "pooled"}},
    "blacklist": "",
    "bigwig": {"per_sample": False, "normalize": "RPGC", "bin": 25},
    "motif": {"enabled": False, "homer_genome": "", "size": "given",
              "background": "", "extra": ""},
    "diffbind": {"enabled": False, "contrasts": [], "analysis": "DESeq2",
                 "summit_flank": 250, "use_controls": False,
                 "fdr": 0.05, "foldchange": 1.0, "batch_correction": True},
    "qc": {"nsc_rsc": False, "frip": True, "deeptools": True,
           "tss": False, "organelle": False,
           "organelle_patterns": ["chrc", "chrm"]},
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

# --- v0.5 keys: peak.replicate / blacklist / qc.tss / qc.organelle ---
vc_case("replicate block is optional (legacy configs stay valid)",
        lambda c: c["peak"].pop("replicate"), [], expect_error=False)
vc_case("bad idr_rank reported",
        lambda c: c["peak"]["replicate"].__setitem__("idr_rank", "pvalue"),
        ["idr_rank"])
vc_case("idr_threshold out of range reported",
        lambda c: c["peak"]["replicate"].__setitem__("idr_threshold", 5),
        ["peak.replicate.idr_threshold"])
vc_case("replicate qvalue non-numeric reported",
        lambda c: c["peak"]["replicate"].__setitem__("qvalue", "x"),
        ["peak.replicate.qvalue"])
vc_case("consensus_min_replicates below 2 reported",
        lambda c: c["peak"]["replicate"].__setitem__("consensus_min_replicates", 1),
        ["consensus_min_replicates"])
vc_case("bad frip_on reported",
        lambda c: c["peak"]["replicate"].__setitem__("frip_on", "idr"),
        ["frip_on"])
vc_case("frip_on=consensus without enabled reported",
        lambda c: c["peak"]["replicate"].__setitem__("frip_on", "consensus"),
        ["frip_on=consensus requires peak.replicate.enabled"])
vc_case("frip_on=consensus with enabled accepted",
        lambda c: (c["peak"]["replicate"].__setitem__("enabled", True),
                   c["peak"]["replicate"].__setitem__("frip_on", "consensus")),
        [], expect_error=False)
vc_case("replicate block of wrong type reported",
        lambda c: c["peak"].__setitem__("replicate", "on"),
        ["peak.replicate must be a mapping"])
vc_case("non-boolean qc.tss reported",
        lambda c: c["qc"].__setitem__("tss", "yes"), ["qc.tss"])
vc_case("organelle_patterns of wrong shape reported",
        lambda c: c["qc"].__setitem__("organelle_patterns", "chrc"),
        ["organelle_patterns"])
vc_case("non-string blacklist reported",
        lambda c: c.__setitem__("blacklist", 3), ["blacklist must be a string"])
# --- v0.5 phase 3/4/5 keys: bigwig / motif / diffbind ---
vc_case("bigwig/motif/diffbind blocks are optional (legacy configs stay valid)",
        lambda c: (c.pop("bigwig"), c.pop("motif"), c.pop("diffbind")),
        [], expect_error=False)
vc_case("bad bigwig.normalize reported",
        lambda c: c["bigwig"].__setitem__("normalize", "RPKM"),
        ["bigwig.normalize"])
vc_case("non-positive bigwig.bin reported",
        lambda c: c["bigwig"].__setitem__("bin", 0), ["bigwig.bin"])
vc_case("non-boolean bigwig.per_sample reported",
        lambda c: c["bigwig"].__setitem__("per_sample", "yes"), ["bigwig.per_sample"])
vc_case("bad peak.bigwig_measure reported",
        lambda c: c["peak"].__setitem__("bigwig_measure", "log2"),
        ["peak.bigwig_measure"])
vc_case("motif enabled without homer_genome reported",
        lambda c: c["motif"].__setitem__("enabled", True),
        ["motif.homer_genome is required"])
vc_case("motif enabled with homer_genome accepted",
        lambda c: (c["motif"].__setitem__("enabled", True),
                   c["motif"].__setitem__("homer_genome", "hg38")),
        [], expect_error=False)
vc_case("bad diffbind.analysis reported",
        lambda c: c["diffbind"].__setitem__("analysis", "deseq"),
        ["diffbind.analysis"])
vc_case("negative diffbind.summit_flank reported",
        lambda c: c["diffbind"].__setitem__("summit_flank", -1),
        ["diffbind.summit_flank"])
vc_case("diffbind.fdr out of range reported",
        lambda c: c["diffbind"].__setitem__("fdr", 2), ["diffbind.fdr"])
vc_case("diffbind.foldchange below 1 reported",
        lambda c: c["diffbind"].__setitem__("foldchange", 0.5), ["diffbind.foldchange"])
vc_case("diffbind block of wrong type reported",
        lambda c: c.__setitem__("diffbind", "on"), ["diffbind must be a mapping"])

print("== 5. Wildcard constraint regex ==")
rx = _group_regex(["myc_vs_IgG", "atac.leaf"])
check("regex: exact match (with escaping)",
      re.fullmatch(rx, "myc_vs_IgG") and re.fullmatch(rx, "atac.leaf"))
check("regex: unlisted groups and variants do not match",
      not re.fullmatch(rx, "other") and not re.fullmatch(rx, "atacXleaf"))
check("regex: empty list never matches", re.fullmatch(_group_regex([]), "anything") is None)

print("== 5b. Replicate/QC helper functions (real source) ==")
check("pairs: unordered pair enumeration order and count",
      _unordered_pairs(["a", "b", "c"]) == [("a", "b"), ("a", "c"), ("b", "c")])
check("pairs: empty and singleton lists yield no pair",
      _unordered_pairs([]) == [] and _unordered_pairs(["a"]) == [])
slug = idr_pair_slug("rep1", "rep2")
check("idr slug: round-trips through the parser",
      slug == "rep1__vs__rep2" and parse_idr_pair_slug(slug) == ("rep1", "rep2"))
check("group routing: atac/faire always narrow",
      _group_calls_narrow({"seqtype": "atac", "peak_type": "none"})
      and _group_calls_narrow({"seqtype": "faire", "peak_type": "none"}))
check("group routing: chip by peak_type",
      _group_calls_narrow({"seqtype": "chip", "peak_type": "narrow"})
      and not _group_calls_narrow({"seqtype": "chip", "peak_type": "broad"}))
check("organelle: short patterns match exactly, not as substrings",
      _is_organelle_contig("ChrC", ["chrc", "chrm"])
      and _is_organelle_contig("chrM", ["chrc", "chrm"])
      and not _is_organelle_contig("chrUn_ptg0001l", ["pt", "mt"]))
check("organelle: long patterns match as substrings",
      _is_organelle_contig("mitochondrion_genome", ["mitochondr"])
      and _is_organelle_contig("Oschloroplast_fake", ["chloroplast"]))
check("organelle: no match on plain chromosomes",
      not _is_organelle_contig("chr1", ["chrc", "chrm", "pt", "mt"]))

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
          {"keepdup", "qvalue", "broad_cutoff", "atac", "replicate"} <= set(cfg["peak"]))
    check("config: replicate sub-keys present",
          {"enabled", "qvalue", "idr_threshold", "idr_rank",
           "consensus_min_replicates", "frip_on"} <= set(cfg["peak"]["replicate"]))
    check("config: blacklist key present and defaults to disabled",
          cfg.get("blacklist", None) == "")
    check("config: qc switches present",
          {"nsc_rsc", "frip", "deeptools", "tss", "organelle"} <= set(cfg["qc"]))
    check("config: new QC switches default to off",
          cfg["qc"]["tss"] is False and cfg["qc"]["organelle"] is False
          and cfg["peak"]["replicate"]["enabled"] is False)
    check("config: v0.5 phase 3/4/5 keys present and default to off",
          cfg.get("bigwig", {}).get("per_sample") is False
          and cfg.get("motif", {}).get("enabled") is False
          and cfg.get("diffbind", {}).get("enabled") is False
          and cfg["peak"].get("bigwig_measure", "FE") == "FE")
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

print("== 10. v0.5 stage scripts (tss_from_bed / tss_score / organelle_summary / replicate_summary) ==")
_stmp = tempfile.mkdtemp(prefix="v05_scripts_")


def _run_py(script, args, cwd=None):
    return subprocess.run([sys.executable, os.path.join(REPO, "workflow", "scripts", script)] + args,
                          capture_output=True, text=True, timeout=60, cwd=cwd)


# --- tss_from_bed.py: strand-aware derivation, comment skipping, BED4 error ---
_bed = os.path.join(_stmp, "genes.bed")
with open(_bed, "w", newline="\n", encoding="utf-8") as fh:
    fh.write("# comment line skipped\n")
    fh.write('track name="also skipped"\n')
    fh.write("chr1\t1000\t2000\tgene1\t0\t+\n")
    fh.write("chr2\t3000\t4000\tgene2\t0\t-\n")
_tss = os.path.join(_stmp, "tss.bed")
_r = _run_py("tss_from_bed.py", [_bed, _tss])
_tss_lines = open(_tss, encoding="utf-8").read().splitlines() if os.path.exists(_tss) else []
check("tss_from_bed: exits 0, skips comments, 1-bp TSS per strand",
      _r.returncode == 0
      and _tss_lines == ["chr1\t1000\t1001", "chr2\t3999\t4000"],
      f"rc={_r.returncode} lines={_tss_lines}")
_bed4 = os.path.join(_stmp, "bed4.bed")
with open(_bed4, "w", newline="\n", encoding="utf-8") as fh:
    fh.write("chr1\t1000\t2000\tgene1\n")
_r = _run_py("tss_from_bed.py", [_bed4, os.path.join(_stmp, "x.bed")])
check("tss_from_bed: BED4 input rejected with a clear error",
      _r.returncode != 0 and "BED6 required" in _r.stderr, f"rc={_r.returncode}")

# --- tss_score.py: profile aggregation, baseline normalization, max = TSSE ---
# 40 bins: outer 10 bins on each side are the baseline (all 1.0); the bin at
# index 20 spikes to 5.0 (bin 19 stays 1.0) -> TSSE = 5.0, center = mean(1,5) = 3.0
_vals = [1.0] * 40
_vals[20] = 5.0
_mat = os.path.join(_stmp, "matrix.txt")
with open(_mat, "w", newline="\n", encoding="utf-8") as fh:
    fh.write("@" + '{"upstream": [2000], "downstream": [2000]}\n')
    for _ in range(2):  # two identical TSS rows aggregate to the same profile
        fh.write("chr1\t100\t101\tg\t0\t+\t" + "\t".join(str(v) for v in _vals) + "\n")
_sc = os.path.join(_stmp, "score.tsv")
_r = _run_py("tss_score.py", [_mat, "sampleA", _sc])
_sc_line = open(_sc, encoding="utf-8").read().strip() if os.path.exists(_sc) else ""
check("tss_score: TSSE is the normalized profile max, center the TSS bins",
      _r.returncode == 0 and _sc_line == "sampleA\t5.0000\t3.0000\t1.000000\t2",
      f"rc={_r.returncode} line={_sc_line!r}")
_zeromat = os.path.join(_stmp, "zero.txt")
with open(_zeromat, "w", newline="\n", encoding="utf-8") as fh:
    fh.write("@{}\n")
    fh.write("chr1\t100\t101\tg\t0\t+\t" + "\t".join("0" for _ in range(40)) + "\n")
_r = _run_py("tss_score.py", [_zeromat, "sampleB", _sc])
check("tss_score: zero baseline yields NA instead of a division by zero",
      _r.returncode == 0 and open(_sc, encoding="utf-8").read().startswith("sampleB\tNA\tNA\t"),
      f"rc={_r.returncode}")

# --- organelle_summary.py: fractions, '*' row ignored, pattern matching ---
for _s, _rows in [("s1", [("chr1", 100000, 800, 5), ("ChrC", 120000, 150, 1),
                          ("ChrM", 90000, 50, 0), ("*", 0, 0, 40)]),
                  ("s2", [("chr1", 100000, 900, 4), ("chr2", 100000, 100, 2)])]:
    with open(os.path.join(_stmp, f"{_s}_idxstats.tsv"), "w",
              newline="\n", encoding="utf-8") as fh:
        for row in _rows:
            fh.write("\t".join(str(x) for x in row) + "\n")
_os = os.path.join(_stmp, "org.tsv")
_om = os.path.join(_stmp, "org_mqc.tsv")
_r = _run_py("organelle_summary.py",
             ["--patterns", "chrc,chrm", "--out", _os, "--mqc", _om,
              os.path.join(_stmp, "s1_idxstats.tsv"),
              os.path.join(_stmp, "s2_idxstats.tsv")])
_org = open(_os, encoding="utf-8").read().splitlines() if os.path.exists(_os) else []
check("organelle: fractions over mapped reads, '*' ignored, per-contig detail",
      _r.returncode == 0 and len(_org) == 3
      and _org[1].startswith("s1\t1000\t200\t0.2000\t")
      and "ChrC:150" in _org[1] and "ChrM:50" in _org[1]
      and _org[2].startswith("s2\t1000\t0\t0.0000\t-"),
      f"rc={_r.returncode} rows={_org}")
check("organelle: mqc wrapper carries the MultiQC header",
      os.path.exists(_om) and "# id: 'organelle_table'" in open(_om, encoding="utf-8").read())

# --- replicate_summary.py: modes, counts, retained fraction ---
for _rel in ("4.peak/replicates/g1", "4.peak"):
    os.makedirs(os.path.join(_stmp, "results", _rel), exist_ok=True)


def _w(path, n):
    with open(path, "w", newline="\n", encoding="utf-8") as fh:
        fh.writelines(f"chr1\t{i}\t{i + 10}\n" for i in range(n))


_w(os.path.join(_stmp, "results", "4.peak", "replicates", "g1", "a_peaks.narrowPeak"), 10)
_w(os.path.join(_stmp, "results", "4.peak", "replicates", "g1", "b_peaks.narrowPeak"), 20)
_w(os.path.join(_stmp, "results", "4.peak", "g1_IDR_peaks.narrowPeak"), 5)
os.makedirs(os.path.join(_stmp, "results", "4.peak", "replicates", "g2"), exist_ok=True)
_w(os.path.join(_stmp, "results", "4.peak", "replicates", "g2", "h_peaks.broadPeak"), 8)
_w(os.path.join(_stmp, "results", "4.peak", "g2_peaks.broadPeak"), 8)
_rep_samples = os.path.join(_stmp, "samples.csv")
with open(_rep_samples, "w", newline="\n", encoding="utf-8") as fh:
    fh.write("sample_id,role,group,seqtype,layout,peak_type\n")
    fh.write("a,treat,g1,chip,PE,narrow\n")
    fh.write("b,treat,g1,chip,PE,narrow\n")
    fh.write("ctl,control,g1,chip,PE,narrow\n")
    fh.write("h,treat,g2,chip,PE,broad\n")
_rs = os.path.join(_stmp, "rep.tsv")
_rm = os.path.join(_stmp, "rep_mqc.tsv")
_r = _run_py("replicate_summary.py",
             ["--samples", _rep_samples, "--results-dir", os.path.join(_stmp, "results"),
              "--out", _rs, "--mqc", _rm])
_rep = open(_rs, encoding="utf-8").read().splitlines() if os.path.exists(_rs) else []
check("replicate_summary: idr mode rows with counts and retained fraction",
      _r.returncode == 0 and len(_rep) == 3
      and _rep[1] == "g1\tchip\tnarrow\t2\tidr\ta=10;b=20\t5\t0.3333"
      and _rep[2] == "g2\tchip\tbroad\t1\tpooled\th=8\t8\t1.0000",
      f"rc={_r.returncode} rows={_rep}")
check("replicate_summary: mqc wrapper carries the MultiQC header",
      os.path.exists(_rm) and "# id: 'replicate_summary_table'" in open(_rm, encoding="utf-8").read())

# --- diffbind_sheet.py: sheet layout from the extended sample table ---
_dbtmp = os.path.join(_stmp, "db")
os.makedirs(os.path.join(_dbtmp, "results", "4.peak", "replicates", "g1"), exist_ok=True)
os.makedirs(os.path.join(_dbtmp, "results", "4.peak", "replicates", "g4"), exist_ok=True)
for _rel in ("results/3.align/bowtie2",):
    os.makedirs(os.path.join(_dbtmp, _rel), exist_ok=True)


def _touch(path, line="chr1\t1\t10\n"):
    with open(path, "w", newline="\n", encoding="utf-8") as fh:
        fh.write(line)


for _s in ("a", "b", "c", "d"):
    _touch(os.path.join(_dbtmp, "results", "3.align", "bowtie2", f"{_s}_sorted.bam"), "")
    _touch(os.path.join(_dbtmp, "results", "4.peak", "replicates", "g1" if _s in ("a", "b") else "g4",
                        f"{_s}_peaks.narrowPeak"))
_db_samples = os.path.join(_dbtmp, "samples.csv")
with open(_db_samples, "w", newline="\n", encoding="utf-8") as fh:
    fh.write("sample_id,role,group,seqtype,layout,peak_type,condition,batch\n")
    fh.write("a,treat,g1,chip,PE,narrow,WT,b1\n")
    fh.write("b,treat,g1,chip,PE,narrow,WT,b2\n")
    fh.write("ctl,control,g1,chip,PE,narrow,,\n")
    fh.write("c,treat,g4,chip,PE,narrow,mut,b1\n")
    fh.write("d,treat,g4,chip,PE,narrow,mut,b2\n")
_db_sheet = os.path.join(_dbtmp, "sheet.tsv")
_r = _run_py("diffbind_sheet.py",
             ["--samples", _db_samples, "--results-dir",
              os.path.join(_dbtmp, "results"), "--groups", "g1", "g4",
              "--replicates", "--out", _db_sheet])
_sheet = open(_db_sheet, encoding="utf-8").read().splitlines() if os.path.exists(_db_sheet) else []
check("diffbind_sheet: conditions, replicates, per-replicate peaks, batches",
      _r.returncode == 0 and len(_sheet) == 5
      and _sheet[0] == "SampleID\tCondition\tReplicate\tbamReads\tbamControl\tPeakFile\tBatch"
      and _sheet[1].startswith("a\tWT\t1\t") and "replicates/g1/a_peaks.narrowPeak" in _sheet[1]
      and _sheet[1].endswith("\tb1") and _sheet[4].startswith("d\tmut\t2\t"),
      f"rc={_r.returncode} rows={_sheet[:3]}")
_r = _run_py("diffbind_sheet.py",
             ["--samples", _db_samples, "--results-dir",
              os.path.join(_dbtmp, "results"), "--groups", "nope", "g1",
              "--out", _db_sheet])
check("diffbind_sheet: unknown contrast group fails with a clear error",
      _r.returncode != 0 and "not found" in _r.stderr, f"rc={_r.returncode}")
shutil.rmtree(_stmp, ignore_errors=True)

print()
if FAILED:
    print(f"Result: {len(FAILED)} check(s) failed -> {FAILED}")
    sys.exit(1)
print("Result: all passed")
