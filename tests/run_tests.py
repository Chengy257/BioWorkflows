#!/usr/bin/env python3
"""零依赖测试：样本表解析、路由辅助函数、config/envs 完整性。

运行: python tests/run_tests.py（在仓库根目录）
- 从 workflow.smk 提取【真实源码】执行（非副本），保证测试与实现一致
- pyyaml 可选：缺失时跳过 yaml 检查（CI 中会安装）
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
# 从 workflow.smk 提取被测函数（真实源码）
# ---------------------------------------------------------------------
class WorkflowError(Exception):
    """snakemake.exceptions.WorkflowError 的替身（签名兼容：单 str 参数）"""


def extract_workflow_functions():
    with open(os.path.join(REPO, "workflow.smk"), encoding="utf-8") as fh:
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
    # 截到 _group_regex 函数体结束（下一个顶层 def 或注释分隔）
    m = re.search(r"\n(?=def |# -|BAM_TARGETS)", block_b[len("def _group_regex("):])
    if m:
        block_b = block_b[: len("def _group_regex(") + m.start() + 1]

    ns = {"csv": csv, "os": os, "re": re, "WorkflowError": WorkflowError,
          "REPO_DIR": REPO}
    exec(compile(block_a + "\n\n" + block_b, "workflow.smk(extracted)", "exec"), ns)
    return ns


WF = extract_workflow_functions()
load_sample_table = WF["load_sample_table"]
_resolve_sample_table = WF["_resolve_sample_table"]
_group_regex = WF["_group_regex"]


def write_csv(rows):
    """rows[0] 为表头；返回临时文件路径"""
    fd, path = tempfile.mkstemp(suffix=".csv")
    with os.fdopen(fd, "w", newline="", encoding="utf-8") as fh:
        csv.writer(fh).writerows(rows)
    return path


def expect_error(name, rows, needle):
    path = write_csv(rows)
    try:
        load_sample_table(path)
        check(name, False, "未抛出 WorkflowError")
    except WorkflowError as e:
        check(name, needle in str(e), f"报错信息不含 {needle!r}: {e}")
    finally:
        os.unlink(path)


HEADER = ["sample_id", "role", "group", "seqtype", "layout", "peak_type"]

# ---------------------------------------------------------------------
print("== 1. 样本表解析（真实示例文件） ==")
samples, groups, seqtype_of = load_sample_table(
    os.path.join(REPO, "sample_info.example.csv"))
check("example: 样本数 8 且去重保序",
      samples == ["myc", "IgG", "H3K27me3_rep1", "IgG_cuta",
                  "atac_leaf_1", "atac_leaf_2", "faire_root", "Input_faire"],
      str(samples))
check("example: myc_vs_IgG 结构",
      groups["myc_vs_IgG"]["treat"] == ["myc"]
      and groups["myc_vs_IgG"]["control"] == ["IgG"]
      and groups["myc_vs_IgG"]["seqtype"] == "chip"
      and groups["myc_vs_IgG"]["peak_type"] == "narrow")
check("example: atac_leaf 无对照合法",
      groups["atac_leaf"]["control"] == [] and groups["atac_leaf"]["seqtype"] == "atac")
check("example: faire 组 seqtype 映射", seqtype_of["faire_root"] == "faire")

s2, g2, st2 = load_sample_table(os.path.join(REPO, "sample_info.csv"))
check("repo: sample_info.csv 可解析",
      s2 == ["myc", "IgG"] and g2["myc_vs_IgG"]["peak_type"] == "narrow")

print("== 2. 样本表校验错误路径 ==")
expect_error("缺列报错", HEADER[:4] + [["a", "treat", "g", "chip"]], "缺少列")
expect_error("role 非法", [HEADER, ["a", "TREATMENT", "g", "chip", "PE", "narrow"]], "role")
expect_error("seqtype 非法", [HEADER, ["a", "treat", "g", "chipseq", "PE", "narrow"]], "seqtype")
expect_error("layout SE 拒绝", [HEADER, ["a", "treat", "g", "chip", "SE", "narrow"]], "PE")
expect_error("chip 无峰型拒绝", [HEADER, ["a", "treat", "g", "chip", "PE", "none"]], "narrow 或 broad")
expect_error("atac 指定峰型拒绝", [HEADER, ["a", "treat", "g", "atac", "PE", "narrow"]], "peak_type")
expect_error("组内峰型不一致拒绝",
             [HEADER, ["a", "treat", "g", "chip", "PE", "narrow"],
              ["b", "control", "g", "chip", "PE", "broad"]], "一致")
expect_error("组无 treat 拒绝", [HEADER, ["b", "control", "g", "chip", "PE", "narrow"]], "treat")
expect_error("同样本跨 seqtype 拒绝",
             [HEADER, ["a", "treat", "g1", "chip", "PE", "narrow"],
              ["a", "treat", "g2", "atac", "PE", "none"]], "不同 seqtype")
expect_error("空表拒绝", [HEADER], "没有任何数据行")
expect_error("样本名含逗号拒绝", [HEADER, ["a,b", "treat", "g", "chip", "PE", "narrow"]], "非法字符")
expect_error("样本名含双下划线拒绝", [HEADER, ["a__b", "treat", "g", "chip", "PE", "narrow"]], "非法字符")
expect_error("分组名含空格拒绝", [HEADER, ["a", "treat", "g 1", "chip", "PE", "narrow"]], "非法字符")

print("== 3. ATAC 峰调用模式白名单 ==")
with open(os.path.join(REPO, "workflow.smk"), encoding="utf-8") as fh:
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
        check(f"mode={bad!r} 拒绝", False, "未抛出 WorkflowError")
    except WorkflowError:
        check(f"mode={bad!r} 拒绝", True)
for good in ["bampe", "shifted"]:
    ns = {"config": {"peak": {"atac": {"mode": good}}}, "WorkflowError": WorkflowError}
    try:
        exec(compile(mode_block, "mode_check", "exec"), ns)
        check(f"mode={good!r} 放行", True)
    except WorkflowError as e:
        check(f"mode={good!r} 放行", False, str(e))

print("== 4. config 集中校验（validate_config，真实源码提取） ==")
vc_block = wf_src[wf_src.index("def validate_config("):]
_end = vc_block.find("\nvalidate_config(config)")
if _end == -1:
    _end = vc_block.find("\n# ---")
vc_block = vc_block[:_end]
vc_ns = {"WorkflowError": WorkflowError, "os": os, "REPO_DIR": REPO,
         "ASSAYS": ("chip", "cuttag", "atac", "faire")}
exec(compile(vc_block, "validate_config(extracted)", "exec"), vc_ns)
validate_config = vc_ns["validate_config"]

import copy  # noqa: E402

GOOD_CFG = {
    "genome_fa": "/nonexistent/genome.fa", "gtf": "/nonexistent/genes.gtf",
    "bed": "/nonexistent/genes.bed", "chromsize": "/nonexistent/chrom.sizes",
    "genome_size": "3.7e8", "grouplist": "sample_info.csv",
    "threads": 12, "bowtie2_extra": "--very-sensitive", "min_mapq": 30,
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
        check(name, not expect_error, "未按预期抛出 WorkflowError")
    except WorkflowError as e:
        if expect_error:
            check(name, all(n in str(e) for n in needles),
                  f"报错缺少 {needles}: {e}")
        else:
            check(name, False, f"不应报错: {e}")


vc_case("合法完整 config 放行（参考文件缺失仅警告）",
        lambda c: None, [], expect_error=False)
vc_case("缺 trim 键报错", lambda c: c.pop("trim"), ["缺少必需配置键: trim"])
vc_case("dedup 非布尔报错",
        lambda c: c["dedup"].__setitem__("chip", 1), ["dedup.chip"])
vc_case("qc 开关非布尔报错",
        lambda c: c["qc"].__setitem__("frip", "yes"), ["qc.frip"])
vc_case("min_mapq 负数报错",
        lambda c: c.__setitem__("min_mapq", -1), ["min_mapq"])
vc_case("peak.qvalue 超区间报错",
        lambda c: c["peak"].__setitem__("qvalue", 5), ["peak.qvalue"])
vc_case("trim.quality 负数报错",
        lambda c: c["trim"].__setitem__("quality", -1), ["trim.quality"])
vc_case("trim.error_rate 超区间报错",
        lambda c: c["trim"].__setitem__("error_rate", 2), ["trim.error_rate"])


def _two_errors(c):
    c.pop("gtf")
    c["dedup"].pop("atac")


vc_case("多错误一次汇总报出",
        _two_errors, ["共 2 项", "gtf", "dedup.atac"])

print("== 5. 通配符约束正则 ==")
rx = _group_regex(["myc_vs_IgG", "atac.leaf"])
check("regex: 精确匹配（含转义）",
      re.fullmatch(rx, "myc_vs_IgG") and re.fullmatch(rx, "atac.leaf"))
check("regex: 不匹配未列分组与变形",
      not re.fullmatch(rx, "other") and not re.fullmatch(rx, "atacXleaf"))
check("regex: 空列表永不匹配", re.fullmatch(_group_regex([]), "anything") is None)

print("== 6. config 与 envs 完整性 ==")
try:
    import yaml  # noqa: F401
    HAS_YAML = True
except ImportError:
    HAS_YAML = False
    print("  SKIP  pyyaml 未安装，跳过 yaml 解析检查")

if HAS_YAML:
    with open(os.path.join(REPO, "config", "config.yaml"), encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    required_top = ["genome_fa", "gtf", "bed", "chromsize", "genome_size",
                    "grouplist", "threads", "bowtie2_extra", "min_mapq",
                    "dedup", "peak", "qc", "trim"]
    check("config: 顶层键齐全", all(k in cfg for k in required_top),
          str([k for k in required_top if k not in cfg]))
    try:
        validate_config(cfg)
        check("config: 仓库默认 config 通过 validate_config", True)
    except WorkflowError as e:
        check("config: 仓库默认 config 通过 validate_config", False, str(e))
    check("config: dedup 四 assay 齐全",
          set(cfg["dedup"]) == {"chip", "cuttag", "atac", "faire"})
    check("config: peak 子键齐全",
          {"keepdup", "qvalue", "broad_cutoff", "atac"} <= set(cfg["peak"]))
    check("config: qc 开关齐全",
          {"nsc_rsc", "frip", "deeptools"} <= set(cfg["qc"]))

    env_dir = os.path.join(REPO, "envs")
    env_files = sorted(f for f in os.listdir(env_dir) if f.endswith(".yaml"))
    ok = True
    for f in env_files:
        with open(os.path.join(env_dir, f), encoding="utf-8") as fh:
            doc = yaml.safe_load(fh)
        if "dependencies" not in doc:
            ok = False
            print(f"  FAIL  envs/{f} 缺 dependencies")
    check(f"envs: {len(env_files)} 个环境文件均可解析且含 dependencies", ok)

print()
if FAILED:
    print(f"结果: {len(FAILED)} 项失败 -> {FAILED}")
    sys.exit(1)
print("结果: 全部通过")
