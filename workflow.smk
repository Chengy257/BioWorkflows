# =====================================================================
# chip_cuttag_atac_faire —— 统一入口 Snakefile
#
# ChIP-seq / CUT&Tag / ATAC-seq / FAIRE-seq 一体化 Snakemake 流程
# 根据 sample_info.csv 的 seqtype 列自动路由，支持混型项目。
#
# 用法（在数据工作目录运行）：
#   snakemake -s /path/to/repo/workflow.smk --configfile /path/to/repo/config/config.yaml \
#       --use-conda --cores 18 -j 3 -k
# 详见 README.md
# =====================================================================

import csv
import os
import re

from snakemake.exceptions import WorkflowError

configfile: os.path.join(workflow.basedir, "config", "config.yaml")

# 所有 shell 以 bash -eo pipefail 执行：任何一步失败立即中断，管道错误可捕获。
# 不加 set -u：conda 激活脚本在 -u 下会因未绑定变量报错。
shell.executable("/bin/bash")
shell.prefix("set -eo pipefail; ")

REPO_DIR = workflow.basedir
ENVS = os.path.join(REPO_DIR, "envs")

ASSAYS = ("chip", "cuttag", "atac", "faire")

# 样本名/分组名仅允许字母数字._-：逗号会破坏峰列表拼接与 MACS2 多文件参数，
# "__" 是 FRiP 输出的分隔符，空格/制表符会破坏 shell 展开。
_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")

# ---------------------------------------------------------------------
# 样本表解析与校验
# ---------------------------------------------------------------------
REQUIRED_COLUMNS = ("sample_id", "role", "group", "seqtype", "layout", "peak_type")


def _resolve_sample_table(path):
    """样本表路径：绝对路径 > 工作目录相对路径 > 仓库目录相对路径"""
    p = os.path.expanduser(str(path))
    if os.path.isabs(p):
        return p
    if os.path.exists(p):
        return os.path.abspath(p)
    return os.path.join(REPO_DIR, p)


def load_sample_table(path):
    """解析样本表，返回 (样本 id 列表, 分组 dict, 样本->seqtype 映射)。

    分组 dict 结构：group -> {seqtype, peak_type, layout, treat: [], control: []}
    """
    samples = []          # 去重后的样本 id（对照组同样需要比对）
    seqtype_of = {}
    groups = {}

    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        missing = [c for c in REQUIRED_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise WorkflowError(
                f"样本表 {path} 缺少列: {missing}；"
                f"必须包含 {list(REQUIRED_COLUMNS)}，参见 sample_info.example.csv"
            )
        for lineno, row in enumerate(reader, start=2):
            sid = (row["sample_id"] or "").strip()
            role = (row["role"] or "").strip().lower()
            grp = (row["group"] or "").strip()
            seqtype = (row["seqtype"] or "").strip().lower()
            layout = (row["layout"] or "").strip().upper()
            peak_type = (row["peak_type"] or "").strip().lower()

            if not sid or not grp:
                raise WorkflowError(f"样本表第 {lineno} 行：sample_id/group 不能为空")
            for label, value in (("sample_id", sid), ("group", grp)):
                if not _NAME_RE.match(value) or "__" in value:
                    raise WorkflowError(
                        f"样本表第 {lineno} 行：{label}={value!r} 含非法字符，"
                        "仅允许字母数字与 . _ -（不以 - 开头，且不含连续下划线 __）"
                    )
            if role not in ("treat", "control"):
                raise WorkflowError(f"样本表第 {lineno} 行：role 必须是 treat 或 control，当前为 {role!r}")
            if seqtype not in ASSAYS:
                raise WorkflowError(f"样本表第 {lineno} 行：seqtype 必须是 {'/'.join(ASSAYS)}，当前为 {seqtype!r}")
            if layout != "PE":
                raise WorkflowError(
                    f"样本表第 {lineno} 行：当前版本仅支持 PE（双端）数据，layout={layout!r}"
                )
            if peak_type not in ("narrow", "broad", "none"):
                raise WorkflowError(
                    f"样本表第 {lineno} 行：peak_type 必须是 narrow/broad/none，当前为 {peak_type!r}"
                )
            if seqtype in ("atac", "faire") and peak_type != "none":
                raise WorkflowError(
                    f"样本表第 {lineno} 行：atac/faire 的 peak_type 请填 none（峰型由流程固定）"
                )
            if sid in seqtype_of and seqtype_of[sid] != seqtype:
                raise WorkflowError(f"样本表第 {lineno} 行：样本 {sid} 出现在不同 seqtype 的分组中")
            if role == "treat" and peak_type == "none" and seqtype in ("chip", "cuttag"):
                raise WorkflowError(
                    f"样本表第 {lineno} 行：chip/cuttag 处理组样本 {sid} 必须指定 narrow 或 broad"
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
                    f"样本表第 {lineno} 行：分组 {grp} 内各行 seqtype/peak_type/layout 必须一致"
                )
            g[role].append(sid)

    bad_groups = [g for g, v in groups.items() if not v["treat"]]
    if bad_groups:
        raise WorkflowError(f"以下分组缺少 role=treat 的处理样本: {bad_groups}")
    if not samples:
        raise WorkflowError(f"样本表 {path} 没有任何数据行")

    return samples, groups, seqtype_of


SAMPLES, GROUPS, SEQTYPE_OF = load_sample_table(_resolve_sample_table(config["grouplist"]))

try:
    config["threads"] = int(config["threads"])
except (TypeError, ValueError):
    raise WorkflowError(
        f"config threads 必须是整数（不要加引号），当前为 {config['threads']!r}"
    )

# 峰调用模式白名单：mode 只接受 bampe / shifted，笔误静默落入 shifted 的风险需在解析期拦截
if config["peak"]["atac"]["mode"] not in ("bampe", "shifted"):
    raise WorkflowError(
        "config peak.atac.mode 必须是 bampe 或 shifted，当前为 "
        f"{config['peak']['atac']['mode']!r}"
    )


def validate_config(cfg):
    """集中校验 config：必需键、子键、类型与取值范围，一次汇总报出全部问题；
    参考文件存在性只警告不中断（--lint/dry-run 时参考文件常不在本机，硬校验会误伤）。"""
    errors, warnings = [], []
    required = ("genome_fa", "gtf", "bed", "chromsize", "genome_size",
                "grouplist", "threads", "bowtie2_extra", "min_mapq",
                "dedup", "peak", "qc", "trim")
    for key in required:
        if key not in cfg:
            errors.append(f"缺少必需配置键: {key}")
    for key in ("dedup", "qc", "peak", "trim"):
        if key in cfg and not isinstance(cfg[key], dict):
            errors.append(f"{key} 必须是映射（含子键），当前为 {cfg[key]!r}")
    if isinstance(cfg.get("dedup"), dict):
        for assay in ASSAYS:
            v = cfg["dedup"].get(assay)
            if not isinstance(v, bool):
                errors.append(f"dedup.{assay} 必须是 true/false，当前为 {v!r}")
    if isinstance(cfg.get("qc"), dict):
        for key in ("nsc_rsc", "frip", "deeptools"):
            v = cfg["qc"].get(key)
            if not isinstance(v, bool):
                errors.append(f"qc.{key} 必须是 true/false，当前为 {v!r}")
    for key in ("min_mapq",):
        v = cfg.get(key)
        if isinstance(v, bool) or not isinstance(v, int) or v < 0:
            errors.append(f"{key} 必须是 >= 0 的整数，当前为 {v!r}")
    if isinstance(cfg.get("peak"), dict):
        p = cfg["peak"]
        for key in ("qvalue", "broad_cutoff"):
            try:
                if not 0 < float(p[key]) <= 1:
                    errors.append(f"peak.{key} 必须在 (0, 1] 区间，当前为 {p[key]!r}")
            except (TypeError, ValueError):
                errors.append(f"peak.{key} 必须是数值，当前为 {p[key]!r}")
    if isinstance(cfg.get("trim"), dict):
        t = cfg["trim"]
        for key, lo in (("quality", 0), ("stringency", 1)):
            v = t.get(key)
            if isinstance(v, bool) or not isinstance(v, int) or v < lo:
                errors.append(f"trim.{key} 必须是 >= {lo} 的整数，当前为 {v!r}")
        try:
            er = float(t.get("error_rate"))
            if not 0 < er <= 1:
                errors.append(f"trim.error_rate 必须在 (0, 1] 区间，当前为 {er!r}")
        except (TypeError, ValueError):
            errors.append(f"trim.error_rate 必须是数值，当前为 {t.get('error_rate')!r}")
        if not isinstance(t.get("extra", ""), str):
            errors.append("trim.extra 必须是字符串")
    if errors:
        raise WorkflowError(
            f"config 校验失败（共 {len(errors)} 项）:\n  " + "\n  ".join(errors)
        )
    for key in ("genome_fa", "gtf", "bed", "chromsize"):
        p = str(cfg[key])
        if not os.path.isabs(p) and not os.path.exists(p):
            p = os.path.join(REPO_DIR, p)
        if not os.path.exists(p):
            warnings.append(f"参考文件不存在（运行前请确认）: {key} = {cfg[key]}")
    for w in warnings:
        print(f"[config 警告] {w}")


validate_config(config)

# ---------------------------------------------------------------------
# 常用查询函数（各 rules/*.smk 共用）
# ---------------------------------------------------------------------

def assay_needs_dedup(seqtype):
    """按 assay 决定是否 picard 去重；CUT&Tag 默认不去重（保留 PCR 重复）。"""
    return bool(config["dedup"][seqtype])


def sample_bam(sample):
    suffix = "rmdup.bam" if assay_needs_dedup(SEQTYPE_OF[sample]) else "sorted.bam"
    return f"3.align/bowtie2/{sample}_{suffix}"


def group_bams(group, role):
    return [sample_bam(s) for s in GROUPS[group][role]]


def group_control_arg(wc):
    bams = group_bams(wc.group, "control")
    return "-c " + ",".join(bams) if bams else ""


def group_peak_file(group):
    if GROUPS[group]["peak_type"] == "broad":
        return f"4.peak/{group}_peaks.broadPeak"
    return f"4.peak/{group}_peaks.narrowPeak"


def _groups_of(assay, peak_type=None):
    return [g for g, v in GROUPS.items()
            if v["seqtype"] == assay and (peak_type is None or v["peak_type"] == peak_type)]


def _group_regex(groups):
    """将分组名列表编为通配符约束；空列表给出永不匹配的正则。"""
    if not groups:
        return "(?!x)x"
    return "(?:" + "|".join(re.escape(g) for g in groups) + ")"


# ---------------------------------------------------------------------
# 目标汇总
# ---------------------------------------------------------------------

BAM_TARGETS = [f"3.align/bowtie2/{s}_sorted.bam" for s in SAMPLES]
BAM_TARGETS += [f"3.align/bowtie2/{s}_rmdup.bam"
                for s in SAMPLES if assay_needs_dedup(SEQTYPE_OF[s])]

PEAK_TARGETS = [group_peak_file(g) for g in GROUPS]
BW_TARGETS = [f"4.peak/{g}_FE.bw" for g in GROUPS]

QC_TARGETS = ["2.cleandata/fastqc/multiqc/multiqc_report.html"]
if config["qc"]["nsc_rsc"]:
    QC_TARGETS += [f"5.QC/spp/{s}_NSC.txt" for s in SAMPLES]
if config["qc"]["frip"]:
    QC_TARGETS += ["5.QC/frip/FRiP_summary.tsv"]
if config["qc"]["deeptools"]:
    QC_TARGETS += [
        "5.QC_deeptools/multiBamSummary.npz",
        "5.QC_deeptools/heatmap_SpearmanCorr_readCounts.png",
        "5.QC_deeptools/PCA_readCounts.png",
        "5.QC_deeptools/fingerprints.png",
        "5.QC_deeptools/fragmentsize.png",
        "5.QC_deeptools/profile_scaled.png",
    ]


rule all:
    input:
        QC_TARGETS,
        BAM_TARGETS,
        PEAK_TARGETS,
        BW_TARGETS,
        "4.peak/anno_result/Peakanno_PeakDistributions.pdf",


# ---------------------------------------------------------------------
# 各模块规则
# ---------------------------------------------------------------------
include: os.path.join(REPO_DIR, "rules", "upstream.smk")
include: os.path.join(REPO_DIR, "rules", "dedup.smk")
include: os.path.join(REPO_DIR, "rules", "callpeak.smk")
include: os.path.join(REPO_DIR, "rules", "annotation.smk")
if config["qc"]["nsc_rsc"]:
    include: os.path.join(REPO_DIR, "rules", "spp_qc.smk")
if config["qc"]["frip"]:
    include: os.path.join(REPO_DIR, "rules", "frip.smk")
if config["qc"]["deeptools"]:
    include: os.path.join(REPO_DIR, "rules", "qc_deeptools.smk")
