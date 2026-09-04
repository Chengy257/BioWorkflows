###############################################
## 通用定义：路径常量、样本解析、fastq 探测、辅助函数
## 本文件只定义 Python 侧内容与通用规则；被 workflow/Snakefile 首个 include
###############################################
import os

SCRIPTS = os.path.join(WORKFLOW_DIR, "scripts")
ENVS = os.path.join(WORKFLOW_DIR, "envs")

## 运行产物根目录（相对运行目录，config: results_dir）
RD = config.get("results_dir", "results").rstrip("/") + "/"


def R(path=""):
    """拼接 results_dir 下的产物路径（规则 input/output 在解析期求值）"""
    return f"{RD}{path}"


def res(key, default=None):
    """物种资源读取：顶层显式键优先（Snakefile 已将 species.yaml 预设展开到顶层）"""
    return config.get(key, default)


def lres(key):
    """lncRNA 附属配置：config['lncrna'][key] 优先，回落到顶层键；仍缺失则报错"""
    v = (config.get("lncrna") or {}).get(key)
    if v in (None, "") and key in config:
        v = config[key]
    if v in (None, ""):
        raise ValueError(f"config 缺少 lncrna.{key}（见 config.yaml 的 lncrna 段）")
    return v


def get_samples():
    """样本表（config: SampleListFile，相对运行目录；表头须含 id,group）"""
    ids = []
    with open(config["SampleListFile"], "r") as samples_list:
        next(samples_list)
        for line in samples_list:
            line = line.strip().split(",")
            ids.append(line[0])
    return ids


SAMPLES = get_samples()


def _raw_reads(sample):
    """探测样本原始 fastq：先 PE（<id>_1/<id>_2），后 SE（<id>.fastq.gz / <id>.fq.gz）。
    使用精确文件名匹配，避免样本 id 互为前缀时 glob 误配。"""
    for pat1, pat2 in (("1.rawdata/{0}_1.fastq.gz", "1.rawdata/{0}_2.fastq.gz"),
                       ("1.rawdata/{0}_1.fq.gz", "1.rawdata/{0}_2.fq.gz")):
        if os.path.exists(pat1.format(sample)) and os.path.exists(pat2.format(sample)):
            return "PE", [pat1.format(sample), pat2.format(sample)]
    for pat in ("1.rawdata/{0}.fastq.gz", "1.rawdata/{0}.fq.gz"):
        if os.path.exists(pat.format(sample)):
            return "SE", [pat.format(sample)]
    raise ValueError(
        "sample {0}: 1.rawdata/ 下未找到原始 fastq"
        "（支持 {0}_1.fastq.gz+{0}_2.fastq.gz、{0}_1.fq.gz+{0}_2.fq.gz PE，"
        "或 {0}.fastq.gz / {0}.fq.gz SE）".format(sample))


def get_fastq(wildcards):
    """runSTAR 的输入依赖（trim 后文件），随文库类型自动切换。"""
    layout, _ = _raw_reads(wildcards.sample)
    s = wildcards.sample
    if layout == "PE":
        return [R(f"2.cleandata/trim/{s}_1_val_1.fq.gz"),
                R(f"2.cleandata/trim/{s}_2_val_2.fq.gz")]
    return [R(f"2.cleandata/trim/{s}_trimmed.fq.gz")]


def trimmed_reads(sample):
    """STAR --readFilesIn 参数：PE 两个文件，SE 一个文件。"""
    layout, _ = _raw_reads(sample)
    if layout == "PE":
        return f"{RD}2.cleandata/trim/{sample}_1_val_1.fq.gz {RD}2.cleandata/trim/{sample}_2_val_2.fq.gz"
    return f"{RD}2.cleandata/trim/{sample}_trimmed.fq.gz"


def is_paired_end(wildcards):
    return "True" if _raw_reads(wildcards.sample)[0] == "PE" else "False"


def trim_reports(sample):
    """某样本的 trim 报告（trim 规则统一重命名为规范名）。"""
    layout, _ = _raw_reads(sample)
    if layout == "PE":
        return [R(f"2.cleandata/trim/{sample}_1_trimming_report.txt"),
                R(f"2.cleandata/trim/{sample}_2_trimming_report.txt")]
    return [R(f"2.cleandata/trim/{sample}_trimming_report.txt")]


def all_trim_reports(wildcards=None):
    """全部样本的 trim 报告（聚合规则输入）。"""
    out = []
    for s in SAMPLES:
        out += trim_reports(s)
    return out


## STAR 参数：组装管线使用更严格的剪接过滤（与原 AS 管线一致），其余用标准参数；
## 两者均可经 config 的 star_extra_args 追加
STAR_ARGS_STANDARD = ("--twopassMode Basic --genomeLoad NoSharedMemory "
                      "--quantMode GeneCounts --outSAMattrIHstart 0 ")
STAR_ARGS_ASSEMBLY = ("--twopassMode Basic --outFilterType BySJout --alignIntronMin 20 "
                      "--alignIntronMax 5000 --alignMatesGapMax 5000 "
                      "--outFilterMatchNminOverLread 0.66 --outFilterScoreMinOverLread 0.66 "
                      "--winAnchorMultimapNmax 70 --seedSearchStartLmax 45 --outSAMattrIHstart 0 "
                      "--outSAMstrandField intronMotif --genomeLoad NoSharedMemory "
                      "--quantMode TranscriptomeSAM GeneCounts ")
STAR_ARGS = (STAR_ARGS_ASSEMBLY if PIPELINE == "as" else STAR_ARGS_STANDARD) + \
            config.get("star_extra_args", "")


###############################################
## 软件版本记录（P2-9）：conda 环境清单 + snakemake 版本 + git 提交
###############################################
rule software_versions:
    output:
        R("software_versions.yaml"),
    log:
        R("logs/software_versions.log.txt"),
    shell:
        "python {SCRIPTS}/collect_versions.py --envs {ENVS} --workflow {WORKFLOW_DIR} --out {output} >> {log} 2>&1"
