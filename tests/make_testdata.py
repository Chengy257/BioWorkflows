#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""生成微型 ChIP/CUT&Tag/ATAC 回归测试数据集（dry-run 回归与 --real-run 验证共用）。

纯标准库实现、固定随机种子（默认 42），同参数两次生成逐字节一致
（gzip 头部 mtime 固定为 0）。数据不入库，在测试机上即时生成。
产物（写入 --outdir 指定目录，路径均为相对数据工作目录的布局）：

    ref/genome.fa             2 条 100 kb 染色体（chr1/chr2）
    ref/genes.gtf             每染色体 30 个模拟基因（gene + exon 两级行，
                              exon 带 gene_id/transcript_id，可过 ChIPseeker）
    ref/genes.bed             BED6 基因模型（chrom/start/end/name/score/strand）
    ref/chrom.sizes           染色体长度表（bigWig 用）
    1.rawdata/{sample}_1.fq.gz 与 _2.fq.gz
                              PE 模拟 reads（长度 50，片段 150-300bp，
                              序列真实取自参考基因组——--real-run 可过 bowtie2）
    samples.csv               5 样本 6 列样本表（sample_id,role,group,seqtype,layout,peak_type）
    config.yaml               测试 config（相对路径、threads: 2、qc 三开关全 true）

样本：chip_treat_rep1/rep2 + chip_control（narrow 组 g1）、
      atac_treat_rep1/rep2（atac 组 g2）。
富集设计：chr1 上预置 3 个 2kb 峰区，treat 样本 70% 片段取自峰区，
control 全基因组均匀采样。

用法:
    python tests/make_testdata.py --outdir <目录> [--reads 50000] [--seed 42]
"""
import argparse
import gzip
import io
import os
import random

# ---------------------------------------------------------------------
# 固定参数（与实施计划 Task 4.2 一致，改前需同步 run_tests.py 断言）
# ---------------------------------------------------------------------
CHROM_LEN = 100000        # 单条染色体长度（bp）
N_CHROM = 2               # 染色体条数（chr1/chr2）
GENES_PER_CHROM = 30      # 每染色体基因数
READ_LEN = 50             # reads 长度
FRAG_MIN, FRAG_MAX = 150, 300   # 片段长度范围（bp）
GENE_START = 1000         # 首个基因起点（0-based）
GENE_SPACING = 3000       # 基因间距
GENE_LEN = 2000           # 基因跨度
PEAK_REGIONS = [(20000, 22000), (50000, 52000), (80000, 82000)]  # chr1 富集峰区（0-based 半开）
PEAK_PROB = 0.7           # treat 样本取自峰区的概率（control 均匀采样）
ERROR_RATE = 0.005        # 替换型测序错误率（少量错误避免病态完全重复）
SEED = 42                 # 默认随机种子

# 样本集：(sample_id, role, group, seqtype, peak_type)
# 命名与组名均满足 workflow 的 _NAME_RE（字母数字._-，无连续下划线 __）
SAMPLES = [
    ("chip_treat_rep1", "treat",   "g1", "chip", "narrow"),
    ("chip_treat_rep2", "treat",   "g1", "chip", "narrow"),
    ("chip_control",    "control", "g1", "chip", "narrow"),
    ("atac_treat_rep1", "treat",   "g2", "atac", "none"),
    ("atac_treat_rep2", "treat",   "g2", "atac", "none"),
]

BASES = "ACGT"
_COMP = str.maketrans("ACGT", "TGCA")

# 测试 config：键集覆盖 validate_config 全部必需键；路径相对数据工作目录。
# qc.nsc_rsc=true 在 CI dry-run 中只构建 DAG 不执行 spp，安全。
CONFIG_YAML = """\
# =====================================================================
# 合成测试数据集专用 config（由 tests/make_testdata.py 生成，勿手工编辑）
# 所有路径相对数据工作目录；用法: bash run.sh -P <工作目录> -c <本文件> -n
# =====================================================================

# ---------- 参考基因组（合成 2 x 100kb） ----------
genome_fa: "ref/genome.fa"     # bowtie2 建索引输入
gtf: "ref/genes.gtf"           # 峰注释用（ChIPseeker makeTxDbFromGFF 可解析）
bed: "ref/genes.bed"           # 预留（deeptools 用）
chromsize: "ref/chrom.sizes"   # bigWig 生成用
genome_size: "180000"          # MACS2 -g 有效基因组大小：两条 100kb 染色体

# ---------- 样本与资源 ----------
grouplist: "samples.csv"       # 生成器同目录输出的样本表
threads: 2

# ---------- 质量修剪（trim_galore） ----------
trim:
  quality: 25
  stringency: 3
  error_rate: 0.1
  extra: ""

# ---------- 比对 ----------
bowtie2_extra: "--end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700"
min_mapq: 30

# ---------- 去重 ----------
dedup:
  chip: true
  cuttag: false
  atac: true
  faire: true

# ---------- 峰调用 ----------
peak:
  keepdup: all
  qvalue: 0.05
  broad_cutoff: 0.05
  atac:
    mode: bampe
    shift: -100
    extsize: 200

# ---------- 注释与信号窗口 ----------
region_flank: 3000

# ---------- QC（三开关全开，覆盖全部 QC 分支的 dry-run） ----------
qc:
  nsc_rsc: true    # dry-run 不执行 spp 故安全；--real-run 需 phantompeakqualtools
  frip: true
  deeptools: true
"""


def revcomp(s):
    """反向互补（R2 取片段另一端的反向互补链）。"""
    return s.translate(_COMP)[::-1]


def mutate(seq, rng):
    """按 ERROR_RATE 引入替换型测序错误（保持序列可比对回参考基因组）。"""
    if ERROR_RATE <= 0:
        return seq
    out = []
    for b in seq:
        if rng.random() < ERROR_RATE:
            out.append(rng.choice(BASES))
        else:
            out.append(b)
    return "".join(out)


def gzip_text(path):
    """固定 mtime 的 gzip 文本写句柄（保证同参数两次生成逐字节一致，LF 行尾）。"""
    raw = gzip.GzipFile(path, "wb", compresslevel=6, mtime=0)
    return io.TextIOWrapper(raw, encoding="ascii", newline="\n")


def build_reference(seed):
    """返回 (染色体序列 dict, 基因结构列表)。

    染色体序列由固定种子生成；基因位置为确定性布局（等距、正负链交替），
    不消耗随机数——因此 genome.fa 与 --reads 参数无关，可独立复现。
    """
    rng = random.Random(seed)
    chroms = {f"chr{i + 1}": "".join(rng.choices(BASES, k=CHROM_LEN))
              for i in range(N_CHROM)}
    genes = []
    for i in range(N_CHROM * GENES_PER_CHROM):
        chrom = f"chr{i // GENES_PER_CHROM + 1}"
        idx = i % GENES_PER_CHROM
        strand = "+" if i % 2 == 0 else "-"
        gid = f"gene{i + 1}"
        s0 = GENE_START + idx * GENE_SPACING          # 基因起点（0-based）
        exons = [(s0, s0 + 800), (s0 + 1200, s0 + GENE_LEN)]  # 两外显子（0-based 半开）
        genes.append({"id": gid, "chrom": chrom, "strand": strand,
                      "start": s0, "end": s0 + GENE_LEN, "exons": exons})
    return chroms, genes


def write_reference(outdir, chroms, genes):
    """写参考四件套 ref/{genome.fa, genes.gtf, genes.bed, chrom.sizes}（LF 行尾）。"""
    ref_dir = os.path.join(outdir, "ref")
    os.makedirs(ref_dir, exist_ok=True)

    with open(os.path.join(ref_dir, "genome.fa"), "w",
              encoding="ascii", newline="\n") as fh:
        for c in sorted(chroms):
            fh.write(f">{c}\n")
            s = chroms[c]
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")

    # GTF：最简 gene + exon 两级行；exon 带 transcript_id 供 makeTxDbFromGFF 分组
    with open(os.path.join(ref_dir, "genes.gtf"), "w",
              encoding="ascii", newline="\n") as fh:
        for g in genes:
            attr = f'gene_id "{g["id"]}";'
            fh.write(f'{g["chrom"]}\ttest\tgene\t{g["start"] + 1}\t{g["end"]}\t'
                     f'.\t{g["strand"]}\t.\t{attr}\n')
            exon_attr = attr + f' transcript_id "{g["id"]}.1";'
            for s, e in g["exons"]:
                fh.write(f'{g["chrom"]}\ttest\texon\t{s + 1}\t{e}\t'
                         f'.\t{g["strand"]}\t.\t{exon_attr}\n')

    # BED6：chrom/start/end/name/score/strand（ChIPseeker/rtracklayer 兼容）
    with open(os.path.join(ref_dir, "genes.bed"), "w",
              encoding="ascii", newline="\n") as fh:
        for g in genes:
            fh.write(f'{g["chrom"]}\t{g["start"]}\t{g["end"]}\t'
                     f'{g["id"]}\t0\t{g["strand"]}\n')

    with open(os.path.join(ref_dir, "chrom.sizes"), "w",
              encoding="ascii", newline="\n") as fh:
        for c in sorted(chroms):
            fh.write(f"{c}\t{len(chroms[c])}\n")


def sample_fragment(chroms, rng, role):
    """按角色采样一个基因组片段：treat 以 PEAK_PROB 概率取自 chr1 峰区，
    其余（含 control）全基因组均匀采样。序列必为参考基因组真实子串。"""
    frag_len = rng.randint(FRAG_MIN, FRAG_MAX)
    if role == "treat" and rng.random() < PEAK_PROB:
        chrom = "chr1"
        rs, re_ = rng.choice(PEAK_REGIONS)
        start = rng.randint(rs, re_ - frag_len)   # 片段完整落入峰区
    else:
        chrom = f"chr{rng.randint(1, N_CHROM)}"
        start = rng.randint(0, CHROM_LEN - frag_len)
    return chroms[chrom][start:start + frag_len]


def write_fastqs(outdir, chroms, genes, reads_per_sample, seed):
    """按样本生成 PE reads（1.rawdata/{sample}_1.fq.gz 与 _2.fq.gz）。

    每样本独立 rng（seed + 样本序号），与样本生成顺序解耦；
    R2 取片段另一端的反向互补，质量行固定为 'I' 重复。
    """
    raw_dir = os.path.join(outdir, "1.rawdata")
    os.makedirs(raw_dir, exist_ok=True)

    for k, (sid, role, _grp, _seqtype, _pt) in enumerate(SAMPLES):
        rng = random.Random(seed + k)
        fq1 = gzip_text(os.path.join(raw_dir, f"{sid}_1.fq.gz"))
        fq2 = gzip_text(os.path.join(raw_dir, f"{sid}_2.fq.gz"))
        try:
            for n in range(reads_per_sample):
                frag = sample_fragment(chroms, rng, role)
                r1 = mutate(frag[:READ_LEN], rng)
                r2 = revcomp(mutate(frag[-READ_LEN:], rng))
                q = "I" * READ_LEN
                fq1.write(f"@r{n:07d} 1:N:0:1\n{r1}\n+\n{q}\n")
                fq2.write(f"@r{n:07d} 2:N:0:1\n{r2}\n+\n{q}\n")
        finally:
            fq1.close()
            fq2.close()


def write_samples(outdir):
    """写 6 列样本表 samples.csv（表头与 workflow REQUIRED_COLUMNS 一致）。"""
    path = os.path.join(outdir, "samples.csv")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("sample_id,role,group,seqtype,layout,peak_type\n")
        for sid, role, grp, seqtype, pt in SAMPLES:
            fh.write(f"{sid},{role},{grp},{seqtype},PE,{pt}\n")


def write_config(outdir):
    """写测试 config.yaml（相对路径，供 run.sh -c 使用）。"""
    path = os.path.join(outdir, "config.yaml")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(CONFIG_YAML)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True,
                    help="输出根目录（必填；目录已存在则复用）")
    ap.add_argument("--reads", type=int, default=50000,
                    help="每样本 PE 读对数（默认 50000；CI 回归用 2000）")
    ap.add_argument("--seed", type=int, default=SEED,
                    help=f"随机种子（默认 {SEED}，保证可复现）")
    args = ap.parse_args()
    if args.reads < 1:
        ap.error("--reads 必须是正整数")

    os.makedirs(args.outdir, exist_ok=True)
    chroms, genes = build_reference(args.seed)
    write_reference(args.outdir, chroms, genes)
    write_fastqs(args.outdir, chroms, genes, args.reads, args.seed)
    write_samples(args.outdir)
    write_config(args.outdir)

    print(f"[make_testdata] 染色体 {N_CHROM} × {CHROM_LEN}bp，基因 {len(genes)} 个，"
          f"样本 {len(SAMPLES)} × {args.reads} PE 读对（seed={args.seed}）")
    print(f"[make_testdata] 产物根目录: {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()
