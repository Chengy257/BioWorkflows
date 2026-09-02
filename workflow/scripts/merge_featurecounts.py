#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""合并各样本 featureCounts 定量结果。

用法:
    python merge_featurecounts.py <定量目录>

<定量目录>/ 下需有每样本的 <sample>.count 与 <sample>.log
（由 scripts/run_featurecounts.R 产生；.count 5 列: id/effLength/counts/fpkm/tpm）。

输出（写入同目录）:
    count.matrix.tsv            基因 × 样本 raw counts 矩阵（runDESeq2 输入）
    GeneExpression_TPM.xls      基因 × 样本 TPM 矩阵
    GeneExpression_FPKM.xls     基因 × 样本 FPKM 矩阵
    GeneCount_Assigned_logs.xls featureCounts 各状态（Assigned/Unassigned_*）统计

替代原合并 shell 脚本（其依赖的 njoin.sh / transposition.sh 不在仓库中），
行为对齐原输出格式。
"""
import glob
import os
import sys


def sample_name(path):
    """<sample>.count -> 样本名；'-' 替换为 '_'，与 runDESeq2 的 id 处理保持一致"""
    base = os.path.basename(path)
    for ext in (".count", ".log"):
        if base.endswith(ext):
            return base[:-len(ext)].replace("-", "_")
    return base.replace("-", "_")


def read_count_file(path):
    """返回 [(gene_id, counts, fpkm, tpm), ...]"""
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header[:5] != ["id", "effLength", "counts", "fpkm", "tpm"]:
            sys.exit(f"[ERROR] {path} 表头异常: {header[:5]}")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            rows.append((parts[0], parts[2], parts[3], parts[4]))
    return rows


def read_log_file(path):
    """featureCounts $stat 表 -> {status: 数量}；跳过表头行（Status/...）"""
    stats = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or parts[0] in ("Status", ""):
                continue
            stats.setdefault(parts[0], parts[1])
    return stats


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    out_dir = sys.argv[1].rstrip("/\\") or "."

    count_files = sorted(glob.glob(os.path.join(out_dir, "*.count")))
    if not count_files:
        sys.exit(f"[ERROR] {out_dir} 下未找到 *.count 文件")

    ## ---- 合并 count/fpkm/tpm 矩阵（以首个样本的基因顺序为基准）----
    gene_order = []
    values = {}  # gene_id -> {sample: (counts, fpkm, tpm)}
    for cf in count_files:
        s = sample_name(cf)
        for gid, counts, fpkm, tpm in read_count_file(cf):
            if gid not in values:
                values[gid] = {}
                gene_order.append(gid)
            values[gid][s] = (counts, fpkm, tpm)
    samples = [sample_name(cf) for cf in count_files]

    def write_matrix(fname, idx):
        with open(os.path.join(out_dir, fname), "w") as fh:
            fh.write("id\t" + "\t".join(samples) + "\n")
            for gid in gene_order:
                row = [values[gid].get(s, ("NA", "NA", "NA"))[idx] for s in samples]
                fh.write(gid + "\t" + "\t".join(row) + "\n")

    write_matrix("count.matrix.tsv", 0)
    write_matrix("GeneExpression_FPKM.xls", 1)
    write_matrix("GeneExpression_TPM.xls", 2)

    ## ---- 合并 featureCounts 状态统计（行=样本，列=状态）----
    log_files = sorted(glob.glob(os.path.join(out_dir, "*.log")))
    status_order = []
    stat_values = {}  # sample -> {status: value}
    for lf in log_files:
        s = sample_name(lf)
        stat_values[s] = read_log_file(lf)
        for k in stat_values[s]:
            if k not in status_order:
                status_order.append(k)
    with open(os.path.join(out_dir, "GeneCount_Assigned_logs.xls"), "w") as fh:
        fh.write("id\t" + "\t".join(status_order) + "\n")
        for s in samples:
            row = [stat_values.get(s, {}).get(k, "NA") for k in status_order]
            fh.write(s + "\t" + "\t".join(row) + "\n")

    print(f"[merge_featurecounts] {len(gene_order)} genes × {len(samples)} samples -> {out_dir}")


if __name__ == "__main__":
    main()
