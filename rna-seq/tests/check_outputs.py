#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""测试输出断言器（对应路线图 3.1：关键输出断言 + mapping_stat 对齐校验 + DEG 方向校验）。

用法:
    python3 check_outputs.py --work <测试项目目录> --pipeline <pipeline> [--reads N]

按 pipeline 逐项断言 results/ 下的关键产物；任何断言失败打印 [FAIL] 并以退出码 1 结束。
"""
import argparse
import csv
import os
import sys

CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok), detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" — {detail}" if (detail and not ok) else ""))


def read_tsv(path):
    with open(path, encoding="utf-8") as fh:
        return [line.rstrip("\n").split("\t") for line in fh if line.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    ap.add_argument("--pipeline", default="deg", choices=["upstream", "deg", "as"])
    ap.add_argument("--reads", type=int, default=50000, help="每样本预期 PE 读对数（trim 报告总数校验）")
    args = ap.parse_args()

    res = os.path.join(args.work, "results")
    samples = ["ctrl_1", "ctrl_2", "drugA_1", "drugA_2"]

    check("results 根目录存在", os.path.isdir(res), res)

    # ---- 通用产物 ----
    check("software_versions.yaml", os.path.isfile(os.path.join(res, "software_versions.yaml")))
    check("multiqc 报告", os.path.isfile(os.path.join(res, "multiqc", "multiqc_report.html")))
    check("mapping_stat.xls", os.path.isfile(os.path.join(res, "3.align", "mapping_stat.xls")))
    for s in samples:
        bam = os.path.join(res, "3.align", f"{s}_Aligned.sortedByCoord.out.bam")
        bai = bam + ".bai"
        check(f"BAM+BAI {s}", os.path.isfile(bam) and os.path.isfile(bai))

    # ---- mapping_stat 对齐与数值校验 ----
    ms_path = os.path.join(res, "3.align", "mapping_stat.xls")
    if os.path.isfile(ms_path):
        rows = read_tsv(ms_path)
        header = ["ID", "Total_Reads", "Clean_Reads", "Uniquely_mapped", "Uniquely_mapped_ratio"]
        check("mapping_stat 表头", rows and rows[0] == header, str(rows[:1]))
        data = rows[1:]
        check("mapping_stat 行数 = 样本数", len(data) == len(samples), f"{len(data)} 行")
        check("mapping_stat 每行 5 列（对齐校验）", all(len(r) == 5 for r in data),
              str([r for r in data if len(r) != 5]))
        ids = {r[0] for r in data}
        check("mapping_stat 样本 id 完整（无截断）", ids == set(samples), str(sorted(ids)))
        totals_ok = all(r[1].isdigit() for r in data)
        check("mapping_stat Total_Reads 为数值（无 NA）", totals_ok, str([r[1] for r in data]))
        if totals_ok:
            check("mapping_stat Total_Reads = 输入读对数",
                  all(int(r[1]) == args.reads for r in data),
                  str({r[0]: r[1] for r in data}))
        ratios = [r[4].rstrip("%") for r in data]
        check("mapping_stat 比对率可解析", all(v.replace(".", "", 1).isdigit() for v in ratios),
              str(ratios))

    # ---- 上游定量 ----
    if args.pipeline in ("upstream", "deg"):
        cm = os.path.join(res, "4.expression", "count.matrix.tsv")
        check("count.matrix.tsv", os.path.isfile(cm))
        if os.path.isfile(cm):
            rows = read_tsv(cm)
            check("count.matrix 样本列 = 样本表", rows and rows[0] == ["id"] + samples,
                  str(rows[0] if rows else None))
            check("count.matrix 基因数 >= 40", len(rows) - 1 >= 40, f"{max(len(rows) - 1, 0)} 行")
            check("count.matrix 值均为整数", all(c.lstrip('-').isdigit()
                                                 for r in rows[1:] for c in r[1:]))
        check("GeneExpression_TPM.xls", os.path.isfile(os.path.join(res, "4.expression", "GeneExpression_TPM.xls")))
        check("GeneExpression_FPKM.xls", os.path.isfile(os.path.join(res, "4.expression", "GeneExpression_FPKM.xls")))
        check("GeneCount_Assigned_logs.xls", os.path.isfile(os.path.join(res, "4.expression", "GeneCount_Assigned_logs.xls")))

    # ---- as 管线 ----
    if args.pipeline == "as":
        check("merged.gtf", os.path.isfile(os.path.join(res, "4.assembly", "stringtie", "merged.gtf")))
        for s in samples:
            check(f"isoform 定量 {s}", os.path.isfile(os.path.join(res, "4.assembly", "isoform", f"{s}.tab")))

    # ---- DEG 管线 ----
    if args.pipeline == "deg":
        for f in ["5.DEG/flag.log", "5.DEG/GO_KEGG_enrich/flag.log", "6.DEGcompare/flag.log"]:
            check(f"flag {f}", os.path.isfile(os.path.join(res, *f.split("/"))))
        der_dir = os.path.join(res, "5.DEG", "Diff_Expr_Analysis_Results")
        deg_tsv = os.path.join(der_dir, "drugA_vs_control_DESeq2.output.tsv")
        check("DEG 结果表（拼写修正后目录）", os.path.isfile(deg_tsv), deg_tsv)
        for f in ["drugA_vs_control_VolcanoPlot.pdf", "drugA_vs_control_MAPlot.pdf",
                  "DESeq2.normalized.vst.PCA_plot.pdf", "DESeq2.normalized.vst.Pearson_heatmap.pdf",
                  "sessionInfo.txt"]:
            check(f"图表/记录 {f}", os.path.isfile(os.path.join(der_dir, f)))
        if os.path.isfile(deg_tsv):
            rows = read_tsv(deg_tsv)
            header = rows[0]
            check("DEG 表头含 gene_id", header and header[0] == "gene_id", str(header[:3]))
            try:
                i_lfc = header.index("log2FoldChange")
                i_type = header.index("type")
            except ValueError:
                i_lfc = i_type = None
            if i_lfc is not None:
                by_gene = {r[0]: r for r in rows[1:]}
                truth_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                          "data", "truth_degenes.tsv")
                if os.path.isfile(truth_path):
                    truth = {r[0]: r[1] for r in read_tsv(truth_path)[1:]}
                    ok_up = all(float(by_gene[g][i_lfc]) > 1.5
                                for g, d in truth.items() if d == "up" and g in by_gene)
                    ok_down = all(float(by_gene[g][i_lfc]) < -1.5
                                  for g, d in truth.items() if d == "down" and g in by_gene)
                    check("DEG 方向校验：预期上调基因 LFC > 1.5", ok_up,
                          str({g: by_gene[g][i_lfc] for g, d in truth.items() if d == "up" and g in by_gene}))
                    check("DEG 方向校验：预期下调基因 LFC < -1.5", ok_down,
                          str({g: by_gene[g][i_lfc] for g, d in truth.items() if d == "down" and g in by_gene}))
                else:
                    check("truth_degenes.tsv 存在", False, truth_path)
                types = {r[i_type] for r in rows[1:]}
                check("DEG type 取值合法", types <= {"Up", "Down", "Unsig"}, str(types))

    n_fail = sum(1 for _, ok, _ in CHECKS if not ok)
    print(f"\n[check_outputs] {len(CHECKS)} 项断言，{n_fail} 项失败")
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
