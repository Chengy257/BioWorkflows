#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Test output asserter (key output assertions + mapping_stat consistency checks + DEG direction checks).

Usage:
    python3 check_outputs.py --work <test project directory> --pipeline <pipeline> [--reads N]

Asserts the key outputs under results/ one by one per pipeline; any failed assertion prints [FAIL]
and the script exits with status 1.
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
    ap.add_argument("--reads", type=int, default=50000, help="expected PE read pairs per sample (verified against trim report totals)")
    args = ap.parse_args()

    res = os.path.join(args.work, "results")
    samples = ["ctrl_1", "ctrl_2", "drugA_1", "drugA_2"]

    check("results root directory exists", os.path.isdir(res), res)

    # ---- common outputs ----
    check("software_versions.yaml", os.path.isfile(os.path.join(res, "software_versions.yaml")))
    check("multiqc report", os.path.isfile(os.path.join(res, "multiqc", "multiqc_report.html")))
    check("mapping_stat.xls", os.path.isfile(os.path.join(res, "3.align", "mapping_stat.xls")))
    for s in samples:
        bam = os.path.join(res, "3.align", f"{s}_Aligned.sortedByCoord.out.bam")
        bai = bam + ".bai"
        check(f"BAM+BAI {s}", os.path.isfile(bam) and os.path.isfile(bai))

    # ---- mapping_stat consistency and numeric checks ----
    ms_path = os.path.join(res, "3.align", "mapping_stat.xls")
    if os.path.isfile(ms_path):
        rows = read_tsv(ms_path)
        header = ["ID", "Total_Reads", "Clean_Reads", "Uniquely_mapped", "Uniquely_mapped_ratio"]
        check("mapping_stat header", rows and rows[0] == header, str(rows[:1]))
        data = rows[1:]
        check("mapping_stat row count = sample count", len(data) == len(samples), f"{len(data)} rows")
        check("mapping_stat every row has 5 columns (alignment check)", all(len(r) == 5 for r in data),
              str([r for r in data if len(r) != 5]))
        ids = {r[0] for r in data}
        check("mapping_stat sample ids complete (no truncation)", ids == set(samples), str(sorted(ids)))
        totals_ok = all(r[1].isdigit() for r in data)
        check("mapping_stat Total_Reads numeric (no NA)", totals_ok, str([r[1] for r in data]))
        if totals_ok:
            check("mapping_stat Total_Reads = input read pairs",
                  all(int(r[1]) == args.reads for r in data),
                  str({r[0]: r[1] for r in data}))
        ratios = [r[4].rstrip("%") for r in data]
        check("mapping_stat alignment ratios parseable", all(v.replace(".", "", 1).isdigit() for v in ratios),
              str(ratios))

    # ---- upstream quantification ----
    if args.pipeline in ("upstream", "deg"):
        cm = os.path.join(res, "4.expression", "count.matrix.tsv")
        check("count.matrix.tsv", os.path.isfile(cm))
        if os.path.isfile(cm):
            rows = read_tsv(cm)
            check("count.matrix sample columns = sample table", rows and rows[0] == ["id"] + samples,
                  str(rows[0] if rows else None))
            check("count.matrix gene count >= 40", len(rows) - 1 >= 40, f"{max(len(rows) - 1, 0)} rows")
            check("count.matrix values are all integers", all(c.lstrip('-').isdigit()
                                                 for r in rows[1:] for c in r[1:]))
        check("GeneExpression_TPM.xls", os.path.isfile(os.path.join(res, "4.expression", "GeneExpression_TPM.xls")))
        check("GeneExpression_FPKM.xls", os.path.isfile(os.path.join(res, "4.expression", "GeneExpression_FPKM.xls")))
        check("GeneCount_Assigned_logs.xls", os.path.isfile(os.path.join(res, "4.expression", "GeneCount_Assigned_logs.xls")))

    # ---- as pipeline ----
    if args.pipeline == "as":
        check("merged.gtf", os.path.isfile(os.path.join(res, "4.assembly", "stringtie", "merged.gtf")))
        for s in samples:
            check(f"isoform quantification {s}", os.path.isfile(os.path.join(res, "4.assembly", "isoform", f"{s}.tab")))

    # ---- DEG pipeline ----
    if args.pipeline == "deg":
        for f in ["5.DEG/flag.log", "5.DEG/GO_KEGG_enrich/flag.log", "6.DEGcompare/flag.log"]:
            check(f"flag {f}", os.path.isfile(os.path.join(res, *f.split("/"))))
        der_dir = os.path.join(res, "5.DEG", "Diff_Expr_Analysis_Results")
        deg_tsv = os.path.join(der_dir, "drugA_vs_control_DESeq2.output.tsv")
        check("DEG result table (post-spelling-fix directory)", os.path.isfile(deg_tsv), deg_tsv)
        for f in ["drugA_vs_control_VolcanoPlot.pdf", "drugA_vs_control_MAPlot.pdf",
                  "DESeq2.normalized.vst.PCA_plot.pdf", "DESeq2.normalized.vst.Pearson_heatmap.pdf",
                  "sessionInfo.txt"]:
            check(f"plot/record {f}", os.path.isfile(os.path.join(der_dir, f)))
        if os.path.isfile(deg_tsv):
            rows = read_tsv(deg_tsv)
            header = rows[0]
            check("DEG header contains gene_id", header and header[0] == "gene_id", str(header[:3]))
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
                    check("DEG direction check: expected up-regulated genes have LFC > 1.5", ok_up,
                          str({g: by_gene[g][i_lfc] for g, d in truth.items() if d == "up" and g in by_gene}))
                    check("DEG direction check: expected down-regulated genes have LFC < -1.5", ok_down,
                          str({g: by_gene[g][i_lfc] for g, d in truth.items() if d == "down" and g in by_gene}))
                else:
                    check("truth_degenes.tsv exists", False, truth_path)
                types = {r[i_type] for r in rows[1:]}
                check("DEG type values are legal", types <= {"Up", "Down", "Unsig"}, str(types))

    n_fail = sum(1 for _, ok, _ in CHECKS if not ok)
    print(f"\n[check_outputs] {len(CHECKS)} assertions, {n_fail} failed")
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
