#!/usr/bin/env python3
"""Summarize organelle (chloroplast/mitochondrial) read fractions from
samtools idxstats tables (one per sample, produced by the organelle stage).

Contig classification (keep in sync with _is_organelle_contig in
workflow/rules/common.smk): patterns of <= 3 chars must equal the contig
name case-insensitively (so human ALT contigs are not swallowed); longer
patterns match as substrings. The '*' catch-all row idxstats appends for
unmapped reads is ignored; fractions use mapped reads only.

Usage:
    organelle_summary.py --patterns chrc,chrm,pt,mt --out summary.tsv \
        --mqc summary_mqc.tsv sample1_idxstats.tsv [sample2_... ...]
"""
import argparse
import os


def is_organelle_contig(name, patterns):
    lowered = name.lower()
    for pattern in patterns:
        p = str(pattern).lower()
        if len(p) <= 3:
            if lowered == p:
                return True
        elif p in lowered:
            return True
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--patterns", required=True,
                    help="comma-separated contig-name patterns")
    ap.add_argument("--out", required=True, help="summary TSV output path")
    ap.add_argument("--mqc", required=True, help="MultiQC custom-content TSV output path")
    ap.add_argument("idxstats", nargs="+", help="per-sample idxstats TSV files")
    args = ap.parse_args()

    patterns = [p for p in (x.strip() for x in args.patterns.split(",")) if p]
    header = ["sample", "mapped_reads", "organelle_reads", "organelle_fraction",
              "matched_contigs"]
    rows = []
    for path in args.idxstats:
        sample = os.path.basename(path)
        if sample.endswith("_idxstats.tsv"):
            sample = sample[:-len("_idxstats.tsv")]
        total = 0
        org = 0
        matched = []
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 4 or cols[0] == "*":
                    continue
                name, mapped = cols[0], int(cols[2])
                total += mapped
                if is_organelle_contig(name, patterns):
                    org += mapped
                    matched.append(f"{name}:{mapped}")
        fraction = (org / total) if total else 0.0
        rows.append([sample, total, org, f"{fraction:.4f}",
                     ";".join(matched) if matched else "-"])

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="\n", encoding="utf-8") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    with open(args.mqc, "w", newline="\n", encoding="utf-8") as fh:
        fh.write("# id: 'organelle_table'\n")
        fh.write("# section_name: 'Organelle read fraction (chloroplast/mitochondrion)'\n")
        fh.write("# format: 'tsv'\n")
        fh.write("# plot_type: 'table'\n")
        fh.write("# pconfig: {'id': 'organelle_table', 'title': 'Organelle reads'}\n")
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    print(f"[organelle_summary] {len(rows)} sample(s) summarized -> {args.out}")


if __name__ == "__main__":
    main()
