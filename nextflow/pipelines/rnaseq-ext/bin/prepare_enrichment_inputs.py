#!/usr/bin/env python3
"""Prepare clusterProfiler enrichment inputs from a DA differential table.

Reads one nf-core/differentialabundance `*.deseq2.results.tsv` table and
emits, for one contrast:
    <outdir>/<contrast>_UP.DEGs.txt      significant up-regulated gene ids
    <outdir>/<contrast>_DOWN.DEGs.txt    significant down-regulated gene ids
    <outdir>/<contrast>_ALL.DEGs.txt     union of UP and DOWN
    <outdir>/<contrast>_ranked.tsv       full ranked list (gene_id + rank column)

Significance = padj <= padj_max AND |log2FC| >= log2fc_min (the legacy DEG
thresholds; both configurable - no silent absolute configuration).

Usage:
    python prepare_enrichment_inputs.py <de_table.tsv> <outdir> <contrast> \
        --id-col gene_id --log2fc-col log2FoldChange --padj-col padj \
        --stat-col stat --padj-max 0.05 --log2fc-min 1.0
"""
import argparse
import csv
import os
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("de_table")
    p.add_argument("outdir")
    p.add_argument("contrast")
    p.add_argument("--id-col", default="gene_id")
    p.add_argument("--log2fc-col", default="log2FoldChange")
    p.add_argument("--padj-col", default="padj")
    p.add_argument("--stat-col", default="stat")
    p.add_argument("--padj-max", type=float, default=0.05)
    p.add_argument("--log2fc-min", type=float, default=1.0)
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    with open(args.de_table) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {args.id_col, args.log2fc_col, args.padj_col}
        missing = required - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"[ERROR] {args.de_table}: missing column(s) {sorted(missing)} "
                     f"(available: {reader.fieldnames}); adjust the enrichment.*_col params")
        rows = list(reader)

    up, down, ranked = [], [], []
    for row in rows:
        gid = (row[args.id_col] or "").strip()
        if not gid:
            continue
        try:
            lfc = float(row[args.log2fc_col])
        except (TypeError, ValueError):
            lfc = float("nan")
        try:
            padj = float(row[args.padj_col])
        except (TypeError, ValueError):
            padj = 1.0  # rows without an adjusted p-value are never significant
        sig = padj <= args.padj_max
        if sig and lfc >= args.log2fc_min:
            up.append(gid)
        if sig and lfc <= -args.log2fc_min:
            down.append(gid)
        try:
            stat = float(row[args.stat_col])
        except (TypeError, ValueError):
            stat = lfc  # fall back to the fold-change ranking
        ranked.append((gid, stat))

    ranked.sort(key=lambda t: t[1], reverse=True)
    stem = os.path.join(args.outdir, args.contrast)
    for direction, ids in (("UP", up), ("DOWN", down), ("ALL", up + down)):
        with open(f"{stem}_{direction}.DEGs.txt", "w") as fh:
            fh.write("\n".join(ids) + ("\n" if ids else ""))
    with open(f"{stem}_ranked.tsv", "w") as fh:
        fh.write("gene_id\tstat\n")
        for gid, stat in ranked:
            fh.write(f"{gid}\t{stat}\n")

    print(f"[prepare_enrichment_inputs] {args.contrast}: UP={len(up)} DOWN={len(down)} "
          f"ranked={len(ranked)} -> {args.outdir}")


if __name__ == "__main__":
    main()
