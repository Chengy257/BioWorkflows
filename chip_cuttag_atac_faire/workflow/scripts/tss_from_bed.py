#!/usr/bin/env python3
"""Derive a strand-aware 1-bp TSS BED from a BED6 gene model.

'+' genes contribute their start coordinate, '-' genes their end-1
(0-based half-open 1-bp intervals). Lines starting with '#', 'track' or
'browser' are skipped. The input must have >= 6 columns with a strand
column (the workflow's `bed` config key is a BED6 gene model; for a GTF-only
deployment derive a BED6 first).

Usage:
    tss_from_bed.py <genes.bed> <out.bed>
"""
import sys


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: tss_from_bed.py <genes.bed> <out.bed>")
    src, dst = sys.argv[1], sys.argv[2]
    n = 0
    with open(src, encoding="utf-8") as fin, \
            open(dst, "w", newline="\n", encoding="utf-8") as fout:
        for lineno, line in enumerate(fin, start=1):
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                sys.exit(f"{src} line {lineno}: BED6 required "
                         "(6 tab-separated columns with a strand column)")
            chrom, start, end, strand = cols[0], int(cols[1]), int(cols[2]), cols[5]
            tss = end - 1 if strand == "-" else start
            fout.write(f"{chrom}\t{tss}\t{tss + 1}\n")
            n += 1
    if n == 0:
        sys.exit(f"{src}: no TSS derived (empty gene model?)")
    print(f"[tss_from_bed] {n} TSS written -> {dst}")


if __name__ == "__main__":
    main()
