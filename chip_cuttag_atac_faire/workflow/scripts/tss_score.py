#!/usr/bin/env python3
"""Compute the TSS enrichment score from a deeptools computeMatrix
reference-point matrix (one matrix per sample, produced by the tss stage).

Definition (kept close to the ENCODE idea, on the aggregated profile):
  - per-bin mean signal across all TSS rows builds the sample profile;
  - the flanking baseline is the mean of the outermost 10 bins on each side;
  - TSSE = max(profile / baseline); the normalized value of the two bins
    straddling the TSS (center_norm) is also reported.
A zero/empty baseline yields TSSE = NA instead of a division by zero.

Usage:
    tss_score.py <matrix.gz> <sample> <out.tsv>
"""
import gzip
import sys

FLANK_BINS = 10  # outermost bins on each side forming the baseline


def main():
    if len(sys.argv) != 4:
        sys.exit("Usage: tss_score.py <matrix.gz> <sample> <out.tsv>")
    path, sample, out = sys.argv[1], sys.argv[2], sys.argv[3]

    opener = gzip.open if path.endswith(".gz") else open
    sums = None
    n = 0
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("@"):
                continue  # computeMatrix header line
            cols = line.rstrip("\n").split("\t")
            vals = cols[6:]  # 6 leading region columns, then one column per bin
            if not vals:
                continue
            if sums is None:
                sums = [0.0] * len(vals)
            for i, v in enumerate(vals):
                try:
                    sums[i] += float(v)
                except ValueError:
                    sums[i] += 0.0
            n += 1
    if not n or not sums:
        sys.exit(f"{path}: no matrix rows found")

    profile = [s / n for s in sums]
    nbins = len(profile)
    if nbins >= 2 * FLANK_BINS:
        baseline = (sum(profile[:FLANK_BINS]) + sum(profile[-FLANK_BINS:])) / (2 * FLANK_BINS)
    else:
        baseline = sum(profile) / nbins
    if baseline > 0:
        norm = [p / baseline for p in profile]
        tsse = f"{max(norm):.4f}"
        center = (norm[nbins // 2 - 1] + norm[nbins // 2]) / 2
        center_s = f"{center:.4f}"
        baseline_s = f"{baseline:.6f}"
    else:
        tsse = center_s = "NA"
        baseline_s = f"{baseline:.6f}"

    with open(out, "w", newline="\n", encoding="utf-8") as fh:
        fh.write(f"{sample}\t{tsse}\t{center_s}\t{baseline_s}\t{n}\n")
    print(f"[tss_score] {sample}: TSSE={tsse} over {n} TSS rows")


if __name__ == "__main__":
    main()
