#!/usr/bin/env python3
"""Build the reproducible peak set of a group from its pairwise IDR results.

Deliberate, documented simplification of the ENCODE rescue/self-consistency
scheme (inherited from the legacy workflow): every pair of replicates gets
one IDR run (the idr tool's *idrValues.txt, a narrowPeak-like table whose
score column carries the IDR value); the per-group final set is the union of
the pairwise intervals, kept where at least --min-replicates DISTINCT
replicates are supported.

Pair identity comes from the file name, which the pipeline writes as
"<group>__<repA>_vs_<repB>.idrValues.txt": an interval supported by a pair's
result is supported by both replicates of that pair. Overlapping intervals
are merged per chromosome before support counting.

Output (narrowPeak layout): chrom, start, end, name(<group>_<n>), score
(= distinct-replicate support, capped at 1000), strand '.', signalValue
(max IDR score overlapping the merged interval), pValue/qValue 0, peak 0.

Usage:
    idr_union.py --group H3K4me3 --min-replicates 2 \
        --out H3K4me3_idr_consensus.narrowPeak \
        'H3K4me3__rep1_vs_rep2.idrValues.txt' ...
"""
import argparse
import os
import re

PAIR_RE = re.compile(r"^(?P<group>.+)__(?P<a>.+)_vs_(?P<b>.+)\.idrValues\.txt$")


def parse_pair(path):
    name = os.path.basename(path)
    m = PAIR_RE.match(name)
    if not m:
        raise SystemExit(
            f"idr_union: cannot parse replicate pair from file name {name!r} "
            f"(expected '<group>__<repA>_vs_<repB>.idrValues.txt')")
    return m.group("a"), m.group("b")


def read_intervals(path):
    intervals = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                continue
            chrom, start, end = cols[0], int(cols[1]), int(cols[2])
            try:
                score = float(cols[4])
            except (ValueError, IndexError):
                score = 0.0
            intervals.append((chrom, start, end, score))
    return intervals


def merge(intervals):
    """Merge (chrom, start, end) intervals per chromosome; returns
    {chrom: [(start, end), ...]} sorted and non-overlapping."""
    buckets = {}
    for chrom, start, end in intervals:
        buckets.setdefault(chrom, []).append((start, end))
    merged = {}
    for chrom, ivs in buckets.items():
        ivs.sort()
        norm = []
        for start, end in ivs:
            if norm and start <= norm[-1][1]:
                norm[-1] = (norm[-1][0], max(norm[-1][1], end))
            else:
                norm.append((start, end))
        merged[chrom] = norm
    return merged


def overlap_len(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--group", required=True)
    ap.add_argument("--min-replicates", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("idr_files", nargs="+")
    args = ap.parse_args()

    # interval -> (set of supporting replicates, max IDR score)
    support = {}
    for path in args.idr_files:
        rep_a, rep_b = parse_pair(path)
        for chrom, start, end, score in read_intervals(path):
            key = (chrom, start, end)
            entry = support.setdefault(key, [set(), 0.0])
            entry[0].update((rep_a, rep_b))
            entry[1] = max(entry[1], score)

    merged = merge(list(support.keys()))

    rows = []
    for chrom, ivs in sorted(merged.items()):
        for start, end in ivs:
            reps = set()
            best = 0.0
            for (c, s, e), (who, score) in support.items():
                if c == chrom and overlap_len(start, end, s, e) > 0:
                    reps |= who
                    best = max(best, score)
            if len(reps) >= args.min_replicates:
                rows.append((chrom, start, end, best, len(reps)))

    rows.sort(key=lambda r: (r[0], r[1], r[2]))
    with open(args.out, "w", newline="\n", encoding="utf-8") as out:
        out.write("\t".join(
            ["chrom", "start", "end", "name", "score", "strand",
             "signalValue", "pValue", "qValue", "peak"]) + "\n")
        for i, (chrom, start, end, best, nrep) in enumerate(rows, start=1):
            out.write("\t".join([
                chrom, str(start), str(end),
                f"{args.group}_{i}",
                str(min(1000, nrep)),
                ".",
                f"{best:g}",
                "0", "0", "0",
            ]) + "\n")
    print(f"idr_union: {len(rows)} reproducible peaks for {args.group}")


if __name__ == "__main__":
    raise SystemExit(main())
