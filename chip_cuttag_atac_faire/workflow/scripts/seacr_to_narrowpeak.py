#!/usr/bin/env python3
"""Convert a SEACR peak bed to the standard 10-column narrowPeak contract.

SEACR (Sparse Enrichment Analysis for CUT&RUN, Yo et al. 2021,
https://github.com/yeolab/SEACR) writes one 6-column bed per group:

    chrom  start  end  AUC(total signal)  max_signal  max_signal_region

Every downstream consumer of this workflow (FRiP, ChIPseeker annotation,
blacklist filtering, motif enrichment, DiffBind) expects the UCSC
narrowPeak/broadPeak 10-column layout, so the seacr_callpeak rule converts
the file right after the SEACR run:

    chrom  start  end  name  score  strand  signalValue  pValue  qValue  peak

with:
  name         <prefix>_<n> (1-based peak number; the prefix is the group name)
  score        int(min(1000, round(AUC))), clamped to 0-1000 (narrowPeak spec:
               0-1000; SEACR's AUC is an unbounded area so large peaks
               saturate at 1000)
  strand       "." (SEACR emits no strand)
  signalValue  the AUC column, passed through as a float
  pValue/qValue 0 (SEACR emits no per-peak p/q values)
  peak         0 (SEACR emits no summit; the region start stands in)

Usage:
    python seacr_to_narrowpeak.py <input.seacr.bed> <output.narrowPeak> <name_prefix>

Exits non-zero with a clear message on malformed input (a line with fewer
than 4 whitespace-separated columns or non-numeric coordinates/signal).
"""
import math
import sys


def convert(src_path, dst_path, prefix):
    rows = []
    with open(src_path, encoding="utf-8") as src:
        for lineno, line in enumerate(src, start=1):
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(
                    f"{src_path} line {lineno}: expected >= 4 columns "
                    f"(chrom start end AUC ...), got {len(fields)}: {line!r}")
            chrom = fields[0]
            try:
                start, end = int(fields[1]), int(fields[2])
                auc = float(fields[3])
            except ValueError as exc:
                raise ValueError(
                    f"{src_path} line {lineno}: non-numeric coordinates/signal "
                    f"({exc})") from exc
            if math.isnan(auc) or math.isinf(auc):
                raise ValueError(
                    f"{src_path} line {lineno}: AUC must be a finite number, "
                    f"got {fields[3]!r}")
            if start < 0 or end <= start:
                raise ValueError(
                    f"{src_path} line {lineno}: invalid interval {start}-{end}")
            # narrowPeak score is capped at 1000; SEACR's AUC is an unbounded
            # area under the coverage curve, so it saturates there
            score = max(0, min(1000, round(auc)))
            rows.append("\t".join((
                chrom, str(start), str(end),
                f"{prefix}_{len(rows) + 1}",
                str(score), ".",
                str(auc),
                "0", "0", "0",
            )))
    with open(dst_path, "w", newline="\n", encoding="utf-8") as dst:
        for row in rows:
            dst.write(row + "\n")
    return len(rows)


def main(argv):
    if len(argv) != 4:
        print(__doc__.strip(), file=sys.stderr)
        return 2
    src_path, dst_path, prefix = argv[1], argv[2], argv[3]
    if not prefix or any(c.isspace() for c in prefix):
        print(f"invalid name prefix: {prefix!r} (must be non-empty without "
              "whitespace)", file=sys.stderr)
        return 2
    try:
        count = convert(src_path, dst_path, prefix)
    except (OSError, ValueError) as exc:
        print(f"seacr_to_narrowpeak: {exc}", file=sys.stderr)
        return 1
    print(f"seacr_to_narrowpeak: wrote {count} peaks to {dst_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
