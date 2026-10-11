#!/usr/bin/env python3
"""Print one spike-in scale factor per idxstats file, one "sample<TAB>factor"
pair per line: factor = 1e6 / max(spike_mapped, 1).

Spike contigs are identified with the same pattern convention as
spikein_summary.py (<= 3 chars: exact case-insensitive match; longer:
substring). Consumed by the per-sample bigWig rescaling process.

Usage:
    spikein_factor.py --patterns spike sample1.idxstats [sample2.idxstats ...]
"""
import argparse
import os


def is_spike_contig(name, patterns):
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
    ap.add_argument("--patterns", required=True)
    ap.add_argument("idxstats", nargs="+")
    args = ap.parse_args()

    patterns = [p for p in (x.strip() for x in args.patterns.split(",")) if p]
    for path in args.idxstats:
        sample = os.path.basename(path)
        for suffix in ("_idxstats.tsv", ".idxstats"):
            if sample.endswith(suffix):
                sample = sample[:-len(suffix)]
                break
        spike = 0
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 4 or cols[0] == "*":
                    continue
                if is_spike_contig(cols[0], patterns):
                    spike += int(cols[2])
        factor = 1e6 / max(spike, 1)
        print(f"{sample}\t{factor:.6f}")


if __name__ == "__main__":
    raise SystemExit(main())
