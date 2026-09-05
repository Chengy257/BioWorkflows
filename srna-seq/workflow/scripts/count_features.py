#!/usr/bin/env python3
"""Count primary alignments per reference feature in one bowtie1 SAM file."""
import argparse


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("sam")
    ap.add_argument("-o", "--out", required=True)
    args = ap.parse_args()
    counts = {}
    total = 0
    with open(args.sam) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            fields = line.split("\t")
            flag = int(fields[1])
            if flag & 0x100:          # secondary alignment
                continue
            rname = fields[2]
            if rname == "*":
                continue
            counts[rname] = counts.get(rname, 0) + 1
            total += 1
    with open(args.out, "w") as out:
        out.write("feature\treads\n")
        for feature in sorted(counts):
            out.write(f"{feature}\t{counts[feature]}\n")
        out.write(f"__mapped_total\t{total}\n")


if __name__ == "__main__":
    main()
