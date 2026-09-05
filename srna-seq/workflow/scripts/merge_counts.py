#!/usr/bin/env python3
"""Merge per-sample per-class counts into count and RPM matrices.

RPM = reads per million mapped inside that class and sample
(count / class mapped total * 1e6; 0 when the class mapped no reads).
"""
import argparse
import os


def read_counts(path):
    table = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            feature, count = line.rstrip("\n").split("\t")
            table[feature] = int(count)
    return table


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--indir", required=True, help="results/4.expression root")
    ap.add_argument("--classes", nargs="+", required=True)
    ap.add_argument("--samples", nargs="+", required=True)
    args = ap.parse_args()
    combined = {}
    for klass in args.classes:
        per_sample = {s: read_counts(os.path.join(args.indir, klass, f"{s}_counts.txt"))
                      for s in args.samples}
        features = sorted({f for t in per_sample.values()
                           for f in t if f != "__mapped_total"})
        klass_dir = os.path.join(args.indir, klass)
        os.makedirs(klass_dir, exist_ok=True)
        with open(os.path.join(klass_dir, f"{klass}_counts.tsv"), "w") as cf, \
             open(os.path.join(klass_dir, f"{klass}_RPM.tsv"), "w") as rf:
            header = "feature\t" + "\t".join(args.samples) + "\n"
            cf.write(header)
            rf.write(header)
            for feature in features:
                counts = [per_sample[s].get(feature, 0) for s in args.samples]
                cf.write(feature + "\t" + "\t".join(map(str, counts)) + "\n")
                rpm = []
                for sample, count in zip(args.samples, counts):
                    mapped = per_sample[sample].get("__mapped_total", 0)
                    rpm.append(f"{count / mapped * 1e6:.3f}" if mapped else "0.000")
                rf.write(feature + "\t" + "\t".join(rpm) + "\n")
        combined[klass] = (features, per_sample)
    with open(os.path.join(args.indir, "all_classes_counts.tsv"), "w") as af:
        af.write("class\tfeature\t" + "\t".join(args.samples) + "\n")
        for klass in args.classes:
            features, per_sample = combined[klass]
            for feature in features:
                af.write(f"{klass}\t{feature}\t"
                         + "\t".join(str(per_sample[s].get(feature, 0))
                                     for s in args.samples) + "\n")


if __name__ == "__main__":
    main()
