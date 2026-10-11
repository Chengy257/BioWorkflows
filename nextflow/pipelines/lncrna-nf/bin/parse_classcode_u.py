#!/usr/bin/env python3
"""Extract the novel (gffcompare classcode u) intergenic transcripts from a
gffcompare *.annotated.gtf.

Legacy semantics (rna-seq lncRNA_functions.sh runGFFcompare): keep the
transcripts whose gffcompare class_code attribute is 'u' (intergenic, no
overlap with any reference transcript).

Usage:
    parse_classcode_u.py --gtf annotated.gtf --out-gtf novel.gtf --out-ids ids.txt
"""
import argparse


def parse_attrs(field):
    attrs = {}
    for item in field.strip().rstrip(";").split(";"):
        item = item.strip()
        if not item or " " not in item:
            continue
        key, _, value = item.partition(" ")
        attrs[key] = value.strip('"')
    return attrs


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out-gtf", required=True)
    ap.add_argument("--out-ids", required=True)
    args = ap.parse_args()

    ids = []
    with open(args.gtf) as fh, open(args.out_gtf, "w") as out:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attrs = parse_attrs(parts[8])
            if attrs.get("class_code") == "u":
                tid = attrs.get("transcript_id", "")
                if tid and tid not in ids:
                    ids.append(tid)
                out.write(line)

    with open(args.out_ids, "w") as out:
        for tid in ids:
            out.write(tid + "\n")
    print(f"kept {len(ids)} classcode-u transcripts")


if __name__ == "__main__":
    main()
