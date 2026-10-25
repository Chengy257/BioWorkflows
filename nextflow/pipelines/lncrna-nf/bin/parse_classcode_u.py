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

    rows = []
    with open(args.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            tid = parse_attrs(parts[8]).get("transcript_id", "")
            rows.append((line, tid, parts[8]))
    # gffcompare writes class_code ONLY on transcript rows; exon rows carry
    # the same transcript_id. Keep the classcode-u transcripts AND their
    # exon rows (dropping exons zeroes the downstream exonic-length filter
    # and degrades every transcript to a single-exon span).
    keep = {tid for _line, tid, field in rows
            if tid and parse_attrs(field).get("class_code") == "u"}
    with open(args.out_gtf, "w") as out, open(args.out_ids, "w") as ids_out:
        for line, tid, _field in rows:
            if tid in keep:
                out.write(line)
        for tid in sorted(keep):
            ids_out.write(tid + "\n")
    print(f"kept {len(keep)} classcode-u transcripts")


if __name__ == "__main__":
    main()
