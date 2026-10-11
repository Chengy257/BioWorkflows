#!/usr/bin/env python3
"""Intersect the coding-potential filters into the noncoding candidate set.

Decision-record D3 rule (replacing the legacy CPC2 AND CNCI AND length):
a novel transcript is noncoding iff
  - CPC2 calls it noncoding, AND
  - CPAT calls it noncoding, AND
  - its exonic length is >= --min-length.

Input formats:
  - CPC2 tabular: header line then rows; the transcript id is column 1 and
    the verdict column 8 ("coding"/"noncoding") - the legacy awk `$8=="noncoding"`.
  - CPAT tabular (with -g <run id>): a header line containing sequence
    columns then rows; column 1 = sequence id, the last column = the
    coding probability ("coding_prob"/"ORF" layout varies) - we accept the
    documented CPAT 3.x layout where column 7 is 'coding_prob' and treat a
    probability < 0.5 as noncoding (the tool's own recommended cutoff).

Usage:
    coding_potential_intersect.py --cpc2 cpc2.txt --cpat cpat.out.tsv \
        --gtf novel.gtf --min-length 200 --out-ids ids.txt --out-gtf gtf
"""
import argparse


def read_cpc2(path):
    noncoding = set()
    with open(path) as fh:
        header_seen = False
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if not header_seen and parts[0].lower().startswith("id"):
                header_seen = True
                continue
            if len(parts) >= 8 and parts[7].strip().lower() == "noncoding":
                noncoding.add(parts[0].strip())
    return noncoding


def read_cpat(path):
    noncoding = set()
    with open(path) as fh:
        prob_col = None
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if prob_col is None:
                for i, name in enumerate(parts):
                    if name.lower().startswith("coding_prob"):
                        prob_col = i
                        break
                continue
            if prob_col is not None and len(parts) > prob_col:
                try:
                    if float(parts[prob_col]) < 0.5:
                        noncoding.add(parts[0].strip())
                except ValueError:
                    continue
    return noncoding


def exonic_lengths(gtf_path):
    lengths = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue
            tid = None
            for item in parts[8].split(";"):
                item = item.strip()
                if item.startswith("transcript_id"):
                    tid = item.split(" ", 1)[1].strip('"')
                    break
            if tid is None:
                continue
            lengths.setdefault(tid, 0)
            lengths[tid] += int(parts[4]) - int(parts[3]) + 1
    return lengths


def filter_gtf(gtf_path, keep_ids, out_path):
    kept = set()
    with open(gtf_path) as fh, open(out_path, "w") as out:
        for line in fh:
            if line.startswith("#"):
                continue
            for item in line.split("\t")[8].split(";"):
                item = item.strip()
                if item.startswith("transcript_id"):
                    tid = item.split(" ", 1)[1].strip('"')
                    if tid in keep_ids:
                        out.write(line)
                        kept.add(tid)
                    break
    return kept


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cpc2", required=True)
    ap.add_argument("--cpat", required=True)
    ap.add_argument("--gtf", required=True, help="the classcode-u novel GTF")
    ap.add_argument("--min-length", type=int, default=200)
    ap.add_argument("--out-ids", required=True)
    ap.add_argument("--out-gtf", required=True)
    args = ap.parse_args()

    cpc2 = read_cpc2(args.cpc2)
    cpat = read_cpat(args.cpat)
    lengths = exonic_lengths(args.gtf)

    noncoding = cpc2 & cpat
    noncoding = {tid for tid in noncoding
                 if lengths.get(tid, 0) >= args.min_length}

    kept = filter_gtf(args.gtf, noncoding, args.out_gtf)
    with open(args.out_ids, "w") as out:
        for tid in sorted(kept):
            out.write(tid + "\n")
    print(f"noncoding set: {len(kept)} transcripts "
          f"(cpc2 {len(cpc2)}, cpat {len(cpat)}, length >= {args.min_length})")


if __name__ == "__main__":
    main()
