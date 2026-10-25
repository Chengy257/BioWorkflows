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
  - CPAT 3.x <prefix>.ORF_prob.best.tsv: one row per sequence (its best
    ORF); header with a Coding_Probability column, which we treat as
    noncoding below 0.5 (the tool's own recommended cutoff). Sequences in
    <prefix>.no_ORF.txt (no ORF at all) are noncoding per CPAT itself.

Usage:
    coding_potential_intersect.py --cpc2 cpc2.txt --cpat cpat.ORF_prob.best.tsv \
        [--cpat-no-orf cpat.no_ORF.txt] \
        --gtf novel.gtf --min-length 200 --out-ids ids.txt --out-gtf gtf
"""
import argparse
import os
import sys


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


def read_cpat(path, no_orf_path=""):
    """CPAT 3.x best-ORF table: exactly one row per sequence, so a sequence
    is noncoding iff its best-ORF coding probability is < 0.5. A malformed
    table (no coding-probability header) is a hard error, not an empty set."""
    noncoding = set()
    prob_col = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.lstrip("#").split("\t")
            if prob_col is None:
                for i, name in enumerate(parts):
                    if name.strip().lower().startswith("coding_prob"):
                        prob_col = i
                        break
                if prob_col is None:
                    sys.exit(f"{path}: no coding-probability column in the "
                             f"CPAT header: {parts}")
                continue
            if len(parts) <= prob_col:
                sys.exit(f"{path}: short row, expected >= {prob_col + 1} "
                         f"columns: {parts}")
            tid = parts[0].strip()
            try:
                prob = float(parts[prob_col])
            except ValueError:
                sys.exit(f"{path}: non-numeric coding probability for "
                         f"{tid}: {parts[prob_col]}")
            if prob < 0.5:
                noncoding.add(tid)
    # sequences with no ORF at all: CPAT's own verdict is non-coding
    if no_orf_path and os.path.exists(no_orf_path):
        with open(no_orf_path) as fh:
            for line in fh:
                line = line.strip()
                if line:
                    noncoding.add(line)
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
    ap.add_argument("--cpat-no-orf", default="")
    ap.add_argument("--gtf", required=True, help="the classcode-u novel GTF")
    ap.add_argument("--min-length", type=int, default=200)
    ap.add_argument("--out-ids", required=True)
    ap.add_argument("--out-gtf", required=True)
    args = ap.parse_args()

    cpc2 = read_cpc2(args.cpc2)
    cpat = read_cpat(args.cpat, args.cpat_no_orf)
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
