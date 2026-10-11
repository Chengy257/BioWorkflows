#!/usr/bin/env python3
"""Parse a pfam_scan.pl output table into the set of transcript ids that
carry a domain hit at or below the E-value threshold.

Legacy semantics (rna-seq lncrna.smk): keep the hit rows whose sequence
E-value (column 13) is < the threshold; the transcript id is column 1.
Comment lines (#) are skipped; the `-as` (align) layout is accepted as-is.

Usage:
    pfam_hits.py --pfam pfam_scan.output --evalue 1e-5 --out pfam.hitIDs.txt
"""
import argparse


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pfam", required=True)
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    hits = set()
    with open(args.pfam) as fh, open(args.out, "w") as out:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 13:
                continue
            try:
                evalue = float(parts[12])
            except ValueError:
                continue
            if evalue < args.evalue:
                tid = parts[0]
                if tid not in hits:
                    hits.add(tid)
                    out.write(tid + "\n")
    print(f"pfam hits: {len(hits)} transcripts at E < {args.evalue}")


if __name__ == "__main__":
    main()
