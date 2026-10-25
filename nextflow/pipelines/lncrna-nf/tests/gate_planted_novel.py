#!/usr/bin/env python3
"""Scientific gate for the lncrna-nf real-container fixture round.

Consumes the planted-loci manifest (novel.expected.gtf, written by
tests/make_fixtures.py) and the final lncRNA set produced by a real run:

- SENSITIVITY: every planted NONCODING locus (novA/novB) must be recovered
  by >= 1 final-set transcript covering >= 80% of its expected exonic bases
  (locus-based matching: StringTie assigns its own MSTRG ids).
- SPECIFICITY: no planted CODING locus (novC) may contribute ANY transcript
  to the final set.

Usage:
    gate_planted_novel.py --expected reference/novel.expected.gtf \
        --final results-lncrna/final/final_lncRNA.gtf
"""
import argparse
import sys


def attr(field, key):
    for item in field.split(";"):
        item = item.strip()
        if item.startswith(key):
            return item.split(" ", 1)[1].strip('"')
    return None


def parse_gtf(path, need_class=False):
    """tid -> dict(chrom, exons=[(s,e)], cls) from exon rows."""
    tx = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                if len(parts) >= 9 and parts[2] == "transcript" and need_class:
                    tid = attr(parts[8], "transcript_id")
                    if tid is not None:
                        rec = tx.setdefault(tid, {"chrom": parts[0], "exons": [],
                                                  "cls": None})
                        rec["cls"] = attr(parts[8], "planted_class")
                continue
            tid = attr(parts[8], "transcript_id")
            if tid is None:
                continue
            rec = tx.setdefault(tid, {"chrom": parts[0], "exons": [],
                                      "cls": None})
            rec["exons"].append((int(parts[3]), int(parts[4])))
    return tx


def exonic_bases(rec):
    bases = set()
    for s, e in rec["exons"]:
        bases.update(range(s, e + 1))
    return bases


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--expected", required=True,
                    help="novel.expected.gtf (planted-loci manifest)")
    ap.add_argument("--final", required=True, help="final_lncRNA.gtf of the run")
    ap.add_argument("--min-recall-frac", type=float, default=0.8)
    args = ap.parse_args()

    expected = parse_gtf(args.expected, need_class=True)
    final = parse_gtf(args.final)
    final_sets = {tid: exonic_bases(rec) for tid, rec in final.items()}

    failures = []
    print(f"final set: {len(final)} transcripts")
    for tid, rec in sorted(expected.items()):
        if rec["cls"] is None:
            sys.exit(f"{args.expected}: transcript {tid} has no planted_class")
        want = exonic_bases(rec)
        if rec["cls"] == "noncoding":
            best = max((len(want & got) for got in final_sets.values()),
                       default=0)
            frac = best / len(want)
            ok = frac >= args.min_recall_frac
            print(f"SENSITIVITY {tid}: best exonic recall {frac:.2f} "
                  f"(need >= {args.min_recall_frac}) -> {'PASS' if ok else 'FAIL'}")
            if not ok:
                failures.append(f"{tid} not recovered (recall {frac:.2f})")
        elif rec["cls"] == "coding":
            hits = [ftid for ftid, got in final_sets.items() if want & got]
            print(f"SPECIFICITY {tid}: final-set overlap {'NONE' if not hits else hits} "
                  f"-> {'PASS' if not hits else 'FAIL'}")
            if hits:
                failures.append(f"coding locus {tid} leaked into the final set: {hits}")

    if failures:
        print("GATE FAILED:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("GATE PASSED (all planted noncoding recovered, all planted coding excluded)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
