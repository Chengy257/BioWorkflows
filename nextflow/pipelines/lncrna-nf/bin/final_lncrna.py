#!/usr/bin/env python3
"""Build the final lncRNA set.

Rule (decision record D4/D5): from the coding-potential-filtered noncoding
GTF, remove every transcript carrying a Pfam or DIAMOND/NR homology hit
(union of the hit-id lists); emit the final lncRNA GTF (sorted by
coordinate), its FASTA via the genome, and the reference-plus-novel merged
GTF (sorted by coordinate) for the expression stage.

Usage:
    final_lncrna.py --gtf noncoding.gtf --reference ref.gtf --genome genome.fa \
        [--pfam 'pfam/*'] [--diamond 'diamond/*'] \
        --out-gtf final_lncRNA.gtf --out-fa final_lncRNA.fa \
        --out-ref ref_withLncRNA.gtf --out-hits homology.hitIDs.txt
"""
import argparse
import glob
import sys


def transcript_id(field):
    for item in field.split(";"):
        item = item.strip()
        if item.startswith("transcript_id"):
            return item.split(" ", 1)[1].strip('"')
    return None


def read_hit_ids(pattern):
    """Collect hit transcript ids from the staged hit files. Both layouts
    occur: one id per line (pfam.hitIDs.txt) or a tabular BLAST table whose
    column 1 is the query id (DIAMOND outfmt 6 'qseqid sseqid')."""
    ids = set()
    for path in sorted(glob.glob(pattern)):
        for line in open(path):
            line = line.rstrip("\n")
            if not line:
                continue
            tid = line.split("\t")[0].strip()
            if tid:
                ids.add(tid)
    return ids


def sort_gtf(lines):
    def key(line):
        parts = line.split("\t")
        return (parts[0], int(parts[3]), 0 if parts[2] == "transcript" else 1)
    return sorted(lines, key=key)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--gtf", required=True, help="coding-potential-filtered noncoding GTF")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--pfam", default="")
    ap.add_argument("--diamond", default="")
    ap.add_argument("--out-gtf", required=True)
    ap.add_argument("--out-fa", required=True)
    ap.add_argument("--out-ref", required=True)
    ap.add_argument("--out-hits", required=True)
    args = ap.parse_args()

    hits = set()
    if args.pfam:
        hits |= read_hit_ids(args.pfam)
    if args.diamond:
        hits |= read_hit_ids(args.diamond)
    with open(args.out_hits, "w") as out:
        for tid in sorted(hits):
            out.write(tid + "\n")

    # filter the noncoding GTF and collect transcript sequences
    kept_lines = []
    kept_tids = []
    seen = set()
    for line in open(args.gtf):
        if line.startswith("#"):
            continue
        tid = transcript_id(line.split("\t")[8]) if len(line.split("\t")) > 8 else None
        if tid is None or tid in hits:
            continue
        kept_lines.append(line)
        if tid not in seen:
            seen.add(tid)
            kept_tids.append(tid)
    with open(args.out_gtf, "w") as out:
        out.writelines(sort_gtf(kept_lines))

    # reference + novel merged GTF, coordinate-sorted
    ref_lines = [l for l in open(args.reference) if not l.startswith("#")]
    with open(args.out_ref, "w") as out:
        out.writelines(sort_gtf(ref_lines + kept_lines))

    # FASTA extraction: two-pass genome read, sequences spliced from exons
    exons = {}
    for line in kept_lines:
        parts = line.rstrip("\n").split("\t")
        if parts[2] != "exon":
            continue
        tid = transcript_id(parts[8])
        exons.setdefault(tid, []).append((parts[0], int(parts[3]) - 1, int(parts[4]), parts[6]))
    genome = {}
    name = None
    chunks = []
    for line in open(args.genome):
        if line.startswith(">"):
            if name is not None:
                genome[name] = "".join(chunks)
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if name is not None:
        genome[name] = "".join(chunks)

    with open(args.out_fa, "w") as out:
        for tid in kept_tids:
            seq_bits = []
            for (chrom, start, end, strand) in sorted(exons.get(tid, []), key=lambda e: e[1]):
                seq = genome.get(chrom, "")
                seq_bits.append(seq[start:end])
            full = "".join(seq_bits)
            if any(e[3] == "-" for e in exons.get(tid, [])):
                full = _revcomp(full)
            out.write(f">{tid}\n")
            for i in range(0, len(full), 60):
                out.write(full[i:i + 60] + "\n")
    print(f"final lncRNA set: {len(kept_tids)} transcripts "
          f"(removed {len(hits)} homology-hit ids)")


def _revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


if __name__ == "__main__":
    sys.exit(main())
