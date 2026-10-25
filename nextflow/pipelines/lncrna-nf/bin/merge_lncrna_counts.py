#!/usr/bin/env python3
"""Merge raw featureCounts per-sample outputs into the lncrna-level
count/TPM/FPKM matrices.

Adapted from the rnaseq-ext merge_featurecounts.py (same raw featureCounts
layout: comment lines start with '#', the header row starts with 'Geneid',
the sample column is the last one) with two lncrna-specific additions: an
optional output filename prefix (the matrices must stay distinct from the
Stage-2 gene-level matrices) and TPM/FPKM computed from the featureCounts
Length column and the per-sample Assigned reads.

Usage:
    python merge_lncrna_counts.py <counts_dir> <id_list.txt> [prefix]

<counts_dir>/ must contain per-sample <sample>.featureCounts.tsv and
<sample>.featureCounts.tsv.summary files. <id_list.txt> lists the expected
sample ids one per line; every id must have both files and no strays are
accepted.

Outputs (written into <counts_dir>):
    <prefix>count.matrix.tsv   gene x sample raw counts (gene order anchored
                               to the first sample)
    <prefix>count_stats.tsv    per-sample assignment statistics
    <prefix>tpm.matrix.tsv     gene x sample TPM
    <prefix>fpkm.matrix.tsv    gene x sample FPKM
"""
import os
import sys


def read_count_file(path):
    """Return (sample, [(gene_id, length, count), ...]) from a raw
    featureCounts file. Layout: Geneid, Chr, Start, End, Strand, Length,
    <sample>; the sample column is the last one."""
    sample = None
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "Geneid":
                sample = parts[-1]
                continue
            if sample is None:
                raise ValueError(f"{path}: data row before the Geneid header")
            rows.append((parts[0], int(parts[5]), int(parts[-1])))
    if sample is None:
        raise ValueError(f"{path}: missing Geneid header row")
    return sample, rows


def header_ok(sample, sid):
    """featureCounts labels the count column with the input BAM basename
    (<id>.markdup.sorted.bam), not the bare sample id; accept both forms but
    nothing that belongs to a different sample."""
    return sample == sid or (sample.startswith(sid + ".") and sample.endswith(".bam"))


def read_summary(path):
    stats = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or parts[0] in ("Status", ""):
                continue
            stats.setdefault(parts[0], parts[1])
    return stats


def main():
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        sys.exit(__doc__)
    counts_dir = sys.argv[1]
    id_list = sys.argv[2]
    prefix = sys.argv[3] if len(sys.argv) == 4 else ""

    ids = [x.strip() for x in open(id_list) if x.strip()]

    per_sample = {}
    gene_order = None
    summaries = {}
    for sid in ids:
        cpath = os.path.join(counts_dir, f"{sid}.featureCounts.tsv")
        spath = os.path.join(counts_dir, f"{sid}.featureCounts.tsv.summary")
        if not os.path.exists(cpath) or not os.path.exists(spath):
            sys.exit(f"missing featureCounts outputs for sample {sid} in {counts_dir}")
        sample, rows = read_count_file(cpath)
        if not header_ok(sample, sid):
            sys.exit(f"{cpath}: header sample column '{sample}' does not match "
                     f"expected '{sid}' or '<{sid}>....bam' (featureCounts "
                     f"labels the column with the input BAM basename)")
        if gene_order is None:
            gene_order = [r[0] for r in rows]
        elif [r[0] for r in rows] != gene_order:
            sys.exit(f"{cpath}: gene order differs from the first sample")
        per_sample[sid] = rows
        summaries[sid] = read_summary(spath)

    # TPM needs the per-sample RPK normalizer; FPKM needs mapped reads.
    rpk_denom = {}
    assigned = {}
    for sid in ids:
        rpk = [c / (length / 1000.0) if length > 0 else 0.0
               for (_g, length, c) in per_sample[sid]]
        rpk_denom[sid] = sum(rpk)
        assigned[sid] = float(summaries[sid].get("Assigned", 0))

    with open(os.path.join(counts_dir, f"{prefix}count.matrix.tsv"), "w") as count_out, \
         open(os.path.join(counts_dir, f"{prefix}tpm.matrix.tsv"), "w") as tpm_out, \
         open(os.path.join(counts_dir, f"{prefix}fpkm.matrix.tsv"), "w") as fpkm_out:
        header = "gene_id\t" + "\t".join(ids) + "\n"
        count_out.write(header)
        tpm_out.write(header)
        fpkm_out.write(header)
        for i, gene in enumerate(gene_order):
            counts = [str(per_sample[sid][i][2]) for sid in ids]
            count_out.write(gene + "\t" + "\t".join(counts) + "\n")
            tpms = []
            fpkms = []
            for sid in ids:
                length, c = per_sample[sid][i][1], float(per_sample[sid][i][2])
                rpk = c / (length / 1000.0) if length > 0 else 0.0
                tpms.append(rpk / rpk_denom[sid] * 1e6 if rpk_denom[sid] > 0 else 0.0)
                fpkms.append(rpk / assigned[sid] * 1e6 if assigned[sid] > 0 else 0.0)
            tpm_out.write(gene + "\t" + "\t".join(f"{v:.6f}" for v in tpms) + "\n")
            fpkm_out.write(gene + "\t" + "\t".join(f"{v:.6f}" for v in fpkms) + "\n")

    with open(os.path.join(counts_dir, f"{prefix}count_stats.tsv"), "w") as out:
        all_keys = []
        for sid in ids:
            for key in summaries[sid]:
                if key not in all_keys:
                    all_keys.append(key)
        out.write("Status\t" + "\t".join(ids) + "\n")
        for key in all_keys:
            out.write(key + "\t" + "\t".join(summaries[sid].get(key, "") for sid in ids) + "\n")

    print(f"wrote matrices for {len(ids)} samples, {len(gene_order)} features into {counts_dir}")


if __name__ == "__main__":
    main()
