#!/usr/bin/env python3
"""Merge raw featureCounts per-sample outputs into gene-count matrices.

Adapted from the legacy rna-seq merge_featurecounts.py for the RAW
nf-core featureCounts module outputs (the legacy version consumed a
pre-parsed 5-column R intermediate).

Usage:
    python merge_featurecounts.py <quant_dir> <id_list.txt>

<quant_dir>/ must contain per-sample <sample>.featureCounts.tsv and
<sample>.featureCounts.tsv.summary files (raw featureCounts format: comment
lines start with '#'; the header row starts with 'Geneid' and the sample
column is the last one), as staged by the subread/featurecounts module.

<id_list.txt> lists the expected sample ids one per line (the adopted sample
sheet order); every id must have a .featureCounts.txt file and no strays are
accepted.

Outputs (written to <quant_dir>):
    count.matrix.tsv   gene x sample raw counts matrix (gene ids from the
                       first sample anchor the row order)
    count_stats.tsv    per-sample featureCounts assignment statistics
                       (Status table from the .summary files)
"""
import os
import sys


def read_count_file(path):
    """Return (sample, [(gene_id, counts), ...]) from a raw featureCounts file.

    Single-sample layout: Geneid, Chr, Start, End, Strand, Length, <sample>;
    the sample column is the last one."""
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
            rows.append((parts[0], parts[-1]))
    if sample is None:
        raise ValueError(f"{path}: missing Geneid header row")
    return sample, rows


def read_summary(path):
    stats = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or parts[0] in ("Status", ""):
                continue
            stats.setdefault(parts[0], parts[1])
    return stats


def header_ok(sample, sid):
    """featureCounts labels the count column with the input BAM basename
    (<id>.markdup.sorted.bam), not the bare sample id; accept both forms but
    nothing that belongs to a different sample."""
    return sample == sid or (sample.startswith(sid + ".") and sample.endswith(".bam"))


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    out_dir = sys.argv[1].rstrip("/\\") or "."
    with open(sys.argv[2]) as fh:
        samples = [l.strip() for l in fh if l.strip()]
    if not samples:
        sys.exit("[ERROR] empty sample id list")

    gene_order = []
    values = {}  # gene_id -> {sample: counts}
    for s in samples:
        path = os.path.join(out_dir, f"{s}.featureCounts.tsv")
        file_samples, rows = read_count_file(path)
        if not header_ok(file_samples, s):
            sys.exit(f"[ERROR] {path}: sample column '{file_samples}' does not match "
                     f"expected id '{s}' or '<{s}>....bam' (featureCounts labels "
                     f"the column with the input BAM basename)")
        for gid, counts in rows:
            if gid not in values:
                values[gid] = {}
                gene_order.append(gid)
            values[gid][s] = counts

    with open(os.path.join(out_dir, "count.matrix.tsv"), "w") as fh:
        fh.write("gene_id\t" + "\t".join(samples) + "\n")
        for gid in gene_order:
            fh.write(gid + "\t" + "\t".join(values[gid].get(s, "NA") for s in samples) + "\n")

    status_order = []
    stat_values = {}
    for s in samples:
        stat_values[s] = read_summary(os.path.join(out_dir, f"{s}.featureCounts.tsv.summary"))
        for k in stat_values[s]:
            if k not in status_order:
                status_order.append(k)
    with open(os.path.join(out_dir, "count_stats.tsv"), "w") as fh:
        fh.write("sample\t" + "\t".join(status_order) + "\n")
        for s in samples:
            fh.write(s + "\t" + "\t".join(stat_values[s].get(k, "NA") for k in status_order) + "\n")

    print(f"[merge_featurecounts] {len(gene_order)} genes x {len(samples)} samples -> {out_dir}")


if __name__ == "__main__":
    main()
