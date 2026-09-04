#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Merge per-sample featureCounts quantification results.

Usage:
    python merge_featurecounts.py <quantification_dir>

<quantification_dir>/ must contain per-sample <sample>.count and <sample>.log files
(produced by scripts/run_featurecounts.R; .count has 5 columns: id/effLength/counts/fpkm/tpm).

Outputs (written to the same directory):
    count.matrix.tsv            gene x sample raw counts matrix (runDESeq2 input)
    GeneExpression_TPM.xls      gene x sample TPM matrix
    GeneExpression_FPKM.xls     gene x sample FPKM matrix
    GeneCount_Assigned_logs.xls featureCounts per-status (Assigned/Unassigned_*) statistics

Replaces the former merge shell script (its dependencies njoin.sh / transposition.sh are not in the repo);
the output format matches the original.
"""
import glob
import os
import sys


def sample_name(path):
    """<sample>.count -> sample name; '-' replaced with '_', consistent with runDESeq2 id handling"""
    base = os.path.basename(path)
    for ext in (".count", ".log"):
        if base.endswith(ext):
            return base[:-len(ext)].replace("-", "_")
    return base.replace("-", "_")


def read_count_file(path):
    """Return [(gene_id, counts, fpkm, tpm), ...]"""
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header[:5] != ["id", "effLength", "counts", "fpkm", "tpm"]:
            sys.exit(f"[ERROR] unexpected header in {path}: {header[:5]}")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            rows.append((parts[0], parts[2], parts[3], parts[4]))
    return rows


def read_log_file(path):
    """featureCounts $stat table -> {status: count}; skips the header row (Status/...)"""
    stats = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or parts[0] in ("Status", ""):
                continue
            stats.setdefault(parts[0], parts[1])
    return stats


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    out_dir = sys.argv[1].rstrip("/\\") or "."

    count_files = sorted(glob.glob(os.path.join(out_dir, "*.count")))
    if not count_files:
        sys.exit(f"[ERROR] no *.count files found in {out_dir}")

    ## ---- Merge count/fpkm/tpm matrices (gene order anchored to the first sample)----
    gene_order = []
    values = {}  # gene_id -> {sample: (counts, fpkm, tpm)}
    for cf in count_files:
        s = sample_name(cf)
        for gid, counts, fpkm, tpm in read_count_file(cf):
            if gid not in values:
                values[gid] = {}
                gene_order.append(gid)
            values[gid][s] = (counts, fpkm, tpm)
    samples = [sample_name(cf) for cf in count_files]

    def write_matrix(fname, idx):
        with open(os.path.join(out_dir, fname), "w") as fh:
            fh.write("id\t" + "\t".join(samples) + "\n")
            for gid in gene_order:
                row = [values[gid].get(s, ("NA", "NA", "NA"))[idx] for s in samples]
                fh.write(gid + "\t" + "\t".join(row) + "\n")

    write_matrix("count.matrix.tsv", 0)
    write_matrix("GeneExpression_FPKM.xls", 1)
    write_matrix("GeneExpression_TPM.xls", 2)

    ## ---- Merge featureCounts status statistics (rows=samples, columns=statuses)----
    log_files = sorted(glob.glob(os.path.join(out_dir, "*.log")))
    status_order = []
    stat_values = {}  # sample -> {status: value}
    for lf in log_files:
        s = sample_name(lf)
        stat_values[s] = read_log_file(lf)
        for k in stat_values[s]:
            if k not in status_order:
                status_order.append(k)
    with open(os.path.join(out_dir, "GeneCount_Assigned_logs.xls"), "w") as fh:
        fh.write("id\t" + "\t".join(status_order) + "\n")
        for s in samples:
            row = [stat_values.get(s, {}).get(k, "NA") for k in status_order]
            fh.write(s + "\t" + "\t".join(row) + "\n")

    print(f"[merge_featurecounts] {len(gene_order)} genes × {len(samples)} samples -> {out_dir}")


if __name__ == "__main__":
    main()
