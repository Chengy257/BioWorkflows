#!/usr/bin/env python3
"""Summarize the replicate-aware peak stage (v0.5).

Reads the 6-column sample table, derives the per-replicate / IDR / consensus
peak file locations from the results layout (keep in sync with the path
helpers in workflow/rules/common.smk), and writes one row per group:

    group, seqtype, peak_type, n_treats, mode,
    replicate_peaks (per-replicate counts, "sample=count" joined by ";"),
    final_peaks, retained_fraction (final / mean replicate count)

plus a MultiQC custom-content copy of the same table.

Usage:
    replicate_summary.py --samples samples.csv --results-dir results \
        --out Replicate_summary.tsv --mqc Replicate_summary_mqc.tsv
"""
import argparse
import csv
import os


def read_groups(samples_csv):
    """Minimal re-parse of the sample table (only the columns this summary
    needs; the authoritative validation lives in common.smk)."""
    groups = {}
    with open(samples_csv, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            sid = (row["sample_id"] or "").strip()
            role = (row["role"] or "").strip().lower()
            grp = (row["group"] or "").strip()
            seqtype = (row["seqtype"] or "").strip().lower()
            peak_type = (row["peak_type"] or "").strip().lower()
            if not grp:
                continue
            g = groups.setdefault(
                grp, {"seqtype": seqtype, "peak_type": peak_type, "treat": []})
            if role == "treat" and sid not in g["treat"]:
                g["treat"].append(sid)
    return groups


def count_lines(path):
    if not os.path.exists(path):
        return None
    with open(path, "rb") as fh:
        return sum(1 for _ in fh)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--samples", required=True, help="6-column sample table")
    ap.add_argument("--results-dir", required=True,
                    help="results root (config results_dir, relative to the working directory)")
    ap.add_argument("--out", required=True, help="summary TSV output path")
    ap.add_argument("--mqc", required=True, help="MultiQC custom-content TSV output path")
    args = ap.parse_args()

    rd = args.results_dir.rstrip("/\\")
    groups = read_groups(args.samples)
    rows = []
    for grp, g in groups.items():
        narrow = g["seqtype"] in ("atac", "faire") or g["peak_type"] == "narrow"
        suffix = "narrowPeak" if narrow else "broadPeak"
        rep_counts = []
        for s in g["treat"]:
            n = count_lines(os.path.join(rd, "4.peak", "replicates", grp,
                                         f"{s}_peaks.{suffix}"))
            rep_counts.append((s, n))
        if len(g["treat"]) >= 2:
            if narrow:
                mode = "idr"
                final = count_lines(os.path.join(rd, "4.peak", f"{grp}_IDR_peaks.narrowPeak"))
            else:
                mode = "consensus"
                final = count_lines(os.path.join(rd, "4.peak", f"{grp}_consensus_peaks.broadPeak"))
        else:
            mode = "pooled"
            final = count_lines(os.path.join(rd, "4.peak", f"{grp}_peaks.{suffix}"))
        known = [n for _, n in rep_counts if n]
        mean_rep = sum(known) / len(known) if known else 0.0
        retained = (final / mean_rep) if (final is not None and mean_rep > 0) else None
        rows.append([grp, g["seqtype"], g["peak_type"], len(g["treat"]), mode,
                     ";".join(f"{s}={n if n is not None else 'NA'}" for s, n in rep_counts),
                     final if final is not None else "NA",
                     f"{retained:.4f}" if retained is not None else "NA"])

    header = ["group", "seqtype", "peak_type", "n_treats", "mode",
              "replicate_peaks", "final_peaks", "retained_fraction"]
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="\n", encoding="utf-8") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    with open(args.mqc, "w", newline="\n", encoding="utf-8") as fh:
        fh.write("# id: 'replicate_summary_table'\n")
        fh.write("# section_name: 'Replicate peak reproducibility'\n")
        fh.write("# format: 'tsv'\n")
        fh.write("# plot_type: 'table'\n")
        fh.write("# pconfig: {'id': 'replicate_summary_table', 'title': 'Replicate peaks'}\n")
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    print(f"[replicate_summary] {len(rows)} group(s) summarized -> {args.out}")


if __name__ == "__main__":
    main()
