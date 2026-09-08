#!/usr/bin/env python3
"""Build a DiffBind sample sheet for one contrast (v0.5 Phase 3).

Reads the workflow sample table, derives per-treat BAM paths and per-sample
peak files from the results layout (replicate calls when --replicates is
passed, pooled group peaks otherwise — keep in sync with common.smk), and
writes the TSV DiffBind's dba(sampleSheet=...) expects:

    SampleID, Condition, Replicate, bamReads, bamControl, PeakFile, Batch

Condition comes from the optional `condition` column (falling back to the
group name); the optional `batch` column becomes the blocking factor. With
--use-controls, groups carrying exactly one control attach it as bamControl
(others get an empty entry — DiffBind tolerates a missing control).

Usage:
    diffbind_sheet.py --samples samples.csv --results-dir results \
        --groups gA gB [--replicates] [--use-controls] --out samplesheet.tsv
"""
import argparse
import csv
import os


def read_rows(samples_csv):
    with open(samples_csv, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--samples", required=True, help="workflow sample table")
    ap.add_argument("--results-dir", required=True,
                    help="results root (config results_dir)")
    ap.add_argument("--groups", required=True, nargs=2,
                    help="the two contrast arm group names")
    ap.add_argument("--replicates", action="store_true",
                    help="use per-replicate peak files (peak.replicate enabled)")
    ap.add_argument("--use-controls", action="store_true", dest="use_controls",
                    help="attach each group's single control as bamControl")
    ap.add_argument("--out", required=True, help="sample sheet TSV output path")
    args = ap.parse_args()

    rd = args.results_dir.rstrip("/\\")
    rows = read_rows(args.samples)
    groups = {}
    for row in rows:
        grp = (row["group"] or "").strip()
        g = groups.setdefault(grp, {"treats": [], "controls": [],
                                    "condition": "", "seqtype": "",
                                    "peak_type": ""})
        g["seqtype"] = (row["seqtype"] or "").strip().lower() or g["seqtype"]
        g["peak_type"] = (row["peak_type"] or "").strip().lower() or g["peak_type"]
        entry = {
            "sample": (row["sample_id"] or "").strip(),
            "condition": (row.get("condition") or "").strip(),
            "batch": (row.get("batch") or "").strip(),
        }
        if (row["role"] or "").strip().lower() == "treat":
            g["treats"].append(entry)
            if entry["condition"]:
                g["condition"] = g["condition"] or entry["condition"]
        else:
            g["controls"].append(entry["sample"])

    out_rows = []
    for grp in args.groups:
        g = groups.get(grp)
        if g is None or not g["treats"]:
            raise SystemExit(f"group {grp!r} not found or has no treats in {args.samples}")
        condition = g["condition"] or grp
        suffix = "narrowPeak" if (g["seqtype"] in ("atac", "faire")
                                  or g["peak_type"] == "narrow") else "broadPeak"
        control = ""
        if args.use_controls and len(g["controls"]) == 1:
            control = os.path.join(rd, "3.align", "bowtie2",
                                   f"{g['controls'][0]}_rmdup.bam")
            if not os.path.exists(control):
                control = os.path.join(rd, "3.align", "bowtie2",
                                       f"{g['controls'][0]}_sorted.bam")
        for rep, entry in enumerate(g["treats"], start=1):
            sid = entry["sample"]
            if args.replicates:
                peaks = os.path.join(rd, "4.peak", "replicates", grp,
                                     f"{sid}_peaks.{suffix}")
            else:
                peaks = os.path.join(rd, "4.peak", f"{grp}_peaks.{suffix}")
            out_rows.append([
                sid, condition, str(rep),
                os.path.join(rd, "3.align", "bowtie2", f"{sid}_rmdup.bam")
                if os.path.exists(os.path.join(rd, "3.align", "bowtie2", f"{sid}_rmdup.bam"))
                else os.path.join(rd, "3.align", "bowtie2", f"{sid}_sorted.bam"),
                control, peaks, entry["batch"],
            ])

    header = ["SampleID", "Condition", "Replicate", "bamReads",
              "bamControl", "PeakFile", "Batch"]
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="\n", encoding="utf-8") as fh:
        fh.write("\t".join(header) + "\n")
        for row in out_rows:
            fh.write("\t".join(row) + "\n")
    print(f"[diffbind_sheet] {len(out_rows)} samples ({args.groups[0]} vs "
          f"{args.groups[1]}) -> {args.out}")


if __name__ == "__main__":
    main()
