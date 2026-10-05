#!/usr/bin/env python3
"""Build a DiffBind sample sheet for one contrast.

Reads the chipseq-ext project sample sheet (CSV with sample_id, group,
condition, control, batch), maps every treat sample of the two contrast
conditions to its per-sample BAM and peak file (staged under bams/ and
peaks/), and writes the TSV that DiffBind's dba(sampleSheet=...) expects:

    SampleID, Condition, Replicate, bamReads, bamControl, Peaks, PeakCaller,
    Batch

(These are DiffBind's own sample-sheet column names; a PeakFile column is
silently ignored by dba(sampleSheet=...), and PeakCaller must be
narrowpeak|broadpeak per file suffix to encode the score-column convention.)

Paths inside the sheet are RELATIVE (bams/<file>, peaks/<file>): DiffBind
resolves them against the sheet's own directory, and the pipeline stages the
BAMs and peaks under the same layout in the analysis task.

Condition comes from the sample sheet's condition column (falling back to the
group name). The optional batch column becomes the blocking factor; with
--use-controls, groups carrying exactly one control attach that control's BAM
as bamControl (others get an empty entry - DiffBind tolerates a missing
control).

Usage:
    diffbind_sheet.py --contrast WT_vs_mut --sheet samplesheet.csv \
        --bam-suffix '.mLb.clN.sorted.bam' --peak-suffix '_peaks.narrowPeak' \
        --peak-type narrowpeak --use-controls F --out sheet.tsv \
        bams/* peaks/*
"""
import argparse
import csv
import os
import sys

COLUMNS = ["SampleID", "Condition", "Replicate", "bamReads", "bamControl",
           "Peaks", "PeakCaller", "Batch"]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--contrast", required=True,
                    help="'<conditionA>_vs_<conditionB>'")
    ap.add_argument("--sheet", required=True,
                    help="project sample sheet (CSV)")
    ap.add_argument("--bam-suffix", required=True)
    ap.add_argument("--peak-suffix", required=True)
    ap.add_argument("--peak-type", required=True, choices=["narrowpeak", "broadpeak"])
    ap.add_argument("--use-controls", required=True, choices=["T", "F"])
    ap.add_argument("--out", required=True)
    ap.add_argument("files", nargs="+",
                    help="staged BAM files followed by peak files")
    args = ap.parse_args()

    cond_a, cond_b = args.contrast.split("_vs_", 1)
    peak_ext = os.path.splitext(args.peak_suffix)[1] or args.peak_suffix

    # staged files are split by directory (bams/* first, peaks/* second);
    # classify by suffix so the call order does not matter
    bam_by_sample = {}
    peak_by_sample = {}
    for path in args.files:
        base = os.path.basename(path)
        if base.endswith(args.bam_suffix):
            bam_by_sample[base[: -len(args.bam_suffix)]] = path
        elif base.endswith(args.peak_suffix):
            peak_by_sample[base[: -len(args.peak_suffix)]] = path

    with open(args.sheet, encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh) if (r.get("sample_id") or "").strip()]

    meta = {}
    group_control = {}
    for row in rows:
        sid = row["sample_id"].strip()
        group = (row.get("group") or "").strip()
        meta[sid] = {
            "id": sid,
            "group": group,
            "condition": (row.get("condition") or group).strip(),
            "control": (row.get("control") or "").strip(),
            "batch": (row.get("batch") or "").strip(),
        }
        if (row.get("control") or "").strip():
            group_control.setdefault(group, row["control"].strip())

    # group -> control id; the control sample may live in its own group
    ctl_for_group = {}
    for group, ctl in group_control.items():
        ctl_for_group[group] = ctl

    out_rows = []
    missing = []
    replicate_counter = {}
    for row in rows:
        sid = row["sample_id"].strip()
        m = meta[sid]
        condition = m["condition"]
        if condition not in (cond_a, cond_b):
            continue
        if m["control"] and group_control.get(m["group"]) == sid:
            # this sample IS the group's control: never a DiffBind sample row
            continue
        bam = bam_by_sample.get(sid)
        peak = peak_by_sample.get(sid)
        if not bam or not peak:
            missing.append(sid)
            continue
        replicate_counter[condition] = replicate_counter.get(condition, 0) + 1
        bam_control = ""
        if args.use_controls == "T":
            ctl_id = ctl_for_group.get(m["group"], "")
            if ctl_id and ctl_id in bam_by_sample:
                bam_control = os.path.relpath(bam_by_sample[ctl_id], ".")
        out_rows.append({
            "SampleID": sid,
            "Condition": condition,
            "Replicate": replicate_counter[condition],
            "bamReads": os.path.relpath(bam, "."),
            "bamControl": bam_control,
            "Peaks": os.path.relpath(peak, "."),
            "PeakCaller": args.peak_type,
            "Batch": m["batch"],
        })

    if missing:
        print(f"diffbind_sheet: missing staged BAM or peak files for: "
              f"{', '.join(sorted(set(missing)))}", file=sys.stderr)
        return 1
    if not out_rows:
        print(f"diffbind_sheet: no samples match conditions "
              f"{cond_a!r} / {cond_b!r}", file=sys.stderr)
        return 1

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="\n", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(out_rows)
    print(f"diffbind_sheet: {len(out_rows)} samples for contrast {args.contrast}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
