#!/usr/bin/env python3
"""Tabulate per-sample read fate across the sRNA cascade.

Writes a wide TSV (one row per sample) and a MultiQC custom-content
long TSV (sample x stage reads) into the QC directory.
"""
import argparse
import os


def fq_read_count(path):
    with open(path, "rb") as fh:
        lines = sum(1 for _ in fh)
    return lines // 4


def trimmed_reads(report):
    with open(report) as fh:
        for line in fh:
            if "Reads written (passing filters)" in line:
                return int(line.split(":")[1].strip().split()[0].replace(",", ""))
    return 0


def mapped_total(counts_path):
    with open(counts_path) as fh:
        next(fh)
        for line in fh:
            feature, count = line.rstrip("\n").split("\t")
            if feature == "__mapped_total":
                return int(count)
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--indir", required=True, help="results/4.expression root")
    # nargs="*": an empty cascade still writes the zero-class summary.
    ap.add_argument("--classes", nargs="*", required=True)
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--trim-reports", nargs="+", required=True)
    ap.add_argument("--raw-counts", type=int, nargs="+", required=True,
                    help="raw reads per sample (same order as --samples)")
    ap.add_argument("--outdir", required=True, help="results/5.QC")
    ap.add_argument("--genome-configured", default="false",
                    help="true emits the genome unmapped column, false writes NA")
    args = ap.parse_args()

    header = ("sample\traw\ttrimmed"
              + "".join(f"\t{c}_mapped\t{c}_unmapped" for c in args.classes)
              + "\tgenome_unmapped")
    rows = []
    long_rows = []
    for sample, raw, report in zip(args.samples, args.raw_counts, args.trim_reports):
        row = [sample, str(raw), str(trimmed_reads(report))]
        stage_labels = ["raw", "trimmed"]
        stage_values = [raw, trimmed_reads(report)]
        for klass in args.classes:
            mapped = mapped_total(os.path.join(args.indir, klass, f"{sample}_counts.txt"))
            unmapped_path = os.path.join(
                args.indir, "..", "3.align", "filter", klass, f"{sample}_unmapped.fq")
            unmapped = fq_read_count(unmapped_path) if os.path.exists(unmapped_path) else -1
            row += [str(mapped), str(unmapped)]
            stage_labels += [f"{klass}_mapped", f"{klass}_unmapped"]
            stage_values += [mapped, unmapped]
        if args.genome_configured.lower() == "true":
            genome_path = os.path.join(
                args.indir, "..", "3.align", "genome", f"{sample}_unmapped.fq")
            row.append(str(fq_read_count(genome_path)) if os.path.exists(genome_path) else "-1")
        else:
            row.append("NA")
        stage_labels.append("genome_unmapped")
        stage_values.append(row[-1])
        rows.append("\t".join(row))
        for label, value in zip(stage_labels, stage_values):
            long_rows.append((sample, label, value))

    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, "cascade_summary.tsv"), "w") as fh:
        fh.write(header + "\n" + "\n".join(rows) + "\n")
    with open(os.path.join(args.outdir, "cascade_summary_mqc.tsv"), "w") as fh:
        fh.write("# id: 'srna_cascade'\n")
        fh.write("# section_name: 'sRNA cascade read fate'\n")
        fh.write("# plot_type: 'bargraph'\n")
        fh.write("# pconfig:\n")
        fh.write("#     id: 'srna_cascade_bg'\n")
        fh.write("#     title: 'sRNA cascade read fate'\n")
        fh.write("sample\tstage\treads\n")
        for sample, stage, reads in long_rows:
            fh.write(f"{sample}\t{stage}\t{reads}\n")


if __name__ == "__main__":
    main()
