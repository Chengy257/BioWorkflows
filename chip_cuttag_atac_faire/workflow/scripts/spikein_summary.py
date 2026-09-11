#!/usr/bin/env python3
"""Spike-in normalization summary: one row per sample from the per-sample
samtools flagstat reports of the spike-in alignment stage (the spike_in
stage).

Per sample: spike_total (flagstat "in total" line), spike_mapped (flagstat
"mapped" line), spike_rate = mapped / total, and scale_factor =
1e6 / max(spike_mapped, 1) — the factor the bigWig signal tracks are
multiplied by when spike_in.scale_bigwigs is enabled (per-sample tracks use
their own factor, group tracks the mean over their treat samples). A sample
whose flagstat is missing or unparseable renders an all-NA row.

Count-line formats (keep in sync with the compiled copies in
workflow/rules/common.smk and workflow/scripts/gates_summary.py): the
flagstat mapped summary line is "N + N mapped (x% : N/A)" on samtools <= 1.17
and "N + N mapped" or "N + N mapped: x%" on newer builds — accept a bare
end-of-line, a space, or an opening parenthesis after the keyword. The
"primary mapped" line carries text between the counts and the keyword, and
the "with itself and mate mapped" line carries text before it, so neither
matches.

Usage:
    spikein_summary.py --samples s1 s2 --spike-dir results/5.QC/spike_in \
        --name lambda --out Spikein_summary.tsv --mqc Spikein_summary_mqc.tsv
"""
import argparse
import os
import re

FLAGSTAT_TOTAL_RE = re.compile(r"^\s*(\d+)\s+\+\s+\d+\s+in total\b")
FLAGSTAT_MAPPED_RE = re.compile(r"^\s*(\d+)\s+\+\s+\d+\s+mapped(?:\s|\(|$)")

HEADER = ["sample", "spike_total", "spike_mapped", "spike_rate", "scale_factor"]


def parse_flagstat_counts(path):
    """(total, mapped) from a samtools flagstat report; (None, None) when the
    file is missing or carries no parseable count lines."""
    total = None
    mapped = None
    try:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                m = FLAGSTAT_TOTAL_RE.match(line)
                if m:
                    total = int(m.group(1))
                    continue
                m = FLAGSTAT_MAPPED_RE.match(line)
                if m:
                    mapped = int(m.group(1))
    except OSError:
        return None, None
    return total, mapped


def sample_row(sample, spike_dir):
    """One summary row; a missing/unparseable flagstat yields NA cells."""
    total, mapped = parse_flagstat_counts(
        os.path.join(spike_dir, f"{sample}_flagstat.txt"))
    if total is None or mapped is None:
        return [sample, "NA", "NA", "NA", "NA"]
    rate = (mapped / total) if total else 0.0
    factor = 1e6 / max(mapped, 1)
    return [sample, total, mapped, f"{rate:.4f}", f"{factor:.4f}"]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--samples", nargs="+", required=True,
                    help="sample ids in output order")
    ap.add_argument("--spike-dir", required=True,
                    help="results/5.QC/spike_in directory (per-sample flagstat inputs)")
    ap.add_argument("--name", required=True,
                    help="spike-in label used in the MultiQC section title")
    ap.add_argument("--out", required=True, help="summary TSV output path")
    ap.add_argument("--mqc", required=True, help="MultiQC custom-content TSV output path")
    args = ap.parse_args()

    rows = [sample_row(s, args.spike_dir) for s in args.samples]

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="\n", encoding="utf-8") as fh:
        fh.write("\t".join(HEADER) + "\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    label = args.name.replace("'", "")
    with open(args.mqc, "w", newline="\n", encoding="utf-8") as fh:
        fh.write("# id: 'spikein_summary_table'\n")
        fh.write(f"# section_name: 'Spike-in normalization ({label})'\n")
        fh.write("# format: 'tsv'\n")
        fh.write("# plot_type: 'table'\n")
        fh.write(f"# pconfig: {{'id': 'spikein_summary_table', 'title': 'Spike-in ({label})'}}\n")
        fh.write("\t".join(HEADER) + "\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    print(f"[spikein_summary] {len(rows)} sample(s) summarized -> {args.out}")


if __name__ == "__main__":
    main()
