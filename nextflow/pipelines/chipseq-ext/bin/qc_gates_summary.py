#!/usr/bin/env python3
"""QC gate summary: one PASS/WARN/FAIL row per sample aggregating the QC
metrics available in the adopted nf-core/chipseq run's published outputs.

Metric sources and gates (a metric whose source file is missing or
unparseable renders NA and does not gate):

    mapping_rate        flagstat/{sample}.mLb.clN.sorted.flagstat
                        (samtools flagstat mapped/total)
                        gate: value >= --mapping-rate-min
    dup_rate            picard MarkDuplicates metrics (--dup-metrics glob)
                        (PERCENT_DUPLICATION)
                        gate: value <= --dup-rate-max
    organelle_fraction  organelle summary table row (--organelle glob),
                        computed there with the configured contig patterns
                        gate: value <= --organelle-max (skipped when empty/NA)

Overall per sample: FAIL if any metric FAILs, WARN if no FAIL but at least
one WARN, PASS otherwise; all-NA renders NA. The pipeline never hard-fails on
these gates (informational, mirroring the legacy behavior).

Flagstat count-line formats accepted (samtools <= 1.17 and newer builds):
"N + N mapped" with any of end-of-line, a space, or an opening parenthesis
after the keyword.
"""
import argparse
import glob
import os
import re

FLAGSTAT_TOTAL_RE = re.compile(r"^\s*(\d+)\s+\+\s+\d+\s+in total\b")
FLAGSTAT_MAPPED_RE = re.compile(r"^\s*(\d+)\s+\+\s+\d+\s+mapped(?:\s|\(|$)")

HEADER = ["sample", "mapping_rate", "dup_rate", "organelle_fraction", "gate"]


def parse_flagstat(path):
    """(total, mapped) from a samtools flagstat report; (None, None) when the
    file carries no parseable count lines."""
    total = None
    mapped = None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            m = FLAGSTAT_TOTAL_RE.match(line)
            if m:
                total = int(m.group(1))
                continue
            m = FLAGSTAT_MAPPED_RE.match(line)
            if m:
                mapped = int(m.group(1))
    return total, mapped


def parse_dup_rate(path):
    """PERCENT_DUPLICATION (0-1) from a picard MarkDuplicates metrics file;
    None when unparseable."""
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    # metrics row: LIBRARY header line, then one data row; the percent
    # duplication is the 9th column (index 8) after the LIBRARY column
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.startswith("LIBRARY\t"):
            cols = lines[i + 1].split("\t") if i + 1 < len(lines) else []
            if len(cols) > 8:
                try:
                    return float(cols[8])
                except ValueError:
                    return None
    return None


def gate(value, direction, threshold):
    """PASS/WARN/FAIL against a threshold; direction is 'min' or 'max'."""
    if value is None or threshold in (None, ""):
        return None
    try:
        threshold = float(threshold)
    except (TypeError, ValueError):
        return None
    if direction == "min":
        if value >= threshold:
            return "PASS"
        return "FAIL"
    if value <= threshold:
        return "PASS"
    return "FAIL"


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--mapping-rate-min", default="0.70")
    ap.add_argument("--dup-rate-max", default="0.50")
    ap.add_argument("--organelle-max", default="",
                    help="organelle fraction gate; empty = informational only")
    ap.add_argument("--dup-metrics", default=None,
                    help="glob for picard duplicate metrics files")
    ap.add_argument("--dup-suffix", default="",
                    help="file-name suffix after the sample id in the metrics "
                         "basenames (stripped before joining)")
    ap.add_argument("--flagstat-suffix", default="",
                    help="file-name suffix after the sample id in the flagstat "
                         "basenames (stripped to derive the sample column)")
    ap.add_argument("--organelle", default=None,
                    help="glob for the organelle summary TSV (joined by sample)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--mqc", required=True)
    ap.add_argument("flagstats", nargs="+", help="per-sample flagstat files")
    args = ap.parse_args()

    def strip_suffix(name, suffix):
        if suffix and name.endswith(suffix):
            return name[: -len(suffix)]
        return re.sub(r"(\.flagstat|\.flagstat\.txt|\.txt)$", "", name)

    # sample -> dup rate; both sides are joined on the exact sample id after
    # stripping the pipeline-declared file-name suffixes
    dup_by_sample = {}
    if args.dup_metrics:
        for path in glob.glob(args.dup_metrics):
            name = os.path.basename(path)
            dup_by_sample[strip_suffix(name, args.dup_suffix)] = parse_dup_rate(path)

    # sample -> organelle fraction (column 4 of the summary table)
    org_by_sample = {}
    if args.organelle:
        for path in glob.glob(args.organelle):
            with open(path, encoding="utf-8") as fh:
                header = fh.readline().rstrip("\n").split("\t")
                idx = {name: i for i, name in enumerate(header)}
                for line in fh:
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) > idx.get("organelle_fraction", 3):
                        try:
                            org_by_sample[cols[0]] = float(cols[idx["organelle_fraction"]])
                        except ValueError:
                            org_by_sample[cols[0]] = None

    def dup_for(sample):
        return dup_by_sample.get(sample)

    def org_for(sample):
        return org_by_sample.get(sample)

    rows = []
    for path in args.flagstats:
        name = os.path.basename(path)
        sample = strip_suffix(name, args.flagstat_suffix)
        total, mapped = parse_flagstat(path)
        mapping_rate = (mapped / total) if (total and mapped is not None) else None
        dup_rate = dup_for(sample)
        org = org_for(sample)

        gates = [
            gate(mapping_rate, "min", args.mapping_rate_min),
            gate(dup_rate, "max", args.dup_rate_max),
            gate(org, "max", args.organelle_max),
        ]
        if all(g is None for g in gates):
            overall = "NA"
        elif "FAIL" in gates:
            overall = "FAIL"
        elif "WARN" in gates:
            overall = "WARN"
        else:
            overall = "PASS"

        def fmt(v):
            return f"{v:.4f}" if v is not None else "NA"

        rows.append([sample, fmt(mapping_rate), fmt(dup_rate), fmt(org), overall])

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    for path in (args.out, args.mqc):
        with open(path, "w", newline="\n", encoding="utf-8") as fh:
            if path == args.mqc:
                fh.write("# id: 'qc_gates_table'\n")
                fh.write("# section_name: 'QC gates (PASS/WARN/FAIL per sample)'\n")
                fh.write("# format: 'tsv'\n")
                fh.write("# plot_type: 'table'\n")
                fh.write("# pconfig: {'id': 'qc_gates_table', 'title': 'QC gates'}\n")
            fh.write("\t".join(HEADER) + "\n")
            for row in rows:
                fh.write("\t".join(str(x) for x in row) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
