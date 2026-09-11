#!/usr/bin/env python3
"""QC gate summary: one PASS/WARN/FAIL row per sample aggregating the
workflow's existing QC metrics (the qc.gates stage).

Metric sources and gates (a metric whose source stage is off, or whose file
is missing/unparseable, renders NA and does not gate):

    mapping_rate        {gates_dir}/{sample}_flagstat.txt (samtools flagstat
                        mapped/total)                gate: value >= mapping_rate_min
    dup_rate            {align_dir}/{sample}_dup_metrics.txt (picard
                        PERCENT_DUPLICATION)         gate: value <= dup_rate_max
    frip                {qc_dir}/frip/{group}__{sample}.frip.tsv
                                                     gate: value >= frip_min
    nsc / rsc           {qc_dir}/spp/{sample}_NSC.txt / _RSC.txt
                                                     gate: value >= nsc_min / rsc_min
    tss_enrichment      {qc_dir}/tss/{sample}_TSSE.txt
                                                     gate: value >= tss_min
    organelle_fraction  {qc_dir}/organelle/Organelle_summary.tsv (the sample's
                        row; computed there with the configured contig patterns)
                                                     gate: value <= organelle_max

Overall per sample: FAIL if any metric FAILs, WARN if no FAIL but at least
one NA, PASS otherwise. Strictly informational — callers never abort the
pipeline on a gate. The sample->group mapping for the FRiP filenames comes
from the sample table, baked into the invocation as sample=group tokens.

Usage:
    gates_summary.py --samples s1=g1 s2=g1 ctl=g1 \
        --align-dir results/3.align/bowtie2 --qc-dir results/5.QC \
        --gates-dir results/5.QC/gates \
        --out results/5.QC/gates/gate_summary.tsv \
        --mqc results/5.QC/gates/gate_summary_mqc.tsv \
        --thresholds mapping_rate_min=0.7 dup_rate_max=0.5 frip_min=0.01 \
            nsc_min=1.05 rsc_min=0.8 tss_min=6.0 organelle_max=0.2
"""
import argparse
import os
import re

FLAGSTAT_TOTAL_RE = re.compile(r"^\s*(\d+)\s+\+\s+\d+\s+in total\b")
# The mapped summary line is "N + N mapped (x% : N/A)" on samtools <= 1.17
# and "N + N mapped" or "N + N mapped: x%" on newer builds — accept a bare
# end-of-line, a space, or an opening parenthesis after the keyword. The
# "primary mapped" line carries text between the counts and the keyword, and
# the "with itself and mate mapped" line carries text before it, so neither
# matches.
FLAGSTAT_MAPPED_RE = re.compile(r"^\s*(\d+)\s+\+\s+\d+\s+mapped(?:\s|\(|$)")

# Built-in defaults (keep in sync with qc.gates.thresholds in
# workflow/rules/common.smk).
DEFAULT_THRESHOLDS = {
    "mapping_rate_min": 0.70,
    "dup_rate_max": 0.50,
    "frip_min": 0.01,
    "nsc_min": 1.05,
    "rsc_min": 0.8,
    "tss_min": 6.0,
    "organelle_max": 0.20,
}

# metric column -> (threshold key, comparison): "gte" fails below the minimum,
# "lte" fails above the maximum.
METRICS = (
    ("mapping_rate", "mapping_rate_min", "gte"),
    ("dup_rate", "dup_rate_max", "lte"),
    ("frip", "frip_min", "gte"),
    ("nsc", "nsc_min", "gte"),
    ("rsc", "rsc_min", "gte"),
    ("tss_enrichment", "tss_min", "gte"),
    ("organelle_fraction", "organelle_max", "lte"),
)

HEADER = ["sample", "group", "mapping_rate", "dup_rate", "frip", "nsc", "rsc",
          "tss_enrichment", "organelle_fraction", "failed", "na", "gate"]


def _to_float(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def _first_token(path):
    """First token of the first non-empty line, or None (single-value files)."""
    try:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                parts = line.split()
                if parts:
                    return parts[0]
    except OSError:
        return None
    return None


def parse_mapping_rate(path):
    """mapped/total from a samtools flagstat report; None when absent."""
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
        return None
    if total and mapped is not None:
        return mapped / total
    return None


def parse_dup_rate(path):
    """PERCENT_DUPLICATION (first data row) from a picard MarkDuplicates
    metrics file; None when dedup was skipped for the sample."""
    try:
        with open(path, encoding="utf-8") as fh:
            lines = fh.read().splitlines()
    except OSError:
        return None
    for i, line in enumerate(lines):
        cols = [c.strip() for c in line.split("\t")]
        if cols and cols[0] == "LIBRARY" and "PERCENT_DUPLICATION" in cols:
            idx = cols.index("PERCENT_DUPLICATION")
            for data in lines[i + 1:]:
                if not data.strip():
                    continue
                dcols = data.split("\t")
                if len(dcols) <= idx:
                    return None
                return _to_float(dcols[idx])
    return None


def parse_frip(path, sample):
    """FRiP column from the per-(group, sample) frip.tsv (header + one row)."""
    try:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if cols and cols[0] == sample and len(cols) >= 5:
                    return _to_float(cols[-1])
    except OSError:
        return None
    return None


def parse_tss_score(path, sample):
    """TSSE column from the sample's row in {sample}_TSSE.txt (a literal NA
    cell renders as None)."""
    try:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if cols and cols[0] == sample and len(cols) >= 2:
                    return _to_float(cols[1])
    except OSError:
        return None
    return None


def parse_organelle_fraction(path, sample):
    """organelle_fraction column from the project-wide Organelle_summary.tsv
    (one row per sample; computed there with the configured patterns)."""
    try:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if cols and cols[0] == sample and len(cols) >= 4:
                    return _to_float(cols[3])
    except OSError:
        return None
    return None


def collect_values(sample, group, align_dir, qc_dir, gates_dir):
    """Resolve every metric from its canonical filename; a missing source
    (stage off or file absent) yields None -> NA."""
    return {
        "mapping_rate": parse_mapping_rate(
            os.path.join(gates_dir, f"{sample}_flagstat.txt")),
        "dup_rate": parse_dup_rate(
            os.path.join(align_dir, f"{sample}_dup_metrics.txt")),
        "frip": parse_frip(
            os.path.join(qc_dir, "frip", f"{group}__{sample}.frip.tsv"), sample)
        if group else None,
        "nsc": _to_float(_first_token(
            os.path.join(qc_dir, "spp", f"{sample}_NSC.txt"))),
        "rsc": _to_float(_first_token(
            os.path.join(qc_dir, "spp", f"{sample}_RSC.txt"))),
        "tss_enrichment": parse_tss_score(
            os.path.join(qc_dir, "tss", f"{sample}_TSSE.txt"), sample),
        "organelle_fraction": parse_organelle_fraction(
            os.path.join(qc_dir, "organelle", "Organelle_summary.tsv"), sample),
    }


def evaluate(value, threshold, comparison):
    """PASS/FAIL for one metric against its threshold (value is not None)."""
    if comparison == "lte":
        return "PASS" if value <= threshold else "FAIL"
    return "PASS" if value >= threshold else "FAIL"


def build_row(sample, group, values, thresholds):
    """One output row: the metric values, the failed/NA metric lists, and the
    overall gate (FAIL on any FAIL, WARN when only NAs remain, PASS else)."""
    cells = [sample, group or "-"]
    failed = []
    na = []
    for name, key, comparison in METRICS:
        value = values[name]
        cells.append(f"{value:.4f}" if value is not None else "NA")
        if value is None:
            na.append(name)
        elif evaluate(value, thresholds[key], comparison) == "FAIL":
            failed.append(name)
    gate = "FAIL" if failed else ("WARN" if na else "PASS")
    cells.append(",".join(failed) or "-")
    cells.append(",".join(na) or "-")
    cells.append(gate)
    return cells


def parse_sample_tokens(tokens):
    """'sample' or 'sample=group' tokens -> ordered (sample, group) pairs
    (a plain name leaves FRiP as NA)."""
    pairs = []
    for token in tokens:
        sample, _, group = token.partition("=")
        pairs.append((sample, group or None))
    return pairs


def parse_thresholds(pairs):
    """key=value CLI pairs -> {key: float} (unknown keys are reported later)."""
    parsed = {}
    for pair in pairs:
        key, sep, value = pair.partition("=")
        if not sep:
            sys_exit(f"--thresholds expects key=value pairs, got {pair!r}")
        number = _to_float(value)
        if number is None:
            sys_exit(f"--thresholds value for {key!r} is not numeric: {value!r}")
        parsed[key] = number
    return parsed


def sys_exit(message):
    raise SystemExit(message)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--samples", nargs="+", required=True,
                    help="sample tokens: 'sample=group' (the group resolves the "
                         "FRiP filename) or plain 'sample' (FRiP stays NA)")
    ap.add_argument("--align-dir", required=True,
                    help="results/3.align/bowtie2 directory (picard dup metrics)")
    ap.add_argument("--qc-dir", required=True,
                    help="results/5.QC directory (frip/spp/tss/organelle sources)")
    ap.add_argument("--gates-dir", required=True,
                    help="results/5.QC/gates directory (flagstat inputs)")
    ap.add_argument("--out", required=True, help="summary TSV output path")
    ap.add_argument("--mqc", required=True, help="MultiQC custom-content TSV output path")
    ap.add_argument("--thresholds", nargs="+", default=[], metavar="KEY=VALUE",
                    help="gate thresholds overriding the built-in defaults, e.g. "
                         "mapping_rate_min=0.7 frip_min=0.01")
    args = ap.parse_args()

    overrides = parse_thresholds(args.thresholds)
    unknown = sorted(k for k in overrides if k not in DEFAULT_THRESHOLDS)
    if unknown:
        print(f"[gates_summary] warning: unknown threshold keys are ignored: {unknown}")
    thresholds = dict(DEFAULT_THRESHOLDS)
    thresholds.update({k: v for k, v in overrides.items() if k in DEFAULT_THRESHOLDS})

    rows = []
    for sample, group in parse_sample_tokens(args.samples):
        values = collect_values(sample, group, args.align_dir, args.qc_dir,
                                args.gates_dir)
        rows.append(build_row(sample, group, values, thresholds))

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="\n", encoding="utf-8") as fh:
        fh.write("\t".join(HEADER) + "\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")

    with open(args.mqc, "w", newline="\n", encoding="utf-8") as fh:
        fh.write("# id: 'gate_summary_table'\n")
        fh.write("# section_name: 'QC gates (PASS/WARN/FAIL)'\n")
        fh.write("# format: 'tsv'\n")
        fh.write("# plot_type: 'table'\n")
        fh.write("# pconfig: {'id': 'gate_summary_table', 'title': 'QC gates'}\n")
        fh.write("\t".join(HEADER) + "\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")

    counts = {}
    for row in rows:
        counts[row[-1]] = counts.get(row[-1], 0) + 1
    summary = ", ".join(f"{gate}={counts[gate]}" for gate in sorted(counts))
    print(f"[gates_summary] {len(rows)} sample(s) summarized ({summary}) -> {args.out}")


if __name__ == "__main__":
    main()
