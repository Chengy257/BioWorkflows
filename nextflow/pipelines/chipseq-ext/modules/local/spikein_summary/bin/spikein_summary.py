#!/usr/bin/env python3
"""Spike-in normalization summary for combined-reference runs (host and
spike-in genomes aligned together), one row per sample from per-sample
samtools idxstats tables.

Per sample: spike_mapped (idxstats mapped reads on contigs matching the
configured patterns), host_mapped (all other mapped reads),
spike_rate = spike_mapped / (host_mapped + spike_mapped), and scale_factor =
1e6 / max(spike_mapped, 1) - the factor signal tracks are multiplied by when
spike-in rescaling is enabled.

Contig classification: patterns of <= 3 chars must equal the contig name
case-insensitively; longer patterns match as substrings (same convention as
organelle_summary.py). A sample whose idxstats is unparseable renders NA.

Usage:
    spikein_summary.py --patterns spike --name lambda --out Spikein_summary.tsv \
        --mqc Spikein_summary_mqc.tsv sample1.idxstats [sample2.idxstats ...]
"""
import argparse
import os

HEADER = ["sample", "spike_mapped", "host_mapped", "spike_rate", "scale_factor"]


def is_spike_contig(name, patterns):
    lowered = name.lower()
    for pattern in patterns:
        p = str(pattern).lower()
        if len(p) <= 3:
            if lowered == p:
                return True
        elif p in lowered:
            return True
    return False


def parse_idxstats(path, patterns):
    """(spike_mapped, host_mapped) from an idxstats table; (None, None) when
    the file is missing or carries no parseable rows."""
    spike = None
    host = 0
    try:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 4 or cols[0] == "*":
                    continue
                mapped = int(cols[2])
                if is_spike_contig(cols[0], patterns):
                    spike = (spike or 0) + mapped
                else:
                    host += mapped
    except OSError:
        return None, None
    return spike, host


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--patterns", required=True,
                    help="comma-separated contig-name patterns for spike contigs")
    ap.add_argument("--name", required=True,
                    help="spike-in label used in the MultiQC section title")
    ap.add_argument("--out", required=True, help="summary TSV output path")
    ap.add_argument("--mqc", required=True, help="MultiQC custom-content TSV output path")
    ap.add_argument("idxstats", nargs="+", help="per-sample idxstats TSV files")
    args = ap.parse_args()

    patterns = [p for p in (x.strip() for x in args.patterns.split(",")) if p]

    rows = []
    for path in args.idxstats:
        sample = os.path.basename(path)
        for suffix in ("_idxstats.tsv", ".idxstats"):
            if sample.endswith(suffix):
                sample = sample[:-len(suffix)]
                break
        spike, host = parse_idxstats(path, patterns)
        if spike is None:
            rows.append([sample, "NA", "NA", "NA", "NA"])
            continue
        total = spike + host
        rate = (spike / total) if total else 0.0
        factor = 1e6 / max(spike, 1)
        rows.append([sample, spike, host, f"{rate:.4f}", f"{factor:.4f}"])

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


if __name__ == "__main__":
    raise SystemExit(main())
