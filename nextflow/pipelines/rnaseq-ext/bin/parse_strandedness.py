#!/usr/bin/env python3
"""Parse an RSeQC infer_experiment.py report into a featureCounts strand code.

Ports the legacy check_strandedness rule (rna-seq workflow/rules/align.smk)
verbatim: the last two fractions of the report are compared with a 0.4
margin, then mapped through the legacy quantification rule
(workflow/rules/quant.smk): secondstrand -> 1, firststrand -> 2, else 0.

Usage:
    python parse_strandedness.py <infer_experiment.out> <output.strand.txt>
"""
import re
import sys

MARGIN = 0.4


def parse_strand(report_text):
    """Return '1' (secondstrand) | '2' (firststrand) | '0' (unstranded)."""
    # infer_experiment.py writes each fraction on its label line
    # ('Fraction of reads explained by "1++,1--": 0.9053'); the legacy rule
    # compares the last two such values with a 0.4 margin.
    values = [float(f) for f in re.findall(r":\s*([0-9.]+)\s*$", report_text, re.MULTILINE)]
    if len(values) < 2:
        raise ValueError(f"could not read two fractions from infer_experiment report: {report_text!r}")
    f1, f2 = values[-2], values[-1]
    if f1 - f2 > MARGIN:
        return "1"  # secondstrand (legacy quant mapping: -s 1)
    if f2 - f1 > MARGIN:
        return "2"  # firststrand (legacy quant mapping: -s 2)
    return "0"


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    text = open(sys.argv[1]).read()
    code = parse_strand(text)
    with open(sys.argv[2], "w") as fh:
        fh.write(code + "\n")


if __name__ == "__main__":
    main()
