#!/usr/bin/env python3
"""Capture the resolved runtime used by the RNA-seq workflow.

Thin wrapper over the shared BioWorkflows version collector
(shared/python/bioworkflows_versions.py).
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_SHARED = os.environ.get("BIO_WORKFLOWS_SHARED") or os.path.abspath(
    os.path.join(_HERE, "..", "..", "..", "shared"))
sys.path.insert(0, os.path.join(_SHARED, "python"))

from bioworkflows_versions import run_cli  # noqa: E402

TOOLS = [
    "python", "rscript", "fastqc", "trim_galore", "multiqc", "star", "samtools",
    "infer_experiment", "stringtie", "gffcompare", "gffread", "bioawk", "diamond",
    "hmmscan", "bedtools", "pfam_scan", "cpc2", "cnci_python",
]

TOOL_DEFAULTS = {
    "python": "python3", "rscript": "Rscript", "star": "STAR",
    "infer_experiment": "infer_experiment.py",
}

DATABASES = {
    "pfam": "DB_PFAM",
    "nr_diamond": "DB_NR_DIAMOND",
    "cnci_dir": "CNCI_DIR",
}

if __name__ == "__main__":
    run_cli("RNASEQ", TOOLS, TOOL_DEFAULTS, DATABASES)
