#!/usr/bin/env python3
"""Capture the resolved runtime used by the seclip-seq workflow.

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
    "python", "umi_tools", "cutadapt", "seqkit", "bgzip", "fastqc",
    "multiqc", "star", "samtools", "pureclip",
]

TOOL_DEFAULTS = {
    "python": "python3",
    "star": "STAR",
    "seqkit": "seqkit",
    "bgzip": "bgzip",
}

DATABASES = {}

if __name__ == "__main__":
    run_cli("SECLIP", TOOLS, TOOL_DEFAULTS, DATABASES)
