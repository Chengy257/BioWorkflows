#!/usr/bin/env python3
"""Capture the resolved runtime used by the chip/CUT&Tag/ATAC/FAIRE workflow.

Thin wrapper over the shared BioWorkflows version collector
(shared/python/bioworkflows_versions.py). Reads the CHIP_* environment
injected by run.sh (unified environment + software.yaml resolution) and
writes a YAML record into the results directory for methods sections and
reproducibility audits.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_SHARED = os.environ.get("BIO_WORKFLOWS_SHARED") or os.path.abspath(
    os.path.join(_HERE, "..", "..", "..", "shared"))
sys.path.insert(0, os.path.join(_SHARED, "python"))

from bioworkflows_versions import run_cli  # noqa: E402

# Keep in sync with WORKFLOW_TOOLS in runtime_config.py: logical tool names.
TOOLS = [
    "python", "snakemake", "rscript", "trim_galore", "bowtie2", "bowtie2-build",
    "samtools", "picard", "macs2", "bedtools", "deeptools", "spp", "multiqc", "fastqc",
]

TOOL_DEFAULTS = {
    "python": "python3", "rscript": "Rscript", "deeptools": "bamCoverage",
    "spp": "run_spp.R", "snakemake": "snakemake",
}

DATABASES = {}

if __name__ == "__main__":
    run_cli("CHIP", TOOLS, TOOL_DEFAULTS, DATABASES)
