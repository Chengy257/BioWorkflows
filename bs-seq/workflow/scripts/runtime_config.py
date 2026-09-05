#!/usr/bin/env python3
"""Resolve and validate the bs-seq software runtime.

Thin wrapper over the shared BioWorkflows runtime framework
(shared/python/bioworkflows_runtime.py). CLI:
  export --config software.yaml         # emit `export KEY=VALUE` lines
  check  --config software.yaml [--scope all|tools]
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_SHARED = os.environ.get("BIO_WORKFLOWS_SHARED") or os.path.abspath(
    os.path.join(_HERE, "..", "..", "..", "shared"))
sys.path.insert(0, os.path.join(_SHARED, "python"))

from bioworkflows_runtime import WorkflowSpec, run_cli  # noqa: E402

DEFAULT_TOOLS = {
    "python": "python3",
    "trim_galore": "trim_galore",
    "multiqc": "multiqc",
    "bismark": "bismark",
    "bismark_genome_preparation": "bismark_genome_preparation",
    "deduplicate_bismark": "deduplicate_bismark",
    "bismark_methylation_extractor": "bismark_methylation_extractor",
    "coverage2cytosine": "coverage2cytosine",
    "bam2nuc": "bam2nuc",
    "bismark2report": "bismark2report",
    "bismark2summary": "bismark2summary",
    "bowtie2": "bowtie2",
    "samtools": "samtools",
    "fastqc": "fastqc",
    # Required by the shared runtime framework (resolve_runtime reads the
    # "rscript" entry and the export emits BSSEQ_RSCRIPT); bs-seq itself has
    # no R stage, so it is kept OUT of pipeline_tools below: the preflight
    # never demands Rscript, and with r_packages empty it is only resolved,
    # never existence-checked.
    "rscript": "Rscript",
}

# Preflight scope: exactly the bs-seq executables (no R).
_PIPELINE_TOOLS = [name for name in DEFAULT_TOOLS if name != "rscript"]

SPEC = WorkflowSpec(
    env_prefix="BSSEQ",
    default_tools=DEFAULT_TOOLS,
    pipeline_tools={"default": _PIPELINE_TOOLS},
    r_packages={},
    tool_defaults={"python": "python3"},
    orgdb={},
)

if __name__ == "__main__":
    raise SystemExit(run_cli(SPEC))
