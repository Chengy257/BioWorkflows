#!/usr/bin/env python3
"""Resolve and validate the srna-seq software runtime.

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
    "rscript": "Rscript",
    "trim_galore": "trim_galore",
    "multiqc": "multiqc",
    "bowtie": "bowtie",
    "bowtie_build": "bowtie-build",
}

# Preflight scope: exactly the srna-seq executables (no R). The shared
# framework reads the "rscript" entry for resolution/export, but srna-seq
# has no R stage, so the preflight never demands Rscript.
_PIPELINE_TOOLS = [name for name in DEFAULT_TOOLS if name != "rscript"]

SPEC = WorkflowSpec(
    env_prefix="SRNA",
    default_tools=DEFAULT_TOOLS,
    pipeline_tools={"default": _PIPELINE_TOOLS},
    r_packages={},
    tool_defaults={"python": "python3"},
    orgdb={},
)

if __name__ == "__main__":
    raise SystemExit(run_cli(SPEC))
