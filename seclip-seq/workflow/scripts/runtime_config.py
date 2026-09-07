#!/usr/bin/env python3
"""Resolve and validate the seclip-seq software runtime.

Thin wrapper over the shared BioWorkflows runtime framework
(shared/python/bioworkflows_runtime.py); the tool table and CLIPper
handling live here. CLI:
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
    "umi_tools": "umi_tools",
    "cutadapt": "cutadapt",
    "seqkit": "seqkit",
    "bgzip": "bgzip",
    "fastqc": "fastqc",
    "multiqc": "multiqc",
    "star": "STAR",
    "samtools": "samtools",
    "pureclip": "pureclip",
}


def _extra_exports(rt):
    return {"SECLIP_TOOL_CLIPPER": rt["paths"].get("clipper", "")}


def _extra_checks(rt, pipeline, errors, warnings):
    # CLIPper is an external legacy install; absence is not an error (the
    # workflow auto-skips it) but is always reported.
    clipper = rt["paths"].get("clipper", "")
    if clipper:
        print(f"[runtime] [OK] clipper: {clipper}")
    else:
        warnings.append("CLIPper not configured (software.yaml paths.clipper); "
                        "CLIPper peak calling will be skipped")


SPEC = WorkflowSpec(
    env_prefix="SECLIP",
    default_tools=DEFAULT_TOOLS,
    pipeline_tools={"default": list(DEFAULT_TOOLS)},
    r_packages={},
    tool_defaults={"python": "python3", "star": "STAR"},
    orgdb={},
    extra_exports=_extra_exports,
    extra_checks=_extra_checks,
)

if __name__ == "__main__":
    raise SystemExit(run_cli(SPEC))
