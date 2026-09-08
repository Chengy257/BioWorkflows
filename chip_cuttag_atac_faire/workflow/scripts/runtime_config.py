#!/usr/bin/env python3
"""Resolve and validate the unified software runtime (chip/CUT&Tag/ATAC/FAIRE).

Thin wrapper over the shared BioWorkflows runtime framework
(shared/python/bioworkflows_runtime.py). CLI:
  export --config software.yaml                     # emit `export KEY=VALUE` lines
  check  --config software.yaml [--scope all|software|r] [--analysis-config ...]
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_SHARED = os.environ.get("BIO_WORKFLOWS_SHARED") or os.path.abspath(
    os.path.join(_HERE, "..", "..", "..", "shared"))
sys.path.insert(0, os.path.join(_SHARED, "python"))

from bioworkflows_runtime import WorkflowSpec, run_cli  # noqa: E402

# Default tool resolution table: keys are the logical tool names, values are
# command names (or default executables); the tools section of
# config/software.yaml can override them by the same key with a command name or
# an absolute path. The insertion order is load-bearing: it fixes the order of
# the emitted CHIP_TOOL_* export lines, so keep it stable.
DEFAULT_TOOLS = {
    "python": "python3",
    "rscript": "Rscript",
    "snakemake": "snakemake",
    "fastqc": "fastqc",
    "trim_galore": "trim_galore",
    "multiqc": "multiqc",
    "bowtie2": "bowtie2",
    "bowtie2-build": "bowtie2-build",
    "samtools": "samtools",
    "picard": "picard",
    "macs2": "macs2",
    "bedtools": "bedtools",
    "deeptools": "bamCoverage",
    "spp": "run_spp.R",
}

# Preflight scope for --scope software: the non-R tools. "rscript" stays in
# DEFAULT_TOOLS (the framework resolves it unconditionally and export emits
# CHIP_RSCRIPT / CHIP_TOOL_RSCRIPT) but is excluded here so the software
# preflight never demands R; the r scope covers Rscript itself, its version,
# and the R packages below.
_PIPELINE_TOOLS = [name for name in DEFAULT_TOOLS if name != "rscript"]

# R package-level check list (spp is checked at the tool level, no package
# check). Kept non-empty on purpose: it makes the shared framework's r scope
# always active for --scope r / --scope all, matching the previous standalone
# resolver.
R_PACKAGES = {"default": ["GenomicFeatures", "ChIPseeker"]}


def _extra_exports(rt):
    """Optional external tools resolved from software.yaml paths: (not part
    of DEFAULT_TOOLS, so the preflight never demands them). `idr` is the
    classic python2 IDR tool needed only when peak.replicate.enabled is true;
    it installs separately, e.g. `conda create -n idr -c bioconda idr=2.0.4`,
    and reaches the rules as CHIP_IDR (common.smk falls back to a bare
    `idr` on PATH)."""
    values = {}
    idr = (rt["paths"] or {}).get("idr")
    if idr:
        values["CHIP_IDR"] = idr
    return values


SPEC = WorkflowSpec(
    env_prefix="CHIP",
    default_tools=DEFAULT_TOOLS,
    pipeline_tools={"default": _PIPELINE_TOOLS},
    r_packages=R_PACKAGES,
    tool_defaults={"python": "python3"},
    # No OrgDb: peak annotation is GTF-based (makeTxDbFromGFF in
    # annoPeak_batch.R), so there is nothing to add here.
    orgdb={},
    extra_exports=_extra_exports,
)

if __name__ == "__main__":
    raise SystemExit(run_cli(SPEC))
