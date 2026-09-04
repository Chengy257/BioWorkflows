#!/usr/bin/env python3
"""Resolve and validate the unified RNA-seq software runtime.

Thin wrapper over the shared BioWorkflows runtime framework
(shared/python/bioworkflows_runtime.py); the tool tables and lncRNA-specific
checks live here. CLI:
  export --config software.yaml         # emit `export KEY=VALUE` lines
  check  --config software.yaml --pipeline upstream|deg|as|lncrna [--scope ...]
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
    "fastqc": "fastqc",
    "trim_galore": "trim_galore",
    "multiqc": "multiqc",
    "star": "STAR",
    "samtools": "samtools",
    "infer_experiment": "infer_experiment.py",
    "stringtie": "stringtie",
    "gffcompare": "gffcompare",
    "gffread": "gffread",
    "bioawk": "bioawk",
    "diamond": "diamond",
    "hmmscan": "hmmscan",
    "bedtools": "bedtools",
    "pfam_scan": "pfam_scan.pl",
    "cpc2": "CPC2.py",
    "cnci_python": "python2",
}

PIPELINE_TOOLS = {
    "upstream": ["fastqc", "trim_galore", "multiqc", "star", "samtools", "infer_experiment", "rscript", "python"],
    "deg": ["fastqc", "trim_galore", "multiqc", "star", "samtools", "infer_experiment", "rscript", "python"],
    "as": ["fastqc", "trim_galore", "multiqc", "star", "samtools", "infer_experiment", "stringtie", "gffcompare", "gffread", "python"],
    "lncrna": ["fastqc", "trim_galore", "multiqc", "star", "samtools", "infer_experiment", "rscript", "python", "stringtie", "gffcompare", "gffread", "bioawk", "diamond", "hmmscan", "bedtools", "pfam_scan", "cpc2", "cnci_python"],
}

R_PACKAGES = {
    "upstream": ["argparser", "Rsubread", "limma", "edgeR", "getopt"],
    "deg": ["argparser", "Rsubread", "limma", "edgeR", "getopt", "DESeq2", "ggplot2", "BiocParallel", "gplots", "RColorBrewer", "amap", "clusterProfiler", "enrichplot", "aPEAR", "svglite", "magrittr", "dplyr", "VennDiagram", "UpSetR", "R.utils"],
    "as": [],
    "lncrna": ["argparser", "Rsubread", "limma", "edgeR", "getopt"],
}


def _extra_exports(rt):
    return {
        "RNASEQ_CNCI_DIR": rt["paths"].get("cnci_dir", ""),
        "RNASEQ_DB_PFAM": rt["databases"].get("pfam", ""),
        "RNASEQ_DB_NR_DIAMOND": rt["databases"].get("nr_diamond", ""),
    }


def _extra_checks(rt, pipeline, errors, warnings):
    if pipeline != "lncrna":
        return
    cnci_dir = rt["paths"].get("cnci_dir", "")
    if not cnci_dir or not os.path.isfile(os.path.join(cnci_dir, "CNCI.py")):
        errors.append("lncRNA requires paths.cnci_dir containing CNCI.py")
    else:
        print(f"[runtime] [OK] cnci_dir: {cnci_dir}")
    for label, path in (("Pfam database", rt["databases"].get("pfam", "")),
                        ("NR DIAMOND database", rt["databases"].get("nr_diamond", ""))):
        if not path or not os.path.exists(path):
            errors.append(f"lncRNA requires {label}: {path or '<not configured>'}")
        else:
            print(f"[runtime] [OK] {label}: {path}")


SPEC = WorkflowSpec(
    env_prefix="RNASEQ",
    default_tools=DEFAULT_TOOLS,
    pipeline_tools=PIPELINE_TOOLS,
    r_packages=R_PACKAGES,
    tool_defaults={"python": "python3", "rscript": "Rscript", "star": "STAR",
                   "infer_experiment": "infer_experiment.py"},
    orgdb={"deg": {"hsa": "org.Hs.eg.db", "osa": "org.Osativa.eg.db"}},
    extra_exports=_extra_exports,
    extra_checks=_extra_checks,
)

if __name__ == "__main__":
    raise SystemExit(run_cli(SPEC))
