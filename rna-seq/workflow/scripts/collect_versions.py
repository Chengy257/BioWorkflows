#!/usr/bin/env python3
"""Capture the resolved runtime used by the RNA-seq workflow."""
import argparse
import datetime
import os
from pathlib import Path
import shutil
import subprocess

import yaml

TOOLS = [
    "python", "rscript", "fastqc", "trim_galore", "multiqc", "star", "samtools",
    "infer_experiment", "stringtie", "gffcompare", "gffread", "bioawk", "diamond",
    "hmmscan", "bedtools", "pfam_scan", "cpc2", "cnci_python",
]


def command_output(args, env=None):
    try:
        proc = subprocess.run(args, capture_output=True, text=True, timeout=20, env=env)
        return (proc.stdout or proc.stderr).strip()
    except Exception:
        return ""


def git_commit(path):
    return command_output(["git", "-C", path, "rev-parse", "--short", "HEAD"])


def resolved_tool(name):
    key = "RNASEQ_TOOL_" + name.upper().replace("-", "_")
    value = os.environ.get(key, "")
    if value:
        return value
    defaults = {"python": "python3", "rscript": "Rscript", "star": "STAR", "infer_experiment": "infer_experiment.py"}
    cmd = defaults.get(name, name)
    return shutil.which(cmd) or cmd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--software-config", required=True)
    ap.add_argument("--workflow", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--snakemake-version", default="")
    args = ap.parse_args()

    data = {
        "generated_at": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "snakemake": args.snakemake_version or command_output(["snakemake", "--version"]) or "unknown",
        "git_commit": git_commit(args.workflow) or "unknown",
        "software_config": os.path.abspath(args.software_config),
        "runtime": {
            "type": os.environ.get("RNASEQ_SOFTWARE_TYPE", "system"),
            "conda_prefix": os.environ.get("RNASEQ_ENV_PREFIX", ""),
            "rscript": os.environ.get("RNASEQ_RSCRIPT", resolved_tool("rscript")),
            "r_libs_user": os.environ.get("R_LIBS_USER", ""),
        },
        "tools": {name: resolved_tool(name) for name in TOOLS},
        "databases": {
            "pfam": os.environ.get("RNASEQ_DB_PFAM", ""),
            "nr_diamond": os.environ.get("RNASEQ_DB_NR_DIAMOND", ""),
            "cnci_dir": os.environ.get("RNASEQ_CNCI_DIR", ""),
        },
    }
    rscript = data["runtime"]["rscript"]
    if shutil.which(rscript) or os.path.isfile(rscript):
        rv = command_output([rscript, "--vanilla", "-e", "cat(as.character(getRversion()))"])
        if rv:
            data["runtime"]["r_version"] = rv
        libs = command_output([rscript, "--vanilla", "-e", "cat(paste(.libPaths(), collapse=':'))"])
        if libs:
            data["runtime"]["r_lib_paths"] = libs.split(":")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as handle:
        yaml.safe_dump(data, handle, sort_keys=False, allow_unicode=True)
    print(f"[collect_versions] -> {args.out}")


if __name__ == "__main__":
    main()
