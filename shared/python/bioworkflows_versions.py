"""Shared software-version capture for the BioWorkflows projects.

Records the actually-resolved runtime (environment type, Conda prefix,
Rscript/R version/R libraries, per-tool executables, external databases)
plus the workflow git commit into a YAML file inside the results directory.
Each project's collect_versions.py declares its tool list and env prefix.
"""

import argparse
import datetime
import os
from pathlib import Path
import shutil
import subprocess

import yaml


def command_output(args, env=None):
    try:
        proc = subprocess.run(args, capture_output=True, text=True, timeout=20, env=env)
        return (proc.stdout or proc.stderr).strip()
    except Exception:
        return ""


def git_commit(path):
    return command_output(["git", "-C", path, "rev-parse", "--short", "HEAD"])


def resolved_tool(env_prefix, tool_defaults, name):
    key = env_prefix + "_TOOL_" + name.upper().replace("-", "_")
    value = os.environ.get(key, "")
    if value:
        return value
    cmd = tool_defaults.get(name, name)
    return shutil.which(cmd) or cmd


def collect(env_prefix, tools, tool_defaults, databases, workflow_dir,
            software_config, snakemake_version, out):
    data = {
        "generated_at": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "snakemake": snakemake_version or command_output(["snakemake", "--version"]) or "unknown",
        "git_commit": git_commit(workflow_dir) or "unknown",
        "software_config": os.path.abspath(software_config),
        "runtime": {
            "type": os.environ.get(f"{env_prefix}_SOFTWARE_TYPE", "system"),
            "conda_prefix": os.environ.get(f"{env_prefix}_ENV_PREFIX", ""),
            "rscript": os.environ.get(f"{env_prefix}_RSCRIPT",
                                      resolved_tool(env_prefix, tool_defaults, "rscript")),
            "r_libs_user": os.environ.get("R_LIBS_USER", ""),
        },
        "tools": {name: resolved_tool(env_prefix, tool_defaults, name) for name in tools},
        "databases": {
            db: os.environ.get(f"{env_prefix}_{env_suffix}", "")
            for db, env_suffix in databases.items()
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

    Path(out).parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as handle:
        yaml.safe_dump(data, handle, sort_keys=False, allow_unicode=True)
    print(f"[collect_versions] -> {out}")


def run_cli(env_prefix, tools, tool_defaults, databases):
    ap = argparse.ArgumentParser()
    ap.add_argument("--software-config", required=True)
    ap.add_argument("--workflow", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--snakemake-version", default="")
    args = ap.parse_args()
    collect(env_prefix, tools, tool_defaults, databases,
            args.workflow, args.software_config, args.snakemake_version, args.out)
