#!/usr/bin/env python3
"""Resolve and validate the unified software runtime (chip/CUT&Tag/ATAC/FAIRE).

Ported from rna-seq v0.8.0's runtime_config.py:
- export subcommand: prints eval-able `export KEY=VALUE` lines for run.sh to inject;
- check subcommand: preflight; prints a human-readable list of missing tools /
  mismatched R versions / missing R packages and exits non-zero.
"""
import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

import yaml

# Default tool resolution table: keys are the logical names from WORKFLOW_TOOLS,
# values are command names (or default executables); the tools section of
# config/software.yaml can override them by the same key with a command name or
# an absolute path.
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

# The chip workflow has no pipeline split, so preflight checks this one tool set;
# scope filtering: software checks the non-R tools (rscript belongs to the r
# scope), r checks the R binary itself, its version, and packages.
WORKFLOW_TOOLS = ["python", "snakemake", "rscript", "trim_galore", "bowtie2", "bowtie2-build", "samtools", "picard", "macs2", "bedtools", "deeptools", "spp", "multiqc", "fastqc"]

# R package-level check list (spp is checked at the tool level, no package check).
R_PACKAGES = ["GenomicFeatures", "ChIPseeker"]


def load_yaml(path):
    with open(path, encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def expand(value):
    return os.path.expanduser(os.path.expandvars(str(value)))


def conda_prefix(env_cfg):
    prefix = expand(env_cfg.get("conda_prefix") or "")
    name = str(env_cfg.get("conda_name") or "").strip()
    if prefix:
        return os.path.abspath(prefix)
    if not name:
        current = os.environ.get("CONDA_PREFIX", "")
        if current:
            return current
        raise RuntimeError("environment.type=conda requires conda_prefix or conda_name (or an active CONDA_PREFIX)")
    conda = shutil.which("conda")
    if not conda:
        raise RuntimeError(f"Cannot resolve Conda environment '{name}': conda was not found in PATH")
    if name == "base":
        proc = subprocess.run([conda, "info", "--base"], capture_output=True, text=True)
        if proc.returncode:
            raise RuntimeError(proc.stderr.strip() or "conda info --base failed")
        return proc.stdout.strip()
    proc = subprocess.run([conda, "env", "list", "--json"], capture_output=True, text=True)
    if proc.returncode:
        raise RuntimeError(proc.stderr.strip() or "conda env list failed")
    envs = json.loads(proc.stdout).get("envs", [])
    matches = [p for p in envs if os.path.basename(p.rstrip(os.sep)) == name]
    if not matches:
        raise RuntimeError(f"Conda environment '{name}' was not found")
    if len(matches) > 1:
        raise RuntimeError(f"Conda environment name '{name}' is ambiguous; use conda_prefix instead")
    return matches[0]


def resolve_runtime(cfg):
    env_cfg = cfg.get("environment") or {}
    env_type = str(env_cfg.get("type") or "system").lower()
    if env_type not in {"system", "conda"}:
        raise RuntimeError("environment.type must be 'system' or 'conda'")
    prefix = ""
    path_parts = []
    if env_type == "conda":
        prefix = conda_prefix(env_cfg)
        if not os.path.isdir(prefix):
            raise RuntimeError(f"Conda prefix does not exist: {prefix}")

    tools_cfg = cfg.get("tools") or {}
    # Tool-level absolute-path overrides apply only to that tool; the override
    # directory is never prepended wholesale to PATH, which would globally
    # shadow unrelated commands from the main environment.
    if prefix:
        path_parts.append(os.path.join(prefix, "bin"))
    path_parts.append(os.environ.get("PATH", ""))
    runtime_path = os.pathsep.join(x for x in path_parts if x)

    resolved = {}
    for name, default in DEFAULT_TOOLS.items():
        requested = expand(tools_cfg.get(name) or default)
        if os.path.isabs(requested):
            resolved[name] = requested
        else:
            resolved[name] = shutil.which(requested, path=runtime_path) or requested

    r_cfg = cfg.get("r") or {}
    r_requested = expand(r_cfg.get("rscript") or resolved["rscript"])
    if os.path.isabs(r_requested):
        rscript = r_requested
    else:
        rscript = shutil.which(r_requested, path=runtime_path) or r_requested
    resolved["rscript"] = rscript

    lib_paths = [expand(p) for p in (r_cfg.get("lib_paths") or []) if str(p).strip()]
    mode = str(r_cfg.get("lib_mode") or "prepend").lower()
    if mode not in {"prepend", "append", "replace"}:
        raise RuntimeError("r.lib_mode must be prepend, append, or replace")
    existing = [p for p in os.environ.get("R_LIBS_USER", "").split(os.pathsep) if p]
    if mode == "prepend":
        libs = lib_paths + existing
    elif mode == "append":
        libs = existing + lib_paths
    else:
        libs = lib_paths
    libs = list(dict.fromkeys(libs))

    return {
        "type": env_type,
        "prefix": prefix,
        "path": runtime_path,
        "tools": resolved,
        "r_libs": libs,
        "r_cfg": r_cfg,
        "strict": bool(env_cfg.get("strict", True)),
    }


def shell_exports(rt, config_path):
    values = {
        "CHIP_SOFTWARE_CONFIG": os.path.abspath(config_path),
        "CHIP_SOFTWARE_TYPE": rt["type"],
        "CHIP_ENV_PREFIX": rt["prefix"],
        "PATH": rt["path"],
        "CHIP_RSCRIPT": rt["tools"]["rscript"],
        "CHIP_PYTHON": rt["tools"]["python"],
    }
    if rt["r_libs"]:
        values["R_LIBS_USER"] = os.pathsep.join(rt["r_libs"])
    for name, value in rt["tools"].items():
        key = "CHIP_TOOL_" + name.upper().replace("-", "_")
        values[key] = value
    for key, value in values.items():
        print(f"export {key}={shlex.quote(str(value))}")


def cmd_exists(command):
    if os.path.isabs(command):
        return os.path.isfile(command) and os.access(command, os.X_OK)
    return shutil.which(command) is not None


def r_eval(rscript, expression):
    return subprocess.run([rscript, "--vanilla", "-e", expression], capture_output=True, text=True)


def check_runtime(rt, scope, analysis_cfg=None):
    # analysis_cfg is a reserved parameter (--analysis-config extension point); the
    # chip workflow centralizes sample-table and config validation at Snakefile
    # parse time, so no extra checks beyond the package check happen here.
    errors = []
    warnings = []
    print(f"[runtime] Environment: {rt['type']}" + (f" ({rt['prefix']})" if rt["prefix"] else ""))

    if scope in {"all", "software"}:
        for name in WORKFLOW_TOOLS:
            if name == "rscript":
                continue  # rscript belongs to the r scope check
            exe = rt["tools"][name]
            if cmd_exists(exe):
                print(f"[runtime] [OK] {name}: {exe}")
            else:
                errors.append(f"Required executable not found: {name} ({exe})")

    if scope in {"all", "r"}:
        rscript = rt["tools"]["rscript"]
        if not cmd_exists(rscript):
            errors.append(f"Rscript not found: {rscript}")
        else:
            ver = r_eval(rscript, "cat(as.character(getRversion()))")
            if ver.returncode:
                errors.append(f"Unable to run Rscript: {ver.stderr.strip()}")
            else:
                detected = ver.stdout.strip()
                expected = str(rt["r_cfg"].get("version") or "").strip()
                policy = str(rt["r_cfg"].get("version_check") or "major_minor").lower()
                print(f"[runtime] [OK] R: {detected} via {rscript}")
                if expected and policy != "off":
                    if policy in {"major_minor", "warn"}:
                        match = detected.split(".")[:2] == expected.split(".")[:2]
                    elif policy == "exact":
                        match = detected == expected
                    else:
                        errors.append(f"Unknown r.version_check policy: {policy}")
                        match = True
                    if not match:
                        msg = f"R version mismatch: expected {expected}, detected {detected}"
                        (warnings if policy == "warn" else errors).append(msg)
            for path in rt["r_libs"]:
                if os.path.isdir(path):
                    print(f"[runtime] [OK] R library: {path}")
                else:
                    (errors if rt["strict"] else warnings).append(f"R library path does not exist: {path}")
            pkgs = list(R_PACKAGES)
            quoted = ",".join(json.dumps(p) for p in pkgs)
            expr = f"p<-c({quoted});m<-p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)];cat(paste(m,collapse='\\n'))"
            pkgcheck = r_eval(rscript, expr)
            if pkgcheck.returncode:
                errors.append(f"R package check failed: {pkgcheck.stderr.strip()}")
            else:
                missing = [x for x in pkgcheck.stdout.splitlines() if x.strip()]
                if missing:
                    errors.append("Missing R package(s): " + ", ".join(missing))
                    sources = rt["r_cfg"].get("package_sources") or {}
                    target_lib = rt["r_libs"][0] if rt["r_libs"] else "<R-library>"
                    r_bin = os.path.join(os.path.dirname(rscript), "R") if os.path.isabs(rscript) else "R"
                    for package in missing:
                        source = expand(sources.get(package) or "")
                        if source:
                            warnings.append(f"Install hint for {package}: {r_bin} CMD INSTALL -l {target_lib} {source}")
                else:
                    print(f"[runtime] [OK] R packages: {len(pkgs)} required package(s)")

    for msg in warnings:
        print(f"[runtime] [WARN] {msg}", file=sys.stderr)
    for msg in errors:
        print(f"[runtime] [ERROR] {msg}", file=sys.stderr)
    # With strict=false the preflight only warns and never blocks (warn-only);
    # with true, any error exits non-zero.
    return 0 if not errors or not rt["strict"] else 1


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    p_export = sub.add_parser("export")
    p_export.add_argument("--config", required=True)
    p_check = sub.add_parser("check")
    p_check.add_argument("--config", required=True)
    p_check.add_argument("--scope", choices=["all", "software", "r"], default="all")
    p_check.add_argument("--analysis-config")  # reserved extension point, disabled by default
    args = parser.parse_args()
    cfg = load_yaml(args.config)
    try:
        rt = resolve_runtime(cfg)
    except RuntimeError as exc:
        print(f"[runtime] [ERROR] {exc}", file=sys.stderr)
        return 2
    if args.command == "export":
        shell_exports(rt, args.config)
        return 0
    analysis = load_yaml(args.analysis_config) if args.analysis_config else {}
    return check_runtime(rt, args.scope, analysis)


if __name__ == "__main__":
    raise SystemExit(main())
