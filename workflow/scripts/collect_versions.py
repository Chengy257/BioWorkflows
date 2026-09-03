#!/usr/bin/env python3
"""记录 chip/CUT&Tag/ATAC/FAIRE 流程实际使用的运行时与工具版本。

移植自 rna-seq v0.8.0 的 collect_versions.py：
- 读 run.sh 注入的 CHIP_SOFTWARE_CONFIG / CHIP_TOOL_* 环境变量（统一环境 + software.yaml 路线）；
- 逐工具解析实际可执行文件路径，snakemake 版本在脚本内经 importlib.metadata 自取；
- 输出 YAML 到 results 的 5.QC/，供论文方法节与复现审计使用。
"""
import argparse
import datetime
import os
from pathlib import Path
import shutil
import subprocess

import yaml

# 与 runtime_config.py 的 WORKFLOW_TOOLS 保持一致：14 个流程工具的逻辑名
TOOLS = [
    "python", "snakemake", "rscript", "trim_galore", "bowtie2", "bowtie2-build",
    "samtools", "picard", "macs2", "bedtools", "deeptools", "spp", "multiqc", "fastqc",
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
    # run.sh 注入的 CHIP_TOOL_<NAME>（绝对路径）优先；否则按默认命令名在 PATH 解析
    key = "CHIP_TOOL_" + name.upper().replace("-", "_")
    value = os.environ.get(key, "")
    if value:
        return value
    defaults = {"python": "python3", "rscript": "Rscript", "spp": "spp.R"}
    cmd = defaults.get(name, name)
    return shutil.which(cmd) or cmd


def snakemake_version():
    # 直接查当前解释器安装的 snakemake 发行元数据；未安装（如裸机调试）记 unknown
    try:
        from importlib.metadata import version
        return version("snakemake")
    except Exception:
        return "unknown"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--software-config", required=True)
    ap.add_argument("--workflow", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    data = {
        "generated_at": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "snakemake": snakemake_version(),
        "git_commit": git_commit(args.workflow) or "unknown",
        "software_config": os.path.abspath(args.software_config),
        "runtime": {
            "type": os.environ.get("CHIP_SOFTWARE_TYPE", "system"),
            "conda_prefix": os.environ.get("CHIP_ENV_PREFIX", ""),
            "rscript": os.environ.get("CHIP_RSCRIPT", resolved_tool("rscript")),
            "r_libs_user": os.environ.get("R_LIBS_USER", ""),
        },
        "tools": {name: resolved_tool(name) for name in TOOLS},
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
