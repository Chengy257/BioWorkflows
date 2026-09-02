#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""收集工作流版本信息，生成 software_versions.yaml（对应审查 P2-9）。

解析 conda 环境 yaml（不依赖 pyyaml），并记录 snakemake 版本与 git 提交号。
用法:
    python collect_versions.py --envs <envs 目录> --workflow <workflow 目录> --out <输出文件>
"""
import argparse
import datetime
import os
import re
import subprocess


def sh(cmd, cwd=None):
    try:
        r = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True,
                           text=True, errors="replace", timeout=60)
        return r.stdout.strip()
    except Exception:
        return ""


def parse_env(path):
    """解析 env yaml 的 dependencies 段：返回 [(包名, 版本|'')]"""
    pkgs = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.lstrip().startswith("#"):
                continue
            m = re.match(r"\s*-\s*([A-Za-z0-9._-]+?)\s*=\s*([A-Za-z0-9._+-]+)\s*$", line)
            if m:
                pkgs.append((m.group(1), m.group(2)))
    return pkgs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--envs", required=True, help="conda 环境定义目录")
    ap.add_argument("--workflow", required=True, help="workflow 目录（用于读取 git 提交）")
    ap.add_argument("--out", required=True, help="输出文件")
    args = ap.parse_args()

    lines = [
        "# software_versions.yaml — 由流程自动生成（记录本次运行所用的环境定义）",
        f"date: '{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}'",
        f"snakemake: '{sh('snakemake --version') or 'unknown'}'",
    ]
    commit = sh("git rev-parse --short HEAD", cwd=args.workflow)
    if commit:
        lines.append(f"git_commit: '{commit}'")

    for env in sorted(os.listdir(args.envs)):
        if not env.endswith(".yaml"):
            continue
        lines.append(f"{env[:-5]}:")
        for pkg, ver in parse_env(os.path.join(args.envs, env)):
            lines.append(f"  {pkg}: '{ver or 'unpinned'}'")

    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"[collect_versions] -> {args.out}")


if __name__ == "__main__":
    main()
