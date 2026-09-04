#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""样本表（sample_info.csv）校验器。

用法:
    python validate_samples.py <sample_info.csv> [control_group]

表头必须含 id,group；可选 layout（PE/SE/auto）与 batch 列。
对照组名由 config 的 control_group 指定（默认 control）。

检查项:
    [error] 缺少 id/group 列；id 重复/为空/含 '-'、空格、'/'；组名为空
    [error] 对照组不在 group 列中
    [warn ] 样本 id 互为前缀；某组样本数 < 2；layout 取值非法
有 error 时退出码为 1。
"""
import csv
import sys
from collections import Counter


def report(errors, warnings):
    for w in warnings:
        print(f"[WARN ] {w}", file=sys.stderr)
    for e in errors:
        print(f"[ERROR] {e}", file=sys.stderr)
    print(f"[validate_samples] {len(errors)} errors, {len(warnings)} warnings", file=sys.stderr)


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    path = sys.argv[1]
    control = sys.argv[2] if len(sys.argv) > 2 else "control"

    errors, warnings = [], []
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        fields = reader.fieldnames or []
        rows = list(reader)

    for col in ("id", "group"):
        if col not in fields:
            errors.append(f"缺少必需列: {col}（当前表头: {fields}）")
    if errors:
        report(errors, warnings)
        sys.exit(1)

    ids, groups = [], set()
    for i, row in enumerate(rows, start=2):  # 数据自第 2 行起
        sid = (row.get("id") or "").strip()
        grp = (row.get("group") or "").strip()
        if not sid:
            errors.append(f"第 {i} 行: id 为空")
            continue
        if sid in ids:
            errors.append(f"第 {i} 行: 样本 id 重复: {sid}")
        ids.append(sid)
        if "-" in sid:
            errors.append(f"第 {i} 行: id 含 '-'（各脚本会对 '-' 改名，易造成错配）: {sid}")
        if any(c.isspace() or c == "/" for c in sid):
            errors.append(f"第 {i} 行: id 含空格或 '/': {sid}")
        if not grp:
            errors.append(f"第 {i} 行: 组名为空（id={sid}）")
        groups.add(grp)
        layout = (row.get("layout") or "").strip().lower()
        if layout and layout not in ("pe", "se", "auto"):
            warnings.append(f"第 {i} 行: layout 应为 PE/SE/auto，当前: {row.get('layout')}")

    if control not in groups:
        errors.append(f"对照组 '{control}' 不在 group 列中（现有组: {sorted(groups)}）")

    for a in ids:
        for b in ids:
            if a != b and b.startswith(a):
                warnings.append(f"样本 id 互为前缀（易误配，建议改名）: {a} vs {b}")

    cnt = Counter((r.get("group") or "").strip() for r in rows)
    for g, n in sorted(cnt.items()):
        if g and n < 2:
            warnings.append(f"组 {g} 仅 {n} 个样本（无重复时 DESeq2 差异分析不可靠）")

    report(errors, warnings)
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()
