#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""按 sample_info.csv 生成 DEG 组间比较任务文件。

用法:
    python getGroups.py <sample_info.csv> <DEG_dir> [control_group]

在 <DEG_dir> 中寻找 <处理组>_vs_<对照组>_{UP,DOWN,ALL}.DEGs.txt（由 enrich.sh 生成），
为所有处理组两两组合写出任务文件：
    6.DEGcompare/combination/<A>_<r1>_<B>_<r2>
内容两行（路径 集合名），供 multi_enrich_GOKEGG.R 消费。
分组列按表头名 "group" 取（不再取最后一列），对照组名可配置。
"""
import csv
import os
import sys


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    infile, deg_dir = sys.argv[1], sys.argv[2]
    control = sys.argv[3] if len(sys.argv) > 3 else "control"

    groups = []
    with open(infile, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            g = (row.get("group") or "").strip()
            if g and g not in groups and g != control:
                groups.append(g)

    os.makedirs("6.DEGcompare/combination", exist_ok=True)
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            for r1 in ("UP", "DOWN", "ALL"):
                for r2 in ("UP", "DOWN", "ALL"):
                    out = os.path.join("6.DEGcompare", "combination",
                                       f"{groups[i]}_{r1}_{groups[j]}_{r2}")
                    with open(out, "w") as f:
                        a = os.path.join(deg_dir, f"{groups[i]}_vs_{control}_{r1}.DEGs.txt")
                        b = os.path.join(deg_dir, f"{groups[j]}_vs_{control}_{r2}.DEGs.txt")
                        print(f"{a} {groups[i]}_{r1}", file=f)
                        print(f"{b} {groups[j]}_{r2}", file=f)
    print(f"[getGroups] {len(groups)} 处理组 -> {len(groups) * (len(groups) - 1) // 2} 组合")


if __name__ == "__main__":
    main()
