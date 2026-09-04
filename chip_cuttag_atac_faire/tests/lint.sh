#!/bin/bash
#########################################################################
# 静态检查套件（移植自 rna-seq v0.8.0，对应实施计划 Phase 4 Task 4.1）：
# bash 语法/shellcheck（可选）、Python 编译、R 解析（可选）、YAML 解析、
# snakemake --lint。缺少的可选工具自动跳过并提示；CI 全量执行。
# 用法: bash tests/lint.sh      （CI 与本地 pre-check 通用）
#########################################################################
set -uo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR" || exit 1

FAIL=0
note() { echo "[lint] $*"; }

# python 解析：优先 python3，退回 python（Windows Git Bash 两者常缺一）
PYTHON_BIN=""
if command -v python3 >/dev/null 2>&1; then PYTHON_BIN=python3
elif command -v python >/dev/null 2>&1; then PYTHON_BIN=python
fi

note "1/6 bash 语法（bash -n）"
for f in run.sh tests/*.sh; do
    if bash -n "$f"; then note "  OK  $f"; else note "  FAIL $f"; FAIL=1; fi
done

note "2/6 shellcheck（未安装则跳过）"
if command -v shellcheck >/dev/null 2>&1; then
    for f in run.sh tests/*.sh; do
        if shellcheck -S warning "$f"; then note "  OK  $f"; else note "  FAIL $f"; FAIL=1; fi
    done
else
    note "  SKIP shellcheck 未安装"
fi

note "3/6 Python 编译"
if [[ -z "$PYTHON_BIN" ]]; then
    note "  FAIL 未找到 python/python3"; FAIL=1
elif "$PYTHON_BIN" -m py_compile tests/*.py workflow/scripts/*.py; then
    note "  OK  tests/*.py workflow/scripts/*.py"
else
    note "  FAIL python 编译"; FAIL=1
fi

note "4/6 R 解析（未安装 R 则跳过）"
if command -v Rscript >/dev/null 2>&1; then
    for f in workflow/scripts/*.R; do
        if Rscript --vanilla -e "invisible(parse(file='$f'))" >/dev/null 2>&1; then
            note "  OK  $f"
        else
            note "  FAIL $f"; FAIL=1
        fi
    done
else
    note "  SKIP Rscript 未安装"
fi

note "5/6 YAML 解析（缺 python/pyyaml 则跳过）"
if [[ -z "$PYTHON_BIN" ]]; then
    note "  SKIP 未找到 python/python3"
elif ! "$PYTHON_BIN" -c "import yaml" >/dev/null 2>&1; then
    note "  SKIP pyyaml 未安装"
else
    for f in config/*.yaml workflow/environment.yaml workflow/profile/*/config.yaml workflow/multiqc_config.yaml; do
        if "$PYTHON_BIN" -c "import sys,yaml; yaml.safe_load(open(sys.argv[1], encoding='utf-8'))" "$f" >/dev/null 2>&1; then
            note "  OK  $f"
        else
            note "  FAIL $f"; FAIL=1
        fi
    done
fi

note "6/6 snakemake --lint（未安装则跳过；单 pipeline 跑一次）"
if command -v snakemake >/dev/null 2>&1; then
    LINT_OUT=/tmp/lint_chip.txt
    snakemake --lint -s workflow/Snakefile --configfile config/config.yaml >"$LINT_OUT" 2>&1 || true
    # 已改为 external-runtime（无 per-rule conda），仅该类告警豁免，其余告警视为失败
    unexpected=$(grep -E '^    \* ' "$LINT_OUT" | grep -v 'Specify a conda environment or container for each rule.:' || true)
    if [[ -z "$unexpected" ]]; then
        note "  OK  workflow/Snakefile (external-runtime conda warnings ignored)"
    else
        note "  FAIL snakemake --lint（见 $LINT_OUT）"; FAIL=1
    fi
else
    note "  SKIP snakemake 未安装"
fi

if [[ "$FAIL" == "0" ]]; then
    note "全部检查通过 ✔"
else
    note "存在失败项 ✘"
fi
exit "$FAIL"
