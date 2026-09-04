#!/bin/bash
#########################################################################
# 静态检查（对应路线图 3.2）：shell 语法/shellcheck、Python 编译、R 解析、
# snakemake --lint（4 种 pipeline）。缺少的可选工具自动跳过并提示。
# 用法: bash tests/lint.sh      （CI 与本地 pre-check 通用）
#########################################################################
set -uo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR"

FAIL=0
note() { echo "[lint] $*"; }

note "1/5 bash 语法（bash -n）"
for f in run.sh workflow/scripts/*.sh tests/*.sh; do
    if bash -n "$f"; then note "  OK  $f"; else note "  FAIL $f"; FAIL=1; fi
done

note "2/5 shellcheck（未安装则跳过）"
if command -v shellcheck >/dev/null 2>&1; then
    for f in run.sh workflow/scripts/*.sh tests/*.sh; do
        if shellcheck -S warning "$f"; then note "  OK  $f"; else note "  FAIL $f"; FAIL=1; fi
    done
else
    note "  SKIP shellcheck 未安装"
fi

note "3/5 Python 编译"
if python3 -m py_compile workflow/scripts/*.py; then
    note "  OK  workflow/scripts/*.py"
else
    note "  FAIL python 编译"; FAIL=1
fi

note "4/5 R 解析（未安装 R 则跳过）"
if command -v Rscript >/dev/null 2>&1; then
    for f in workflow/scripts/*.R; do
        if Rscript -e "invisible(parse('$f'))" >/dev/null 2>&1; then
            note "  OK  $f"
        else
            note "  FAIL $f"; FAIL=1
        fi
    done
else
    note "  SKIP Rscript 未安装"
fi

note "5/5 snakemake --lint（未安装则跳过；4 种 pipeline）"
if command -v snakemake >/dev/null 2>&1; then
    for p in upstream deg as lncrna; do
        if RNASEQ_PIPELINE="$p" snakemake -s workflow/Snakefile --lint >/tmp/lint_$p.txt 2>&1; then
            note "  OK  pipeline=$p"
        else
            note "  FAIL pipeline=$p（见 /tmp/lint_$p.txt）"; FAIL=1
        fi
    done
else
    note "  SKIP snakemake 未安装"
fi

if [[ "$FAIL" == "0" ]]; then
    note "全部检查通过 ✔"
else
    note "存在失败项 ✘"
fi
exit "$FAIL"
