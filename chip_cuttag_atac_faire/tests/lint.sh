#!/bin/bash
#########################################################################
# Static check suite (ported from rna-seq v0.8.0; implementation plan
# Phase 4 Task 4.1): bash syntax/shellcheck (optional), Python compilation,
# R parsing (optional), YAML parsing, snakemake --lint. Missing optional
# tools are skipped with a note; CI runs everything.
# Usage: bash tests/lint.sh      (used by CI and local pre-checks alike)
#########################################################################
set -uo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR" || exit 1

FAIL=0
note() { echo "[lint] $*"; }

# python resolution: prefer python3, fall back to python (Git Bash on
# Windows often has only one of the two)
PYTHON_BIN=""
if command -v python3 >/dev/null 2>&1; then PYTHON_BIN=python3
elif command -v python >/dev/null 2>&1; then PYTHON_BIN=python
fi

note "1/6 bash syntax (bash -n)"
for f in run.sh tests/*.sh; do
    if bash -n "$f"; then note "  OK  $f"; else note "  FAIL $f"; FAIL=1; fi
done

note "2/6 shellcheck (skipped when not installed)"
if command -v shellcheck >/dev/null 2>&1; then
    for f in run.sh tests/*.sh; do
        if shellcheck -S warning "$f"; then note "  OK  $f"; else note "  FAIL $f"; FAIL=1; fi
    done
else
    note "  SKIP shellcheck not installed"
fi

note "3/6 Python compilation"
if [[ -z "$PYTHON_BIN" ]]; then
    note "  FAIL python/python3 not found"; FAIL=1
elif "$PYTHON_BIN" -m py_compile tests/*.py workflow/scripts/*.py; then
    note "  OK  tests/*.py workflow/scripts/*.py"
else
    note "  FAIL python compilation"; FAIL=1
fi

note "4/6 R parsing (skipped when R is not installed)"
if command -v Rscript >/dev/null 2>&1; then
    for f in workflow/scripts/*.R; do
        if Rscript --vanilla -e "invisible(parse(file='$f'))" >/dev/null 2>&1; then
            note "  OK  $f"
        else
            note "  FAIL $f"; FAIL=1
        fi
    done
else
    note "  SKIP Rscript not installed"
fi

note "5/6 YAML parsing (skipped when python/pyyaml is missing)"
if [[ -z "$PYTHON_BIN" ]]; then
    note "  SKIP python/python3 not found"
elif ! "$PYTHON_BIN" -c "import yaml" >/dev/null 2>&1; then
    note "  SKIP pyyaml not installed"
else
    for f in config/*.yaml workflow/environment.yaml workflow/profile/*/config.yaml workflow/multiqc_config.yaml; do
        if "$PYTHON_BIN" -c "import sys,yaml; yaml.safe_load(open(sys.argv[1], encoding='utf-8'))" "$f" >/dev/null 2>&1; then
            note "  OK  $f"
        else
            note "  FAIL $f"; FAIL=1
        fi
    done
fi

note "6/6 snakemake --lint (skipped when not installed; single pipeline, run once)"
if command -v snakemake >/dev/null 2>&1; then
    LINT_OUT=/tmp/lint_chip.txt
    snakemake --lint -s workflow/Snakefile --configfile config/config.yaml >"$LINT_OUT" 2>&1 || true
    # external-runtime setup (no per-rule conda): only that class of warning is
    # exempted; every other warning counts as a failure
    unexpected=$(grep -E '^    \* ' "$LINT_OUT" | grep -v 'Specify a conda environment or container for each rule.:' || true)
    if [[ -z "$unexpected" ]]; then
        note "  OK  workflow/Snakefile (external-runtime conda warnings ignored)"
    else
        note "  FAIL snakemake --lint (see $LINT_OUT)"; FAIL=1
    fi
else
    note "  SKIP snakemake not installed"
fi

if [[ "$FAIL" == "0" ]]; then
    note "All checks passed"
else
    note "Some checks failed"
fi
exit "$FAIL"
