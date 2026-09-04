#!/bin/bash
#########################################################################
# Static checks: shell syntax/shellcheck, Python compilation, and
# snakemake --lint. Missing optional tools are skipped automatically with
# a note. Usage: bash tests/lint.sh (used by both CI and local pre-checks)
#########################################################################
set -uo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR" || exit 1

FAIL=0
note() { echo "[lint] $*"; }

note "1/4 bash syntax (bash -n)"
for f in run.sh tests/*.sh; do
    if bash -n "$f"; then note "  OK  $f"; else note "  FAIL $f"; FAIL=1; fi
done

note "2/4 shellcheck (skipped if not installed)"
if command -v shellcheck >/dev/null 2>&1; then
    for f in run.sh tests/*.sh; do
        if shellcheck -S warning "$f"; then note "  OK  $f"; else note "  FAIL $f"; FAIL=1; fi
    done
else
    note "  SKIP shellcheck not installed"
fi

note "3/4 Python compilation"
if python3 -m py_compile workflow/scripts/*.py tests/make_testdata.py; then
    note "  OK  workflow/scripts/*.py tests/make_testdata.py"
else
    note "  FAIL python compilation"; FAIL=1
fi

note "4/4 snakemake --lint (skipped if not installed)"
if command -v snakemake >/dev/null 2>&1; then
    snakemake -s workflow/Snakefile --lint >/tmp/lint_seclip.txt 2>&1 || true
    unexpected=$(grep -E '^    \* ' /tmp/lint_seclip.txt | grep -v 'Specify a conda environment or container for each rule.:' || true)
    if [[ -z "$unexpected" ]]; then
        note "  OK  snakemake --lint (external-runtime conda warnings ignored)"
    else
        note "  FAIL snakemake --lint (see /tmp/lint_seclip.txt)"; FAIL=1
    fi
else
    note "  SKIP snakemake not installed"
fi

if [[ "$FAIL" == "0" ]]; then
    note "all checks passed"
else
    note "some checks failed"
fi
exit "$FAIL"
