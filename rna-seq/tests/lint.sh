#!/bin/bash
#########################################################################
# Static checks: shell syntax/shellcheck, Python compilation, R parsing,
# and snakemake --lint (4 pipelines). Missing optional tools are skipped automatically with a note.
# Usage: bash tests/lint.sh      (used by both CI and local pre-checks)
#########################################################################
set -uo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR" || exit 1

FAIL=0
note() { echo "[lint] $*"; }

note "1/5 bash syntax (bash -n)"
for f in run.sh workflow/scripts/*.sh tests/*.sh; do
    if bash -n "$f"; then note "  OK  $f"; else note "  FAIL $f"; FAIL=1; fi
done

note "2/5 shellcheck (skipped if not installed)"
if command -v shellcheck >/dev/null 2>&1; then
    for f in run.sh workflow/scripts/*.sh tests/*.sh; do
        if shellcheck -S warning "$f"; then note "  OK  $f"; else note "  FAIL $f"; FAIL=1; fi
    done
else
    note "  SKIP shellcheck not installed"
fi

note "3/5 Python compilation"
if python3 -m py_compile workflow/scripts/*.py; then
    note "  OK  workflow/scripts/*.py"
else
    note "  FAIL python compilation"; FAIL=1
fi

note "4/5 R parsing (skipped if R is not installed)"
if command -v Rscript >/dev/null 2>&1; then
    for f in workflow/scripts/*.R; do
        if Rscript -e "invisible(parse('$f'))" >/dev/null 2>&1; then
            note "  OK  $f"
        else
            note "  FAIL $f"; FAIL=1
        fi
    done
else
    note "  SKIP Rscript not installed"
fi

note "5/5 snakemake --lint (skipped if not installed; 4 pipelines)"
if command -v snakemake >/dev/null 2>&1; then
    for p in upstream deg as lncrna; do
        RNASEQ_PIPELINE="$p" snakemake -s workflow/Snakefile --lint >/tmp/lint_$p.txt 2>&1 || true
        unexpected=$(grep -E '^    \* ' /tmp/lint_$p.txt | grep -v 'Specify a conda environment or container for each rule.:' || true)
        if [[ -z "$unexpected" ]]; then
            note "  OK  pipeline=$p (external-runtime conda warnings ignored)"
        else
            note "  FAIL pipeline=$p (see /tmp/lint_$p.txt)"; FAIL=1
        fi
    done
else
    note "  SKIP snakemake not installed"
fi

if [[ "$FAIL" == "0" ]]; then
    note "all checks passed"
else
    note "some checks failed"
fi
exit "$FAIL"
