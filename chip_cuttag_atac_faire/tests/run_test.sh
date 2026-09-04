#!/bin/bash
#########################################################################
# One-command regression test (implementation plan Phase 4 Task 4.3;
# skeleton ported from the rna-seq script of the same name):
#   generate synthetic test data -> assemble the test working directory
#   -> dry-run (default) or --real-run end-to-end -> assert -> clean up
#
# Usage:
#   bash tests/run_test.sh [--reads N] [--keep] [--real-run] [--help]
#     --reads     PE read pairs per sample, default 50000 (CI passes 2000)
#     --keep      keep tests/data and tests/work (cleaned up by default)
#     --real-run  run end-to-end and assert that outputs exist (default is
#                 a dry-run that only validates DAG integrity; a real run
#                 needs a full analysis environment with bowtie2/fastqc/
#                 trim_galore/macs2/R etc.; used for server-side validation)
# Requires: dry-run only needs snakemake + python3(+pyyaml); --real-run
# needs a full analysis environment.
#########################################################################
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TESTS_DIR="$REPO_DIR/tests"
DATA_DIR="$TESTS_DIR/data"
WORK_DIR="$TESTS_DIR/work"

usage() {
    cat <<'EOF'
One-command regression test: synthetic data -> assemble working directory
-> dry-run (default) / --real-run end-to-end -> assertions

Usage:
  bash tests/run_test.sh [--reads N] [--keep] [--real-run] [--help]
    --reads     PE read pairs per sample, default 50000 (CI passes 2000)
    --keep      keep tests/data and tests/work (cleaned up by default)
    --real-run  run end-to-end and assert that outputs exist (default only
                dry-runs to validate DAG integrity; a real run needs a full
                analysis environment; used for server-side validation)
    -h, --help  show this help

Requires:
  dry-run only needs snakemake + python3(+pyyaml); --real-run needs a full
  analysis environment (bowtie2/fastqc/trim_galore/macs2/deeptools/R etc.,
  see workflow/environment.yaml).
EOF
}

# ---------------------------------------------------------------------
# Argument parsing: dry-run by default; --real-run enables the end-to-end
# branch
# ---------------------------------------------------------------------
READS=50000
KEEP=0
REAL_RUN=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --reads)
            [[ $# -ge 2 ]] || { echo "[ERROR] --reads requires a value" >&2; exit 1; }
            [[ "$2" =~ ^[0-9]+$ ]] || { echo "[ERROR] --reads must be a positive integer: $2" >&2; exit 1; }
            READS="$2"; shift 2 ;;
        --keep)     KEEP=1; shift ;;
        --real-run) REAL_RUN=1; shift ;;
        -h|--help)  usage; exit 0 ;;
        *) echo "[ERROR] Unknown argument: $1 (see --help for usage)" >&2; exit 1 ;;
    esac
done
MODE="dry-run"
[[ "$REAL_RUN" == 1 ]] && MODE="real-run"

echo "[test] 1/5 Checking dependencies (mode=$MODE, reads=$READS)"
command -v snakemake >/dev/null || { echo "[ERROR] snakemake not found" >&2; exit 1; }
command -v python3   >/dev/null || { echo "[ERROR] python3 not found" >&2; exit 1; }
# run.sh needs pyyaml to resolve software.yaml at runtime; fail early with a clear message
python3 -c 'import yaml' >/dev/null 2>&1 || { echo "[ERROR] python3 is missing pyyaml" >&2; exit 1; }

echo "[test] 2/5 Generating synthetic test data (reads=$READS, fixed seed)"
python3 "$TESTS_DIR/make_testdata.py" --outdir "$DATA_DIR" --reads "$READS"

echo "[test] 3/5 Assembling test project tests/work"
rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR"
# Prefer symlinking 1.rawdata (saves time and space); fall back to copying
# when symlinks are unavailable (e.g. default Windows configuration)
if ln -s "$DATA_DIR/1.rawdata" "$WORK_DIR/1.rawdata" 2>/dev/null; then
    echo "[test]   1.rawdata: linked via symlink"
else
    mkdir -p "$WORK_DIR/1.rawdata"
    cp "$DATA_DIR"/1.rawdata/*.fq.gz "$WORK_DIR/1.rawdata/"
    echo "[test]   1.rawdata: copied in (symlink unavailable)"
fi
cp -r "$DATA_DIR/ref" "$WORK_DIR/ref"
cp "$DATA_DIR/samples.csv" "$DATA_DIR/config.yaml" "$WORK_DIR/"

echo "[test] 4/5 Running the workflow ($MODE)"
RUN_ARGS=(-P "$WORK_DIR" -c "$WORK_DIR/config.yaml" -j 4)
# dry-run (default) adds -n: builds the DAG without executing; --real-run
# drops -n for a real end-to-end run
[[ "$REAL_RUN" == 1 ]] || RUN_ARGS+=(-n)
CAPTURE_LOG="$WORK_DIR/run_test.capture.log"
set +e
bash "$REPO_DIR/run.sh" "${RUN_ARGS[@]}" >"$CAPTURE_LOG" 2>&1
STATUS=$?
set -e
if [[ "$STATUS" != 0 ]]; then
    echo "[ERROR] run.sh exited with code $STATUS (tail of output below)" >&2
    tail -40 "$CAPTURE_LOG" >&2 || true
    exit 1
fi

echo "[test] 5/5 Assertions"
FAIL=0
if [[ "$REAL_RUN" == 1 ]]; then
    # ---------- real-run output existence assertions (5 samples: 3 chip + 2 atac) ----------
    shopt -s nullglob
    bams=("$WORK_DIR"/results/3.align/bowtie2/*_sorted.bam)
    shopt -u nullglob
    if [[ ${#bams[@]} -eq 5 ]]; then
        echo "  PASS  results/3.align/bowtie2/*_sorted.bam count ${#bams[@]} (expected 5)"
    else
        echo "  FAIL  results/3.align/bowtie2/*_sorted.bam count ${#bams[@]} (expected 5)"; FAIL=1
    fi
    EXPECTED=(
        "results/4.peak/g1_peaks.narrowPeak"
        "results/4.peak/g2_peaks.narrowPeak"
        "results/5.QC/frip/FRiP_summary.tsv"
        "results/2.cleandata/fastqc/multiqc/multiqc_report.html"
        "results/5.QC/software_versions.yaml"
    )
    for rel in "${EXPECTED[@]}"; do
        if [[ -s "$WORK_DIR/$rel" ]]; then
            echo "  PASS  $rel"
        else
            echo "  FAIL  $rel (missing or empty)"; FAIL=1
        fi
    done
else
    # ---------- dry-run DAG assertions: exit code 0 (checked above) and three key rule names present ----------
    SNAKE_LOG="$WORK_DIR/snakemake.logs.txt"   # run.sh default log (stdout is tee'd to disk)
    for rule in bowtie2_mapping callpeak_narrow callpeak_atac; do
        if grep -q "$rule" "$CAPTURE_LOG" "$SNAKE_LOG" 2>/dev/null; then
            echo "  PASS  DAG contains rule $rule"
        else
            echo "  FAIL  DAG does not contain rule $rule"; FAIL=1
        fi
    done
fi
if [[ "$FAIL" != "0" ]]; then
    echo "[ERROR] Assertions failed; keeping workspace for debugging: $WORK_DIR" >&2
    exit 1
fi

## DAG regeneration (best effort: only runs when graphviz is installed;
## failure does not affect the test result)
if command -v dot >/dev/null 2>&1; then
    (cd "$WORK_DIR" && CHIP_CONFIG="$WORK_DIR/config.yaml" \
        snakemake -s "$REPO_DIR/workflow/Snakefile" --configfile "$WORK_DIR/config.yaml" --dag \
        | dot -Tsvg -o "$REPO_DIR/docs/dag_test.svg") \
        && echo "[test] DAG diagram updated: docs/dag_test.svg" \
        || echo "[test] [WARN] DAG generation failed (skipped; does not affect the test result)"
else
    echo "[test] [WARN] graphviz not installed; skipping DAG regeneration"
fi

echo "[test] All checks passed (mode=$MODE reads=$READS)"
if [[ "$KEEP" != "1" ]]; then
    rm -rf "$WORK_DIR" "$DATA_DIR"
    echo "[test] Cleaned tests/work and tests/data (use --keep to keep them)"
fi
