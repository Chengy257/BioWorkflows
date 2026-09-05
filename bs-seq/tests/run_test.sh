#!/bin/bash
#########################################################################
# One-command regression test (ported from the chip_cuttag_atac_faire
# script of the same name): generate synthetic test data -> assemble the
# test working directory -> dry-run (default) or --real-run end-to-end
# -> assert -> clean up
#
# Usage:
#   bash tests/run_test.sh [--reads N] [--keep] [--real-run] [--help]
#     --reads     PE read pairs per sample, default 50000
#     --keep      keep tests/data and tests/work (cleaned up by default)
#     --real-run  run end-to-end and assert that outputs exist (default is
#                 a dry-run that only validates DAG integrity; a real run
#                 needs a full analysis environment with bismark/bowtie2/
#                 samtools/trim-galore/fastqc/multiqc)
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
  --reads     PE read pairs per sample, default 50000
  --keep      keep tests/data and tests/work (cleaned up by default)
  --real-run  run end-to-end and assert that outputs exist (default only
              dry-runs to validate DAG integrity; a real run needs a full
              analysis environment; used for server validation)
  -h, --help  show this help

Requires:
  dry-run only needs snakemake + python3(+pyyaml); --real-run needs a full
  analysis environment (bismark/bowtie2/samtools/trim-galore/fastqc/
  multiqc, see workflow/environment.yaml).
EOF
}

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
python3 -c 'import yaml' >/dev/null 2>&1 || { echo "[ERROR] python3 is missing pyyaml" >&2; exit 1; }

echo "[test] 2/5 Generating synthetic test data (reads=$READS, fixed seed)"
python3 "$TESTS_DIR/make_testdata.py" --outdir "$DATA_DIR" --reads "$READS"

echo "[test] 3/5 Assembling test project tests/work"
rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR"
if ln -s "$DATA_DIR/1.rawdata" "$WORK_DIR/1.rawdata" 2>/dev/null; then
    echo "[test]   1.rawdata: linked via symlink"
else
    mkdir -p "$WORK_DIR/1.rawdata"
    cp "$DATA_DIR"/1.rawdata/*.fq.gz "$WORK_DIR/1.rawdata/" 2>/dev/null \
        || cp "$DATA_DIR"/1.rawdata/*.fastq.gz "$WORK_DIR/1.rawdata/"
    echo "[test]   1.rawdata: copied in (symlink unavailable)"
fi
cp -r "$DATA_DIR/ref" "$WORK_DIR/ref"
cp "$DATA_DIR/samples.csv" "$DATA_DIR/config.yaml" "$WORK_DIR/"

echo "[test] 4/5 Running the workflow ($MODE)"
RUN_ARGS=(-P "$WORK_DIR" -c "$WORK_DIR/config.yaml" -j 4)
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
    EXPECTED=(
        "results/5.methylation/s1/s1.bam.deduplicated.bismark.cov.gz"
        "results/5.QC/bismark2summary.html"
        "results/5.QC/multiqc/multiqc_report.html"
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
    SNAKE_LOG="$WORK_DIR/snakemake.logs.txt"
    for rule in bismark_align deduplicate methylation_extractor coverage2cytosine; do
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

if command -v dot >/dev/null 2>&1; then
    (cd "$WORK_DIR" && BSSEQ_CONFIG="$WORK_DIR/config.yaml" \
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
