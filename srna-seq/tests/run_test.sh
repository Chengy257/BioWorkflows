#!/bin/bash
#########################################################################
# One-command regression test (ported from the chip_cuttag_atac_faire
# script of the same name): generate synthetic test data -> assemble the
# test working directory -> dry-run (default) or --real-run end-to-end
# -> assert -> clean up
#
# Usage:
#   bash tests/run_test.sh [--reads N] [--keep] [--real-run] [--deg] [--help]
#     --reads     SE reads per sample, default 50000 (CI passes 2000)
#     --keep      keep tests/data and tests/work (cleaned up by default)
#     --real-run  run end-to-end and assert that outputs exist (default is
#                 a dry-run that only validates DAG integrity; a real run
#                 needs a full analysis environment with STAR/umi-tools/
#                 cutadapt/pureclip etc.)
#     --deg       dry-run the optional differential-expression scenario:
#                 sample table gains a group column (s1=control, s2=treat)
#                 and the deg stage is enabled for the miRNA class; asserts
#                 that the DAG contains the deg_deseq2 rule
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
  bash tests/run_test.sh [--reads N] [--keep] [--real-run] [--deg] [--help]
  --reads     SE reads per sample, default 50000
  --keep      keep tests/data and tests/work (cleaned up by default)
  --real-run  run end-to-end and assert that outputs exist (default only
              dry-runs to validate DAG integrity; a real run needs a full
              analysis environment; used for server validation)
  --deg       dry-run the optional differential-expression scenario (group
              column in the sample table, deg enabled for miRNA); asserts
              that the DAG contains the deg_deseq2 rule
  -h, --help  show this help

Requires:
  dry-run only needs snakemake + python3(+pyyaml); --real-run needs a full
  analysis environment (bowtie/trim-galore/fastqc/multiqc, see
  workflow/environment.yaml).
EOF
}

READS=50000
KEEP=0
REAL_RUN=0
DEG=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --reads)
            [[ $# -ge 2 ]] || { echo "[ERROR] --reads requires a value" >&2; exit 1; }
            [[ "$2" =~ ^[0-9]+$ ]] || { echo "[ERROR] --reads must be a positive integer: $2" >&2; exit 1; }
            READS="$2"; shift 2 ;;
        --keep)     KEEP=1; shift ;;
        --real-run) REAL_RUN=1; shift ;;
        --deg)      DEG=1; shift ;;
        -h|--help)  usage; exit 0 ;;
        *) echo "[ERROR] Unknown argument: $1 (see --help for usage)" >&2; exit 1 ;;
    esac
done
if [[ "$DEG" == 1 && "$REAL_RUN" == 1 ]]; then
    echo "[ERROR] --deg is a dry-run-only scenario (the deg stage needs R/DESeq2); drop --real-run" >&2
    exit 1
fi
MODE="dry-run"
[[ "$REAL_RUN" == 1 ]] && MODE="real-run"
[[ "$DEG" == 1 ]] && MODE="${MODE}+deg"

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
    cp "$DATA_DIR"/1.rawdata/*.fq.gz "$WORK_DIR/1.rawdata/"
    echo "[test]   1.rawdata: copied in (symlink unavailable)"
fi
cp -r "$DATA_DIR/ref" "$WORK_DIR/ref"
cp "$DATA_DIR/samples.csv" "$DATA_DIR/config.yaml" "$WORK_DIR/"
if [[ "$DEG" == 1 ]]; then
    # Differential-expression scenario: the sample table gains a group column
    # (single replicate per group -- legal, warns only) and the deg stage is
    # enabled for the miRNA class (configured in the generated cascade).
    printf 'sample_id,group\ns1,control\ns2,treat\n' > "$WORK_DIR/samples.csv"
    cat >> "$WORK_DIR/config.yaml" <<'EOF'
deg:
  enabled: true
  classes: ["miRNA"]
  control_group: "control"
  foldchange: 2
  padj: 0.05
  batch_correction: "F"
  pca_ntop: 2000
EOF
    echo "[test]   deg scenario: group column written, deg stage enabled (miRNA)"
fi

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
        "results/4.expression/miRNA/miRNA_counts.tsv"
        "results/4.expression/miRNA/miRNA_RPM.tsv"
        "results/4.expression/all_classes_counts.tsv"
        "results/5.QC/cascade_summary.tsv"
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
    RULES=(cascade_stage genome_align merge_counts)
    if [[ "$DEG" == 1 ]]; then
        RULES+=(deg_deseq2)
    fi
    for rule in "${RULES[@]}"; do
        if grep -q "$rule" "$CAPTURE_LOG" "$SNAKE_LOG" 2>/dev/null; then
            echo "  PASS  DAG contains rule $rule"
        else
            echo "  FAIL  DAG does not contain rule $rule"; FAIL=1
        fi
    done
    if [[ "$DEG" != 1 ]]; then
        # Default-off contract: the deg stage must be absent from the DAG.
        if grep -q "deg_deseq2" "$CAPTURE_LOG" "$SNAKE_LOG" 2>/dev/null; then
            echo "  FAIL  default DAG must not contain rule deg_deseq2"; FAIL=1
        else
            echo "  PASS  default DAG does not contain rule deg_deseq2"
        fi
    fi
    JOB_TOTAL="$(grep -hE '^total[[:space:]]+[0-9]+' "$CAPTURE_LOG" "$SNAKE_LOG" 2>/dev/null | tail -n1 | awk '{print $NF}')"
    if [[ -n "$JOB_TOTAL" ]]; then
        echo "  INFO  dry-run job total: $JOB_TOTAL"
    else
        echo "  INFO  dry-run job total: not found in the launcher output"
    fi
fi
if [[ "$FAIL" != "0" ]]; then
    echo "[ERROR] Assertions failed; keeping workspace for debugging: $WORK_DIR" >&2
    exit 1
fi

if [[ "$DEG" == 1 ]]; then
    echo "[test] [INFO] --deg scenario: docs/dag_test.svg regeneration skipped (keeps the default-scenario DAG)"
elif command -v dot >/dev/null 2>&1; then
    (cd "$WORK_DIR" && SRNA_CONFIG="$WORK_DIR/config.yaml" \
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
