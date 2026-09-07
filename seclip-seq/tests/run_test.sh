#!/bin/bash
#########################################################################
# One-command regression test (ported from the chip_cuttag_atac_faire
# script of the same name): generate synthetic test data -> assemble the
# test working directory -> dry-run (default) or --real-run end-to-end
# -> assert -> clean up
#
# Usage:
#   bash tests/run_test.sh [--reads N] [--keep] [--real-run] [--consensus]
#                          [--input-control] [--help]
#     --reads      SE reads per sample, default 50000 (CI passes 2000)
#     --keep       keep tests/data and tests/work (cleaned up by default)
#     --real-run   run end-to-end and assert that outputs exist (default is
#                  a dry-run that only validates DAG integrity; a real run
#                  needs a full analysis environment with STAR/umi-tools/
#                  cutadapt/pureclip etc.)
#     --consensus  optional-stage scenario: switch the sample table to the
#                  condition/role form (both samples ip, same condition) and
#                  enable reproducible_peaks + annotate_peaks, then dry-run
#                  and assert the new rules join the DAG (a --consensus
#                  --real-run additionally needs bedtools in the environment)
#     --input-control  optional-stage scenario: 4 samples (2 ip + 2 input
#                  controls, same condition) with reproducible_peaks enabled
#                  including input_control + filter_by_input, then dry-run
#                  and assert the background/flagging rules join the DAG (an
#                  --input-control --real-run additionally needs bedtools)
# Requires: dry-run only needs snakemake + python3(+pyyaml); --real-run
# needs a full analysis environment.
#########################################################################
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TESTS_DIR="$REPO_DIR/tests"
DATA_DIR="$TESTS_DIR/data"
# The assembled working directory defaults to /tmp (Linux filesystem):
# STAR finishes writing some outputs (e.g. the unmapped Fastx) asynchronously,
# and on WSL DrvFs mounts (/mnt/*, 9p protocol) a reader that opens the file
# within seconds of the writer exiting can observe an empty view until the
# host flushes. Keep the default overridable for debugging.
WORK_DIR="${SECLIP_TEST_WORK_DIR:-/tmp/seclip_seq_test_work}"

usage() {
    cat <<'EOF'
One-command regression test: synthetic data -> assemble working directory
-> dry-run (default) / --real-run end-to-end -> assertions

Usage:
  bash tests/run_test.sh [--reads N] [--keep] [--real-run] [--consensus]
                         [--input-control] [--help]
  --reads      SE reads per sample, default 50000
  --keep       keep tests/data and tests/work (cleaned up by default)
  --real-run   run end-to-end and assert that outputs exist (default only
               dry-runs to validate DAG integrity; a real run needs a full
               analysis environment; used for server validation)
  --consensus  optional-stage scenario: condition/role sample table plus
               reproducible_peaks/annotate_peaks enabled, dry-run asserts
               the new rules join the DAG
  --input-control  optional-stage scenario: 2 ip + 2 input samples, one
               condition, reproducible_peaks enabled with input_control and
               filter_by_input, dry-run asserts the background/flagging
               rules join the DAG
  -h, --help   show this help

Requires:
  dry-run only needs snakemake + python3(+pyyaml); --real-run needs a full
  analysis environment (STAR/umi-tools/cutadapt/seqkit/samtools/fastqc/
  multiqc/pureclip, see workflow/environment.yaml); --consensus /
  --input-control --real-run additionally need bedtools.
EOF
}

READS=50000
KEEP=0
REAL_RUN=0
CONSENSUS=0
INPUT_CONTROL=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --reads)
            [[ $# -ge 2 ]] || { echo "[ERROR] --reads requires a value" >&2; exit 1; }
            [[ "$2" =~ ^[0-9]+$ ]] || { echo "[ERROR] --reads must be a positive integer: $2" >&2; exit 1; }
            READS="$2"; shift 2 ;;
        --keep)      KEEP=1; shift ;;
        --real-run)  REAL_RUN=1; shift ;;
        --consensus) CONSENSUS=1; shift ;;
        --input-control) INPUT_CONTROL=1; shift ;;
        -h|--help)   usage; exit 0 ;;
        *) echo "[ERROR] Unknown argument: $1 (see --help)" >&2; exit 1 ;;
    esac
done
if [[ "$CONSENSUS" == 1 && "$INPUT_CONTROL" == 1 ]]; then
    echo "[ERROR] --consensus and --input-control are separate scenarios; pass one at a time" >&2
    exit 1
fi
MODE="dry-run"
[[ "$REAL_RUN" == 1 ]] && MODE="real-run"
[[ "$CONSENSUS" == 1 ]] && MODE="$MODE+consensus"
[[ "$INPUT_CONTROL" == 1 ]] && MODE="$MODE+input_control"

echo "[test] 1/5 Checking dependencies (mode=$MODE, reads=$READS)"
command -v snakemake >/dev/null || { echo "[ERROR] snakemake not found" >&2; exit 1; }
command -v python3   >/dev/null || { echo "[ERROR] python3 not found" >&2; exit 1; }
python3 -c 'import yaml' >/dev/null 2>&1 || { echo "[ERROR] python3 is missing pyyaml" >&2; exit 1; }

echo "[test] 2/5 Generating synthetic test data (reads=$READS, fixed seed)"
GEN_ARGS=(--outdir "$DATA_DIR" --reads "$READS")
[[ "$INPUT_CONTROL" == 1 ]] && GEN_ARGS+=(--with-inputs)
python3 "$TESTS_DIR/make_testdata.py" "${GEN_ARGS[@]}"

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

if [[ "$CONSENSUS" == 1 ]]; then
    # Optional-stage scenario (v0.2): switch the sample table to the
    # condition/role form (both samples in one ip condition) and enable
    # reproducible_peaks + annotate_peaks on top of the generated config.
    # The default (no-flag) run keeps the generated single-column table and
    # both stages disabled.
    cat > "$WORK_DIR/samples.csv" <<'EOF'
sample_id,condition,role
FC_rep1,treatment,ip
FC_rep2,treatment,ip
EOF
    cat >> "$WORK_DIR/config.yaml" <<'EOF'
reproducible_peaks:
  enabled: true
  min_replicates: 2
annotate_peaks:
  enabled: true
EOF
    echo "[test]   --consensus: condition/role sample table + reproducible_peaks/annotate_peaks enabled"
fi

if [[ "$INPUT_CONTROL" == 1 ]]; then
    # Optional-stage scenario (v0.2, input control): 2 ip + 2 input samples in
    # one condition, reproducible_peaks enabled with input_control and
    # filter_by_input. The default (no-flag) run keeps the generated
    # single-column table and the stage disabled.
    cat > "$WORK_DIR/samples.csv" <<'EOF'
sample_id,condition,role
FC_rep1,FC,ip
FC_rep2,FC,ip
FC_in1,FC,input
FC_in2,FC,input
EOF
    cat >> "$WORK_DIR/config.yaml" <<'EOF'
reproducible_peaks:
  enabled: true
  min_replicates: 2
  input_control: true
  filter_by_input: true
EOF
    echo "[test]   --input-control: 2 ip + 2 input samples + reproducible_peaks with input_control/filter_by_input"
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
        "results/5.QC/multiqc/multiqc_report.html"
        "results/5.QC/software_versions.yaml"
        "results/5.callpeak/FC_rep1.pureclip.bed"
        "results/4.rmdup/FC_rep1_readnum.txt"
    )
    if [[ "$CONSENSUS" == 1 ]]; then
        EXPECTED+=(
            "results/6.reproducible_peaks/treatment.consensus.bed"
            "results/6.annotation/FC_rep1.annotation.tsv"
            "results/6.annotation/treatment.consensus.annotation.tsv"
        )
    fi
    if [[ "$INPUT_CONTROL" == 1 ]]; then
        EXPECTED+=(
            "results/6.reproducible_peaks/FC.input_background.bed"
            "results/6.reproducible_peaks/FC.consensus.bed"
            "results/6.reproducible_peaks/FC.consensus.filtered.bed"
        )
    fi
    for rel in "${EXPECTED[@]}"; do
        if [[ -s "$WORK_DIR/$rel" ]]; then
            echo "  PASS  $rel"
        else
            echo "  FAIL  $rel (missing or empty)"; FAIL=1
        fi
    done
else
    SNAKE_LOG="$WORK_DIR/snakemake.logs.txt"
    RULES=(star_align umi_dedup callpeak_pureclip)
    if [[ "$CONSENSUS" == 1 ]]; then
        RULES+=(consensus_peaks gtf_gene_regions annotate_sample_peaks annotate_consensus_peaks)
    fi
    if [[ "$INPUT_CONTROL" == 1 ]]; then
        RULES+=(consensus_peaks_raw input_background flag_input_background filter_input_background)
    fi
    for rule in "${RULES[@]}"; do
        # Word-boundary match so consensus_peaks is not satisfied by the
        # annotate_consensus_peaks table row (underscore is a word character).
        if grep -qE "(^|[^A-Za-z0-9_])${rule}([^A-Za-z0-9_]|$)" "$CAPTURE_LOG" "$SNAKE_LOG" 2>/dev/null; then
            echo "  PASS  DAG contains rule $rule"
        else
            echo "  FAIL  DAG does not contain rule $rule"; FAIL=1
        fi
    done
    total=$(grep -hE '^total[[:space:]]+[0-9]+' "$CAPTURE_LOG" "$SNAKE_LOG" 2>/dev/null | awk '{print $NF}' | tail -1)
    if [[ -n "$total" ]]; then
        echo "  INFO  dry-run job total: $total (baselines: default 23, --consensus 28, --input-control 45)"
    else
        echo "  INFO  dry-run job total: not found in the logs"
    fi
fi
if [[ "$FAIL" != "0" ]]; then
    echo "[ERROR] Assertions failed; keeping workspace for debugging: $WORK_DIR" >&2
    exit 1
fi

if command -v dot >/dev/null 2>&1; then
    (cd "$WORK_DIR" && SECLIP_CONFIG="$WORK_DIR/config.yaml" \
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
