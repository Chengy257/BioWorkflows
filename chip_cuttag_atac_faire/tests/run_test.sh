#!/bin/bash
#########################################################################
# One-command regression test (implementation plan Phase 4 Task 4.3;
# skeleton ported from the rna-seq script of the same name):
#   generate synthetic test data -> assemble the test working directory
#   -> dry-run (default) or --real-run end-to-end -> assert -> clean up
#
# Usage:
#   bash tests/run_test.sh [--reads N] [--keep] [--real-run] \
#       [--replicate] [--qc-full] [--motif] [--diffbind] [--gates] [--spike-in] [--seacr] [--help]
#     --reads     PE read pairs per sample, default 50000 (CI passes 2000)
#     --keep      keep tests/data and tests/work (cleaned up by default)
#     --real-run  run end-to-end and assert that outputs exist (default is
#                 a dry-run that only validates DAG integrity; a real run
#                 needs a full analysis environment with bowtie2/fastqc/
#                 trim_galore/macs2/R etc.; used for server-side validation)
#     --replicate scenario: add a 2-treat broad group and enable the
#                 replicate-aware peak stage (IDR + consensus)
#     --qc-full   scenario: enable tss/organelle QC and the synthetic
#                 blacklist (extended QC)
#     --seacr     scenario: switch the pooled peak caller to SEACR
#                 (peak.caller=seacr; a real run additionally needs the
#                 external SEACR script)
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
  bash tests/run_test.sh [--reads N] [--keep] [--real-run] [--replicate] [--qc-full] [--help]
    --reads     PE read pairs per sample, default 50000 (CI passes 2000)
    --keep      keep tests/data and tests/work (cleaned up by default)
    --real-run  run end-to-end and assert that outputs exist (default only
                dry-runs to validate DAG integrity; a real run needs a full
                analysis environment; used for server-side validation)
    --replicate scenario: 2-treat broad group + replicate-aware peak stage
                (per-replicate calling, IDR, consensus)
    --qc-full   scenario: tss/organelle QC + synthetic blacklist
    --motif     scenario: enable the HOMER motif stage (dummy genome tag;
                dry-run only — a real run needs an external HOMER install)
    --diffbind  scenario: 2-treat narrow group (condition/batch columns) +
                one DiffBind contrast (dry-run; a real run needs DiffBind)
    --gates     scenario: enable the QC gate summary stage (qc.gates;
                per-sample flagstat + PASS/WARN/FAIL gate table)
    --spike-in  scenario: enable the spike-in normalization stage (spike_in;
                synthetic spike-in reference + second-pass alignment)
    --seacr     scenario: pooled peak caller = SEACR (peak.caller=seacr;
                the pooled MACS2 caller rules are replaced and must be
                absent from the DAG; FRiP still consumes the peaks)
    -h, --help  show this help

Requires:
  dry-run only needs snakemake + python3(+pyyaml); --real-run needs a full
  analysis environment (bowtie2/fastqc/trim_galore/macs2/deeptools/R etc.,
  see workflow/environment.yaml; --replicate additionally needs the external
  idr tool, --motif an external HOMER, --diffbind the DiffBind R package,
  --seacr the external SEACR script).
EOF
}

# ---------------------------------------------------------------------
# Argument parsing: dry-run by default; --real-run enables the end-to-end
# branch; --replicate / --qc-full select scenario data/config
# ---------------------------------------------------------------------
READS=50000
KEEP=0
REAL_RUN=0
REPLICATE=0
QC_FULL=0
MOTIF=0
DIFFBIND=0
GATES=0
SPIKE=0
SEACR=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --reads)
            [[ $# -ge 2 ]] || { echo "[ERROR] --reads requires a value" >&2; exit 1; }
            [[ "$2" =~ ^[0-9]+$ ]] || { echo "[ERROR] --reads must be a positive integer: $2" >&2; exit 1; }
            READS="$2"; shift 2 ;;
        --keep)     KEEP=1; shift ;;
        --real-run) REAL_RUN=1; shift ;;
        --replicate) REPLICATE=1; shift ;;
        --qc-full)   QC_FULL=1; shift ;;
        --motif)     MOTIF=1; shift ;;
        --diffbind)  DIFFBIND=1; shift ;;
        --gates)     GATES=1; shift ;;
        --spike-in)  SPIKE=1; shift ;;
        --seacr)     SEACR=1; shift ;;
        -h|--help)  usage; exit 0 ;;
        *) echo "[ERROR] Unknown argument: $1 (see --help for usage)" >&2; exit 1 ;;
    esac
done
MODE="dry-run"
[[ "$REAL_RUN" == 1 ]] && MODE="real-run"
SCENARIO_ARGS=()
[[ "$REPLICATE" == 1 ]] && SCENARIO_ARGS+=(--replicate)
[[ "$QC_FULL" == 1 ]] && SCENARIO_ARGS+=(--qc-full)
[[ "$MOTIF" == 1 ]] && SCENARIO_ARGS+=(--motif)
[[ "$DIFFBIND" == 1 ]] && SCENARIO_ARGS+=(--diffbind)
[[ "$GATES" == 1 ]] && SCENARIO_ARGS+=(--gates)
[[ "$SPIKE" == 1 ]] && SCENARIO_ARGS+=(--spike-in)
[[ "$SEACR" == 1 ]] && SCENARIO_ARGS+=(--seacr)
if [[ "$REPLICATE" == 1 && "$SEACR" == 1 ]]; then
    echo "[ERROR] --replicate and --seacr are mutually exclusive (peak.caller=seacr" \
         "rejects peak.replicate.enabled at parse time)" >&2
    exit 1
fi

echo "[test] 1/5 Checking dependencies (mode=$MODE, reads=$READS, scenario_args=${SCENARIO_ARGS[*]:-none})"
command -v snakemake >/dev/null || { echo "[ERROR] snakemake not found" >&2; exit 1; }
command -v python3   >/dev/null || { echo "[ERROR] python3 not found" >&2; exit 1; }
# run.sh needs pyyaml to resolve software.yaml at runtime; fail early with a clear message
python3 -c 'import yaml' >/dev/null 2>&1 || { echo "[ERROR] python3 is missing pyyaml" >&2; exit 1; }

echo "[test] 2/5 Generating synthetic test data (reads=$READS, fixed seed)"
python3 "$TESTS_DIR/make_testdata.py" --outdir "$DATA_DIR" --reads "$READS" "${SCENARIO_ARGS[@]}"

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
    # ---------- real-run output existence assertions ----------
    # default: 5 samples; --replicate adds the 3-sample broad group g3
    N_EXPECTED=5
    [[ "$REPLICATE" == 1 ]] && N_EXPECTED=8
    shopt -s nullglob
    bams=("$WORK_DIR"/results/3.align/bowtie2/*_sorted.bam)
    shopt -u nullglob
    if [[ ${#bams[@]} -eq $N_EXPECTED ]]; then
        echo "  PASS  results/3.align/bowtie2/*_sorted.bam count ${#bams[@]} (expected $N_EXPECTED)"
    else
        echo "  FAIL  results/3.align/bowtie2/*_sorted.bam count ${#bams[@]} (expected $N_EXPECTED)"; FAIL=1
    fi
    EXPECTED=(
        "results/4.peak/g1_peaks.narrowPeak"
        "results/4.peak/g2_peaks.narrowPeak"
        "results/5.QC/frip/FRiP_summary.tsv"
        "results/2.cleandata/fastqc/multiqc/multiqc_report.html"
        "results/5.QC/software_versions.yaml"
    )
    if [[ "$REPLICATE" == 1 ]]; then
        EXPECTED+=(
            "results/4.peak/g3_peaks.broadPeak"
            "results/4.peak/g1_IDR_peaks.narrowPeak"
            "results/4.peak/g2_IDR_peaks.narrowPeak"
            "results/4.peak/g3_consensus_peaks.broadPeak"
            "results/5.QC/replicate_peaks/Replicate_summary.tsv"
        )
    fi
    if [[ "$QC_FULL" == 1 ]]; then
        EXPECTED+=(
            "results/5.QC/tss/TSSE_summary.tsv"
            "results/5.QC/organelle/Organelle_summary.tsv"
            "results/5.QC/blacklist/blacklist_summary.tsv"
            "results/4.peak/blacklist_filtered/g1_peaks.narrowPeak"
        )
    fi
    if [[ "$MOTIF" == 1 ]]; then
        EXPECTED+=("results/6.motif/g1")
    fi
    if [[ "$DIFFBIND" == 1 ]]; then
        EXPECTED+=(
            "results/6.diffbind/g1__vs__g4/samplesheet.tsv"
            "results/6.diffbind/g1__vs__g4/DB_results.tsv"
        )
    fi
    if [[ "$GATES" == 1 ]]; then
        EXPECTED+=(
            "results/5.QC/gates/gate_summary.tsv"
            "results/5.QC/gates/gate_summary_mqc.tsv"
        )
    fi
    if [[ "$SPIKE" == 1 ]]; then
        EXPECTED+=(
            "results/3.align/spike_in/chip_treat_rep1_sorted.bam"
            "results/5.QC/spike_in/Spikein_summary.tsv"
            "results/5.QC/spike_in/Spikein_summary_mqc.tsv"
        )
    fi
    if [[ "$SEACR" == 1 ]]; then
        EXPECTED+=(
            "results/4.peak/seacr/g1_treat.bg"
            "results/4.peak/seacr/g1_control.bg"
            "results/4.peak/seacr/g2_treat.bg"
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
    # ---------- dry-run DAG assertions: exit code 0 (checked above) and key rule names present ----------
    SNAKE_LOG="$WORK_DIR/snakemake.logs.txt"   # run.sh default log (stdout is tee'd to disk)
    # the pooled caller rules differ per scenario: MACS2 (default) or SEACR
    DAG_RULES=(bowtie2_mapping)
    if [[ "$SEACR" == 1 ]]; then
        DAG_RULES+=(seacr_bedgraph_treat seacr_bedgraph_control seacr_callpeak seacr_bigwig)
    else
        DAG_RULES+=(callpeak_narrow callpeak_atac)
    fi
    if [[ "$REPLICATE" == 1 ]]; then
        DAG_RULES+=(callpeak_narrow_replicate callpeak_broad_replicate idr_pair idr_final broad_consensus replicate_summary)
    fi
    if [[ "$QC_FULL" == 1 ]]; then
        DAG_RULES+=(tss_matrix tss_summary organelle_idxstats organelle_summary blacklist_filter blacklist_summary)
    fi
    if [[ "$MOTIF" == 1 ]]; then
        DAG_RULES+=(motif_enrichment)
    fi
    if [[ "$DIFFBIND" == 1 ]]; then
        DAG_RULES+=(diffbind_sheet diffbind_report)
    fi
    if [[ "$GATES" == 1 ]]; then
        DAG_RULES+=(gates_flagstat qc_gates)
    fi
    if [[ "$SPIKE" == 1 ]]; then
        DAG_RULES+=(spike_bowtie2_index spike_align spike_summary)
    fi
    for rule in "${DAG_RULES[@]}"; do
        if grep -q "$rule" "$CAPTURE_LOG" "$SNAKE_LOG" 2>/dev/null; then
            echo "  PASS  DAG contains rule $rule"
        else
            echo "  FAIL  DAG does not contain rule $rule"; FAIL=1
        fi
    done
    # Exact job-header patterns: the bare "frip" token would also match FRiP
    # file paths in other rules' input listings, and "bigwig" is a substring of
    # seacr_bigwig — match the "rule X:"/"localrule X:" headers instead.
    PRESENT_HEADER_PATTERNS=()
    ABSENT_HEADER_PATTERNS=()
    if [[ "$SEACR" == 1 ]]; then
        PRESENT_HEADER_PATTERNS=("rule frip:")
        # the pooled MACS2 caller rules (and the MACS2 bigwig rule, replaced by
        # seacr_bigwig) must be absent from the SEACR-mode DAG
        ABSENT_HEADER_PATTERNS=("rule bigwig:")
    fi
    for pattern in "${PRESENT_HEADER_PATTERNS[@]}"; do
        if grep -q "$pattern" "$CAPTURE_LOG" "$SNAKE_LOG" 2>/dev/null; then
            echo "  PASS  DAG contains a $pattern job"
        else
            echo "  FAIL  DAG does not contain a $pattern job"; FAIL=1
        fi
    done
    ABSENT_RULES=()
    if [[ "$SEACR" == 1 ]]; then
        ABSENT_RULES+=(callpeak_narrow callpeak_broad callpeak_atac)
    fi
    for rule in "${ABSENT_RULES[@]}"; do
        if grep -q "$rule" "$CAPTURE_LOG" "$SNAKE_LOG" 2>/dev/null; then
            echo "  FAIL  DAG must not contain rule $rule under the SEACR caller"; FAIL=1
        else
            echo "  PASS  DAG does not contain rule $rule"
        fi
    done
    for pattern in "${ABSENT_HEADER_PATTERNS[@]}"; do
        if grep -q "$pattern" "$CAPTURE_LOG" "$SNAKE_LOG" 2>/dev/null; then
            echo "  FAIL  DAG must not contain a $pattern job under the SEACR caller"; FAIL=1
        else
            echo "  PASS  DAG does not contain a $pattern job"
        fi
    done
    # Job count (informational baseline; printed for the CHANGELOG/AGENTS
    # record). Snakemake 7 dry-runs print one rule/localrule block per planned
    # job (printshellcmds is pinned on in every profile); counting both keeps
    # the historical baseline convention (chip default = 47, including `all`).
    N_JOBS=$(grep -cE "^(local)?rule " "$SNAKE_LOG" 2>/dev/null || true)
    echo "  INFO  dry-run planned jobs: ${N_JOBS:-0} (scenario: ${SCENARIO_ARGS[*]:-default})"
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
