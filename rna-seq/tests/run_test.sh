#!/bin/bash
#########################################################################
# One-shot regression test
#   Generate miniature test data -> dry-run -> end-to-end run -> output assertions -> (optional) DAG regeneration
#
# Usage:
#   bash tests/run_test.sh [--pipeline deg|upstream|as] [--reads 50000] [--keep]
#     --pipeline  pipeline to test, default deg (covers all align+quant+deg rules)
#     --reads     PE read pairs per sample, default 50000
#     --keep      keep tests/data and tests/work (cleaned up at the end by default)
# Requires: an already configured unified RNA-seq software environment (snakemake/analysis tools/R packages) and python3;
#           SGE/SLURM is not required.
# Note: the lncrna pipeline depends on external tools (CPC2/CNCI/pfam_scan.pl) and is not covered by the default test.
#########################################################################
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TESTS_DIR="$REPO_DIR/tests"
DATA_DIR="$TESTS_DIR/data"
WORK_DIR="$TESTS_DIR/work"

PIPELINE=deg
READS=50000
KEEP=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --pipeline) PIPELINE="$2"; shift 2 ;;
        --reads)    READS="$2"; shift 2 ;;
        --keep)     KEEP=1; shift ;;
        -h|--help)  grep '^#' "$0" | head -20; exit 0 ;;
        *) echo "[ERROR] unknown argument: $1" >&2; exit 1 ;;
    esac
done
[[ "$PIPELINE" == "lncrna" ]] && { echo "[ERROR] the lncrna pipeline depends on external tools and is not covered by automated tests (choose upstream/deg/as)" >&2; exit 1; }

echo "[test] 1/5 checking dependencies"
command -v snakemake >/dev/null || { echo "[ERROR] snakemake not found" >&2; exit 1; }
command -v python3   >/dev/null || { echo "[ERROR] python3 not found" >&2; exit 1; }

echo "[test] 2/5 generating miniature test data (reads=$READS, fixed seed)"
python3 "$TESTS_DIR/make_testdata.py" --outdir "$DATA_DIR" --reads "$READS"

echo "[test] 3/5 assembling the test project tests/work"
rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR/1.rawdata"
cp "$DATA_DIR"/rawdata/*.gz "$WORK_DIR/1.rawdata/"
cp "$DATA_DIR"/reference/genome.fa "$DATA_DIR"/reference/genes.gtf \
   "$DATA_DIR"/reference/genes.bed "$DATA_DIR"/reference/annotation_full.tsv "$WORK_DIR/"
cp "$DATA_DIR"/samples.csv "$WORK_DIR/samples.csv"

## Test-specific configuration: the small genome requires a lowered STAR SA index parameter; CI uses hsa's standard OrgDb environment.
cat > "$WORK_DIR/config.yaml" <<EOF
pipeline: "$PIPELINE"
results_dir: "results"
SampleListFile: "samples.csv"
control_group: "control"
threads: 4
FoldChange: "2"
padj: "0.05"
pca_ntop: 20000
batch_correction: "F"
species: "hsa"
bed: "genes.bed"
genome: "genome.fa"
gtf: "genes.gtf"
annotation_tsv: "annotation_full.tsv"
star_extra_args: "--genomeSAindexNbases 7"
lncrna:
  gtf: ""
  gtf_PcGs: "genes.gtf"
  threads: 4
EOF
cat > "$WORK_DIR/software.yaml" <<EOF
# CI/regression tests run inside one pre-provisioned environment.
environment:
  type: system
  strict: true
r:
  rscript: "Rscript"
  version: ""
  version_check: off
  lib_paths: []
  lib_mode: prepend
tools: {}
paths: {}
databases: {}
EOF

echo "[test] 4/5 dry-run"
bash "$REPO_DIR/run.sh" -p "$PIPELINE" -P "$WORK_DIR" -c "$WORK_DIR/config.yaml" \
    --software "$WORK_DIR/software.yaml" -j 2 --dry-run --quiet

echo "[test] 5/5 end-to-end run (reusing the current unified software environment)"
bash "$REPO_DIR/run.sh" -p "$PIPELINE" -P "$WORK_DIR" -c "$WORK_DIR/config.yaml" \
    --software "$WORK_DIR/software.yaml" -j 2

echo "[test] output assertions"
python3 "$TESTS_DIR/check_outputs.py" --work "$WORK_DIR" --pipeline "$PIPELINE" --reads "$READS"

## DAG regeneration (P2-4, best effort)
if command -v dot >/dev/null 2>&1; then
    (cd "$WORK_DIR" && RNASEQ_PIPELINE="$PIPELINE" RNASEQ_CONFIG="$WORK_DIR/config.yaml" \
        snakemake -s "$REPO_DIR/workflow/Snakefile" --dag | dot -Tsvg -o "$REPO_DIR/docs/dag_$PIPELINE.svg") \
        && echo "[test] DAG figure updated: docs/dag_$PIPELINE.svg" \
        || echo "[test] [WARN] DAG generation failed (skipped; does not affect the test result)"
else
    echo "[test] [WARN] graphviz not installed; skipping DAG regeneration"
fi

echo "[test] all checks passed (pipeline=$PIPELINE reads=$READS)"
if [[ "$KEEP" != "1" ]]; then
    rm -rf "$WORK_DIR" "$DATA_DIR"
    echo "[test] cleaned tests/work and tests/data (use --keep to retain them)"
fi
