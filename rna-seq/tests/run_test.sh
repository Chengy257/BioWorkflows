#!/bin/bash
#########################################################################
# 一键回归测试（对应路线图 3.1）
#   生成微型测试数据 → dry-run → 端到端运行 → 输出断言 → （可选）DAG 再生成
#
# 用法:
#   bash tests/run_test.sh [--pipeline deg|upstream|as] [--reads 50000] [--keep]
#     --pipeline  测试模式，默认 deg（覆盖 align+quant+deg 全部规则）
#     --reads     每样本 PE 读对数，默认 50000
#     --keep      保留 tests/data 与 tests/work（默认结束时清理）
# 依赖: 已配置好的统一 RNA-seq 软件环境（含 snakemake/分析工具/R 包）与 python3；SGE/SLURM 非必需。
# 说明: lncrna 管线依赖外部工具（CPC2/CNCI/pfam_scan.pl），不在默认测试范围。
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
        *) echo "[ERROR] 未知参数: $1" >&2; exit 1 ;;
    esac
done
[[ "$PIPELINE" == "lncrna" ]] && { echo "[ERROR] lncrna 依赖外部工具，不在自动测试范围（upstream/deg/as 可选）" >&2; exit 1; }

echo "[test] 1/5 检查依赖"
command -v snakemake >/dev/null || { echo "[ERROR] 未找到 snakemake" >&2; exit 1; }
command -v python3   >/dev/null || { echo "[ERROR] 未找到 python3" >&2; exit 1; }

echo "[test] 2/5 生成微型测试数据（reads=$READS, seed 固定）"
python3 "$TESTS_DIR/make_testdata.py" --outdir "$DATA_DIR" --reads "$READS"

echo "[test] 3/5 组装测试项目 tests/work"
rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR/1.rawdata"
cp "$DATA_DIR"/rawdata/*.gz "$WORK_DIR/1.rawdata/"
cp "$DATA_DIR"/reference/genome.fa "$DATA_DIR"/reference/genes.gtf \
   "$DATA_DIR"/reference/genes.bed "$DATA_DIR"/reference/annotation_full.tsv "$WORK_DIR/"
cp "$DATA_DIR"/samples.csv "$WORK_DIR/samples.csv"

## 测试专用配置：小基因组需调低 STAR SA 索引参数；CI 使用 hsa 的标准 OrgDb 环境。
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

echo "[test] 5/5 端到端运行（复用当前统一软件环境）"
bash "$REPO_DIR/run.sh" -p "$PIPELINE" -P "$WORK_DIR" -c "$WORK_DIR/config.yaml" \
    --software "$WORK_DIR/software.yaml" -j 2

echo "[test] 输出断言"
python3 "$TESTS_DIR/check_outputs.py" --work "$WORK_DIR" --pipeline "$PIPELINE" --reads "$READS"

## DAG 再生成（P2-4，尽力而为）
if command -v dot >/dev/null 2>&1; then
    (cd "$WORK_DIR" && RNASEQ_PIPELINE="$PIPELINE" RNASEQ_CONFIG="$WORK_DIR/config.yaml" \
        snakemake -s "$REPO_DIR/workflow/Snakefile" --dag | dot -Tsvg -o "$REPO_DIR/docs/dag_$PIPELINE.svg") \
        && echo "[test] DAG 图已更新: docs/dag_$PIPELINE.svg" \
        || echo "[test] [WARN] DAG 生成失败（跳过，不影响测试结果）"
else
    echo "[test] [WARN] 未安装 graphviz，跳过 DAG 再生成"
fi

echo "[test] 全部通过 ✔ （pipeline=$PIPELINE reads=$READS）"
if [[ "$KEEP" != "1" ]]; then
    rm -rf "$WORK_DIR" "$DATA_DIR"
    echo "[test] 已清理 tests/work 与 tests/data（--keep 可保留）"
fi
