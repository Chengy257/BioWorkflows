#!/bin/bash
#########################################################################
# 一键回归测试（对应实施计划 Phase 4 Task 4.3，骨架移植自 rna-seq 同名脚本）
#   生成合成测试数据 → 组装测试工作目录 → dry-run（默认）或 --real-run 端到端
#   → 断言 → 清理
#
# 用法:
#   bash tests/run_test.sh [--reads N] [--keep] [--real-run] [--help]
#     --reads     每样本 PE 读对数，默认 50000（CI 回归传 2000）
#     --keep      保留 tests/data 与 tests/work（默认结束时清理）
#     --real-run  端到端实跑并断言产物存在（默认仅 dry-run 验证 DAG 完整性；
#                 实跑需含 bowtie2/fastqc/trim_galore/macs2/R 等的完整分析
#                 环境，供服务器实跑验证）
# 依赖: dry-run 仅需 snakemake + python3(+pyyaml)；--real-run 需完整分析环境。
#########################################################################
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TESTS_DIR="$REPO_DIR/tests"
DATA_DIR="$TESTS_DIR/data"
WORK_DIR="$TESTS_DIR/work"

usage() {
    cat <<'EOF'
一键回归测试：合成数据 → 组装工作目录 → dry-run（默认）/ --real-run 端到端 → 断言

用法:
  bash tests/run_test.sh [--reads N] [--keep] [--real-run] [--help]
    --reads     每样本 PE 读对数，默认 50000（CI 回归传 2000）
    --keep      保留 tests/data 与 tests/work（默认结束时清理）
    --real-run  端到端实跑并断言产物存在（默认仅 dry-run 验证 DAG 完整性；
                实跑需完整分析环境，供服务器实跑验证）
    -h, --help  显示本帮助

依赖:
  dry-run 仅需 snakemake + python3(+pyyaml)；--real-run 需完整分析环境
  （bowtie2/fastqc/trim_galore/macs2/deeptools/R 等，见 workflow/environment.yaml）。
EOF
}

# ---------------------------------------------------------------------
# 参数解析：默认 dry-run；--real-run 打开端到端实跑分支
# ---------------------------------------------------------------------
READS=50000
KEEP=0
REAL_RUN=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --reads)
            [[ $# -ge 2 ]] || { echo "[ERROR] --reads 需要一个值" >&2; exit 1; }
            [[ "$2" =~ ^[0-9]+$ ]] || { echo "[ERROR] --reads 必须是正整数: $2" >&2; exit 1; }
            READS="$2"; shift 2 ;;
        --keep)     KEEP=1; shift ;;
        --real-run) REAL_RUN=1; shift ;;
        -h|--help)  usage; exit 0 ;;
        *) echo "[ERROR] 未知参数: $1（--help 查看用法）" >&2; exit 1 ;;
    esac
done
MODE="dry-run"
[[ "$REAL_RUN" == 1 ]] && MODE="real-run"

echo "[test] 1/5 检查依赖（mode=$MODE, reads=$READS）"
command -v snakemake >/dev/null || { echo "[ERROR] 未找到 snakemake" >&2; exit 1; }
command -v python3   >/dev/null || { echo "[ERROR] 未找到 python3" >&2; exit 1; }
# run.sh 运行时解析 software.yaml 需要 pyyaml，提前给出友好报错
python3 -c 'import yaml' >/dev/null 2>&1 || { echo "[ERROR] python3 缺少 pyyaml" >&2; exit 1; }

echo "[test] 2/5 生成合成测试数据（reads=$READS, seed 固定）"
python3 "$TESTS_DIR/make_testdata.py" --outdir "$DATA_DIR" --reads "$READS"

echo "[test] 3/5 组装测试项目 tests/work"
rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR"
# 1.rawdata 优先软链接入（省时省空间）；软链不可用（如 Windows 默认配置）时退回复制
if ln -s "$DATA_DIR/1.rawdata" "$WORK_DIR/1.rawdata" 2>/dev/null; then
    echo "[test]   1.rawdata: 符号链接接入"
else
    mkdir -p "$WORK_DIR/1.rawdata"
    cp "$DATA_DIR"/1.rawdata/*.fq.gz "$WORK_DIR/1.rawdata/"
    echo "[test]   1.rawdata: 复制接入（软链不可用）"
fi
cp -r "$DATA_DIR/ref" "$WORK_DIR/ref"
cp "$DATA_DIR/samples.csv" "$DATA_DIR/config.yaml" "$WORK_DIR/"

echo "[test] 4/5 运行流程（$MODE）"
RUN_ARGS=(-P "$WORK_DIR" -c "$WORK_DIR/config.yaml" -j 4)
# dry-run（默认）加 -n：只构建 DAG 不执行；--real-run 去掉 -n 端到端实跑
[[ "$REAL_RUN" == 1 ]] || RUN_ARGS+=(-n)
CAPTURE_LOG="$WORK_DIR/run_test.capture.log"
set +e
bash "$REPO_DIR/run.sh" "${RUN_ARGS[@]}" >"$CAPTURE_LOG" 2>&1
STATUS=$?
set -e
if [[ "$STATUS" != 0 ]]; then
    echo "[ERROR] run.sh 退出码 $STATUS（输出尾部如下）" >&2
    tail -40 "$CAPTURE_LOG" >&2 || true
    exit 1
fi

echo "[test] 5/5 断言"
FAIL=0
if [[ "$REAL_RUN" == 1 ]]; then
    # ---------- 实跑产物存在性断言（样本数 = 5：chip 组 3 + atac 组 2） ----------
    shopt -s nullglob
    bams=("$WORK_DIR"/results/3.align/bowtie2/*_sorted.bam)
    shopt -u nullglob
    if [[ ${#bams[@]} -eq 5 ]]; then
        echo "  PASS  3.align/bowtie2/*_sorted.bam 共 ${#bams[@]} 个（期望 5）"
    else
        echo "  FAIL  3.align/bowtie2/*_sorted.bam 共 ${#bams[@]} 个（期望 5）"; FAIL=1
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
            echo "  FAIL  $rel（不存在或为空）"; FAIL=1
        fi
    done
else
    # ---------- dry-run DAG 断言：退出码 0（上文已验）且输出含三个关键规则名 ----------
    SNAKE_LOG="$WORK_DIR/snakemake.logs.txt"   # run.sh 默认日志（stdout 经 tee 落盘）
    for rule in bowtie2_mapping callpeak_narrow callpeak_atac; do
        if grep -q "$rule" "$CAPTURE_LOG" "$SNAKE_LOG" 2>/dev/null; then
            echo "  PASS  DAG 含规则 $rule"
        else
            echo "  FAIL  DAG 未含规则 $rule"; FAIL=1
        fi
    done
fi
if [[ "$FAIL" != "0" ]]; then
    echo "[ERROR] 断言失败，保留现场供排查: $WORK_DIR" >&2
    exit 1
fi

## DAG 再生成（尽力而为：本机装了 graphviz 才执行，失败不影响测试结果）
if command -v dot >/dev/null 2>&1; then
    (cd "$WORK_DIR" && CHIP_CONFIG="$WORK_DIR/config.yaml" \
        snakemake -s "$REPO_DIR/workflow/Snakefile" --configfile "$WORK_DIR/config.yaml" --dag \
        | dot -Tsvg -o "$REPO_DIR/docs/dag_test.svg") \
        && echo "[test] DAG 图已更新: docs/dag_test.svg" \
        || echo "[test] [WARN] DAG 生成失败（跳过，不影响测试结果）"
else
    echo "[test] [WARN] 未安装 graphviz，跳过 DAG 再生成"
fi

echo "[test] 全部通过 ✔ （mode=$MODE reads=$READS）"
if [[ "$KEEP" != "1" ]]; then
    rm -rf "$WORK_DIR" "$DATA_DIR"
    echo "[test] 已清理 tests/work 与 tests/data（--keep 可保留）"
fi
