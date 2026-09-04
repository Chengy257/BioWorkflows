#!/bin/bash
#########################################################################
# rna-seq 工作流统一启动脚本（仓库根目录）
# 用法: bash run.sh <pipeline> <project_dir> [config.yaml] [jobs]
#   pipeline     upstream | deg | as | lncrna（大小写不敏感）
#   project_dir  分析项目目录（含 1.rawdata/ 与样本表），不存在则自动创建
#   config       可选；项目配置（默认取项目目录内 config.yaml，仍无则用仓库
#                config/config.yaml 默认配置）
#   jobs         可选；最大并行任务数，默认 10
# pipeline 与配置经环境变量 RNASEQ_PIPELINE / RNASEQ_CONFIG 注入
# workflow/Snakefile；调度器自动选择（有 qsub 用 SGE，否则本地运行），
# 可用环境变量 RUN_PROFILE=sge|default 强制指定。要求 snakemake >= 7。
#########################################################################
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_DIR="$REPO_DIR/workflow"

pipeline=${1:-upstream}
PROJ_DIR=${2:?用法: bash run.sh <upstream|deg|as|lncrna> <project_dir> [config.yaml] [jobs]}
config=${3:-}
jobs=${4:-10}

pipeline=$(echo "$pipeline" | tr '[:upper:]' '[:lower:]')
case "$pipeline" in
    upstream|deg|as|lncrna) ;;
    *) echo "[ERROR] 未知 pipeline: $pipeline（可选 upstream|deg|as|lncrna）" >&2; exit 1 ;;
esac

mkdir -p "$PROJ_DIR"
PROJ_DIR="$(cd "$PROJ_DIR" && pwd)"
cd "$PROJ_DIR"

## 配置：命令行参数 > 项目目录 config.yaml > 仓库默认配置
if [[ -z "$config" && -f "$PROJ_DIR/config.yaml" ]]; then
    config="$PROJ_DIR/config.yaml"
fi
if [[ -n "$config" ]]; then
    config="$(cd "$(dirname "$config")" && pwd)/$(basename "$config")"
fi

## 注入 Snakefile 所需环境变量（解析期确定，不受配置合并顺序影响）
export RNASEQ_PIPELINE="$pipeline"
if [[ -n "$config" ]]; then
    export RNASEQ_CONFIG="$config"
    cfg_for_parse="$config"
else
    cfg_for_parse="$REPO_DIR/config/config.yaml"
fi

## 样本表校验（尽力而为：能从配置解析出样本表时才执行）
sample_list=$(awk -F': *' '/^SampleListFile:/ {gsub(/"/, "", $2); print $2}' "$cfg_for_parse" | tail -1)
control_group=$(awk -F': *' '/^control_group:/ {gsub(/"/, "", $2); print $2}' "$cfg_for_parse" | tail -1)
control_group=${control_group:-control}
if [[ -n "$sample_list" && -f "$sample_list" ]] && command -v python3 >/dev/null 2>&1; then
    echo "[run.sh] 校验样本表 $sample_list ..."
    python3 "$WORKFLOW_DIR/scripts/validate_samples.py" "$sample_list" "$control_group"
else
    echo "[run.sh] 跳过样本表校验（未找到 $sample_list 或无 python3）"
fi

## 调度 profile 选择
if [[ -n "${RUN_PROFILE:-}" ]]; then
    profile="$RUN_PROFILE"
elif command -v qsub >/dev/null 2>&1; then
    profile=sge
else
    profile=default
fi
echo "[run.sh] pipeline=$pipeline profile=$profile jobs=$jobs cfg=${config:-<仓库默认>}"

snakemake -s "$WORKFLOW_DIR/Snakefile" \
    --profile "$WORKFLOW_DIR/profile/$profile" \
    -j "$jobs" \
    2>&1 | tee snakemake.logs.txt
