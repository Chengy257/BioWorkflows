#!/bin/bash
#########################################################################
# rna-seq 工作流统一启动脚本
# 用法: bash run.sh <pipeline> <project_dir> [user_config.yaml] [jobs]
#   pipeline     upstream | deg | as | lncrna
#   project_dir  分析项目目录（含 1.rawdata/ 与样本表），不存在则自动创建
#   user_config  可选；覆盖默认配置（复制 config_user_defined.yaml 模板修改）
#   jobs         可选；最大并行任务数，默认 10
# 调度器自动选择：有 qsub 用 SGE（profile/sge），否则本地运行（profile/default）；
# 可用环境变量 RUN_PROFILE=sge|default 强制指定。要求 snakemake >= 7。
#########################################################################
set -euo pipefail

WORKFLOW_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

pipeline=${1:-upstream}
PROJ_DIR=${2:?用法: bash run.sh <upstream|deg|as|lncrna> <project_dir> [config.yaml] [jobs]}
user_config=${3:-}
jobs=${4:-10}

case "$pipeline" in
    upstream) smk=RNA-seq_upstream_only.smk ;;
    deg)      smk=RNA-seq_up_DEG.smk ;;
    as)       smk=RNA-seq_up_AS.smk ;;
    lncrna)   smk=RNA-seq_up_lncRNA.smk ;;
    *) echo "[ERROR] 未知 pipeline: $pipeline（可选 upstream|deg|as|lncrna）" >&2; exit 1 ;;
esac

mkdir -p "$PROJ_DIR"
cd "$PROJ_DIR"

## 配置优先级：命令行 config > 项目目录 config_user_defined.yaml > 工作流默认配置
if [[ -n "$user_config" ]]; then
    cfg="$user_config"
elif [[ -f config_user_defined.yaml ]]; then
    cfg=config_user_defined.yaml
else
    cfg=""
fi
cfg_for_parse=${cfg:-$WORKFLOW_DIR/config_basic_defaulted.yaml}

## 样本表校验（能从配置解析出样本表时才执行）
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
echo "[run.sh] pipeline=$pipeline profile=$profile jobs=$jobs cfg=${cfg:-<默认>}"

cfg_arg=""
[[ -n "$cfg" ]] && cfg_arg="--configfile $cfg"

snakemake -s "$WORKFLOW_DIR/$smk" \
    --profile "$WORKFLOW_DIR/profile/$profile" \
    $cfg_arg \
    -j "$jobs" \
    2>&1 | tee snakemake.logs.txt
