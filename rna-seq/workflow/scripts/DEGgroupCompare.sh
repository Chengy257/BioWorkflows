#!/bin/bash
#########################################################################
# 组间 DEG 集合比较：生成两两组合任务 → 并行 GO 富集比较 → 交集基因注释
# Usage: DEGgroupCompare.sh <results_dir> <control_group> <annotation_tsv> <sample_info.csv> <threads> [species]
#   species        osa | hsa（默认 osa），转发给 run_deg_compare.R
# 并行执行用 xargs -P 替代原 ParaFly（不再依赖个人软件路径）
#########################################################################

set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RSCRIPT="${RNASEQ_RSCRIPT:-Rscript}"
PYTHON="${RNASEQ_PYTHON:-python3}"

RD=$1
control=$2
anno=$3
samplefile=$4
cpu=$5
spe=${6:-osa}

## 交集注释用的临时文件（mktemp + trap，避免残留在运行目录）
tmp_gene=$(mktemp)
trap 'rm -f "$tmp_gene"' EXIT

mkdir -p "$RD/6.DEGcompare/combination"
"$PYTHON" "$SCRIPT_DIR/getGroups.py" "$samplefile" "$RD/5.DEG/DEGs" "$control" "$RD/6.DEGcompare"

## 生成任务列表（species / orgdb 贯通到每个任务）
rm -f "$RD/6.DEGcompare/jobs"
ls "$RD"/6.DEGcompare/combination/* | while read line;
do
    printf '%q %q %q %q\n' "$RSCRIPT" "$SCRIPT_DIR/run_deg_compare.R" "$line" "$spe" >> "$RD/6.DEGcompare/jobs"
done

## 并行执行（每行一条命令，-P 控制并发数；-r 处理组 <2 时任务列表为空的情况）
## 单个任务失败不中断整体（与 enrich 的容错语义一致），但显式告警
xargs -r -P "$cpu" -d '\n' -I CMD sh -c 'CMD' < "$RD/6.DEGcompare/jobs" \
    || echo "[DEGgroupCompare] [WARN] Some comparison jobs failed; results may be incomplete." >&2

## 交集基因功能注释
mkdir -p "$RD/6.DEGcompare/intersetionGene_anno"
ls "$RD"/6.DEGcompare/*_intersetion.xls | while read id;
do
    prefix=$(basename "$id")
    cat "$id" | fgrep ∩ | cut -f5 | tr " " "\n" | sed 's/,//g' > "$tmp_gene"
    cat "$anno" | fgrep -w -f "$tmp_gene" > "$RD/6.DEGcompare/intersetionGene_anno/${prefix%_*}_intersetionGene_anno.xls"
done

echo "$(date) : All Done!" > "$RD/6.DEGcompare/flag.log"
