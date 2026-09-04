#!/bin/bash
#########################################################################
# 组间 DEG 集合比较：生成两两组合任务 → 并行 GO 富集比较 → 交集基因注释
# Usage: DEGgroupCompare.sh <results_dir> <control_group> <annotation_tsv> <sample_info.csv> <threads>
# 并行执行用 xargs -P 替代原 ParaFly（不再依赖个人软件路径）
#########################################################################

set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

RD=$1
control=$2
anno=$3
samplefile=$4
cpu=$5

mkdir -p "$RD/6.DEGcompare/combination"
python "$SCRIPT_DIR/getGroups.py" "$samplefile" "$RD/5.DEG/DEGs" "$control" "$RD/6.DEGcompare"

## 生成任务列表
rm -f "$RD/6.DEGcompare/jobs"
ls "$RD"/6.DEGcompare/combination/* | while read line;
do
    echo "Rscript $SCRIPT_DIR/run_deg_compare.R $line" >> "$RD/6.DEGcompare/jobs"
done

## 并行执行（每行一条命令，-P 控制并发数；-r 处理组 <2 时任务列表为空的情况）
xargs -r -P "$cpu" -d '\n' -I CMD sh -c 'CMD' < "$RD/6.DEGcompare/jobs"

## 交集基因功能注释
mkdir -p "$RD/6.DEGcompare/intersetionGene_anno"
ls "$RD"/6.DEGcompare/*_intersetion.xls | while read id;
do
    prefix=$(basename "$id")
    cat "$id" | fgrep ∩ | cut -f5 | tr " " "\n" | sed 's/,//g' > tmp.gene
    cat "$anno" | fgrep -w -f tmp.gene > "$RD/6.DEGcompare/intersetionGene_anno/${prefix%_*}_intersetionGene_anno.xls" && rm tmp.gene
done

echo "$(date) : All Done!" > "$RD/6.DEGcompare/flag.log"
