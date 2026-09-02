#!/bin/bash
#########################################################################
# 组间 DEG 集合比较：生成两两组合任务 → 并行 GO 富集比较 → 交集基因注释
# Usage: DEGgroupCompare.sh <control_group> <annotation_tsv> <sample_info.csv> <threads>
# 并行执行用 xargs -P 替代原 ParaFly（不再依赖个人软件路径）
#########################################################################

set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

control=$1
anno=$2
samplefile=$3
cpu=$4

mkdir -p 6.DEGcompare/combination
python "$SCRIPT_DIR/getGroups.py" "$samplefile" 5.DEG/DEGs/ "$control"

## 生成任务列表
rm -f 6.DEGcompare/jobs
ls 6.DEGcompare/combination/* | while read line;
do
    echo "Rscript $SCRIPT_DIR/multi_enrich_GOKEGG.R $line" >> 6.DEGcompare/jobs
done

## 并行执行（每行一条命令，-P 控制并发数）
xargs -P "$cpu" -d '\n' -I CMD sh -c 'CMD' < 6.DEGcompare/jobs

## 交集基因功能注释
mkdir -p 6.DEGcompare/intersetionGene_anno
ls 6.DEGcompare/*_intersetion.xls | while read id;
do
    prefix=$(basename "$id")
    cat "$id" | fgrep ∩ | cut -f5 | tr " " "\n" | sed 's/,//g' > tmp.gene
    cat "$anno" | fgrep -w -f tmp.gene > "6.DEGcompare/intersetionGene_anno/${prefix%_*}_intersetionGene_anno.xls" && rm tmp.gene
done

echo "$(date) : All Done!" > 6.DEGcompare/flag.log
