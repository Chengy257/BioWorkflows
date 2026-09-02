#!/bin/bash
#########################################################################
# File Name: /home/chengyu/workflows/snakemake/rna-seq-workflow/scripts/DEGgroupCompare.sh
# Author: ChengYu
# Description: 
# Created Time: Fri Jul 14 20:41:50 2023
#########################################################################

# /usr/bin/Rscript  /home/chengyu/workflows/snakemake/rna-seq-workflow/scripts/multi_enrich_GOKEGG.R 

# DIR=${1}

mkdir -p 6.DEGcompare/combination/
python /home/chengyu/workflows/snakemake/rna-seq-workflow/scripts/getGroups.py sample_info.csv 5.DEG/DEGs/
# cd ${DIR}
rm -f 6.DEGcompare/jobs
ls 6.DEGcompare/combination/* |while read line ;
do
    echo /usr/bin/Rscript /home/chengyu/workflows/snakemake/rna-seq-workflow/scripts/multi_enrich_GOKEGG.R ${line} >> 6.DEGcompare/jobs
done
/home/chengyu/soft/miniconda3/bin/ParaFly -c 6.DEGcompare/jobs -CPU 10

anno=/home/chengyu/references/osa/rap/20230315/annotation_full.tsv
mkdir -p 6.DEGcompare/intersetionGene_anno
ls 6.DEGcompare/*_intersetion.xls|while read id;
do
    prefix=`basename $id`
    cat $id|fgrep ∩|cut -f5|tr " " "\n"|sed 's/,//g' > tmp.gene
    cat $anno|fgrep -w -f tmp.gene > 6.DEGcompare/intersetionGene_anno/${prefix%_*}_intersetionGene_anno.xls && rm tmp.gene
done

echo `date` ": All Done!" > 6.DEGcompare/flag.log


# /home/chengyu/soft/miniconda3/bin/ParaFly -c 6.DEGcompare/jobs -CPU 4
