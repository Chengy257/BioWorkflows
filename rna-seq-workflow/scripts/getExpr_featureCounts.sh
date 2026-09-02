#!/bin/bash
#########################################################################
# File Name: /home/chengyu/workflows/snakemake/rna-seq-workflow/scripts/getExpr_featureCounts.sh
# Author: ChengYu
# Description: 
# Created Time: Tue Jul 11 17:29:07 2023
#########################################################################

bam=${1}
gtf=${2}
threads=${3}
output=${4}


strandness=0


featureCounts -T $threads -a $gtf -o 

/home/chengyu/soft/miniconda3/envs/lncrna/bin/infer_experiment.py -r /home/chengyu/references/osa/Oryza_sativa.IRGSP-1.0.56.Chr.bed -i ${bam} > temp_strandness
layout=`cat temp_strandness|sed -n '3p'|awk '{print $3}'`
fraction1=`cat temp_strandness|sed -n '5p'|awk '{print $NF}'`
fraction2=`cat temp_strandness|sed -n '6p'|awk '{print $NF}'`
rm -rf temp_strandness

if [ "fraction1" > 0.8 ]
