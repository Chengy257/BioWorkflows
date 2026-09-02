#!/bin/bash
#########################################################################
# File Name: /home/chengyu/workflows/snakemake/chip_cuttag_atac_faire-workflows/scripts/bdgcmp_macs2.sh
# Author: ChengYu
# Description: 
# Created Time: Wed 28 Feb 2024 11:31:43 AM CST
#########################################################################
##
GROUP_name=${1}
DIR=${2}
chromsize="/share/data/reference/osa/chrom.sizes"
## bdgcmp 
macs2 bdgcmp -t ${DIR}/${{GROUP_name}}_treat_pileup.bdg -c ${DIR}/${{GROUP_name}}_control_lambda.bdg -o ${DIR}/${{GROUP_name}}_FE_bdgcmp.bdg -m FE -p 0.00001 
## 
bedtools slop -i ${DIR}/${{GROUP_name}}_FE_bdgcmp.bdg -g ${chromsize} -b 0 |/share/ucsc-tools/bedClip stdin ${chromsize} ${DIR}/${{GROUP_name}}.bdg.clip
## sort 
LC_COLLATE=C sort -k1,1 -k2,2n ${DIR}/${{GROUP_name}}.bdg.clip > ${DIR}/${{GROUP_name}}.bdg.clip.sorted \
&& rm ${DIR}/${{GROUP_name}}.bdg.clip
## bedGraphToBigWig
/share/ucsc-tools/bedGraphToBigWig ${DIR}/${{GROUP_name}}.bdg.clip.sorted ${chromsize} ${DIR}/${{GROUP_name}}_FE_bdgcmp.bw \
&& rm ${DIR}/${{GROUP_name}}.bdg.clip.sorted
echo "["`date`"]: Finished bdgcmp_macs2 programs."
