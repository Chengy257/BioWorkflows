#!/bin/bash
#########################################################################
# File Name: /home/chengyu/workflows/snakemake/rna-seq-workflow/scripts/check_strandness.sh
# Author: ChengYu
# Description: 
# Created Time: Sun 22 Oct 2023 09:55:16 PM CST
#########################################################################



infer_experiment.py -r {input.bed} -i {input.bam} > 3.align/{wildcards.sample}_infer_experiment.out  
cat 3.align/{wildcards.sample}_infer_experiment.out|tail -2|awk '{{print $NF}}'| \
	awk '{{f1=$0;getline;f2=$0; \
		if(f1-f2 > 0.4) {{print "secondstrand"}}  \
		else if(f2-f1 > 0.4) {{print "firststrand"}} \
		else {{print "unstrand"}} }}' > 3.align/{wildcards.sample}.strandedness 
strandedness=`cat 3.align/{wildcards.sample}.strandedness`


strandedness=`cat {input.strand}|head -1|awk '{{print $1}}'`            
## run featureCount_R depends on  different library strandedness type 
if [ ${{strandedness}} == "firststrand" ];then  ## 
	strand="2"
elif [ ${{strandedness}} == "secondstrand" ];then ## 
	strand="1"
elif [ ${{strandedness}} == "unstrand" ];then
	strand="0"
fi



