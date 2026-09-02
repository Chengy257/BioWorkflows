#!/bin/bash
#########################################################################
## Defaulted Settings
# config=/home/chengyu/project/lab/lncRNAwildRice_synteny_230908/OR/config_user_defined.yaml
config=/home/share/workflows/rna-seq/config_basic_defaulted.yaml
smk=/home/chengyu/workflows/snakemake/rna-seq-workflow/RNA-seq_upstream_only.smk
conda_path=/opt/anaconda3/

#########################################################################
## User Defined Settings
DIR=/home/chengyu/project/lab/231107_XTL/mRNA  			## work dir pathway
config=/home/chengyu/project/lab/231107_XTL/mRNA/config_user_defined.yaml	## user defined config file pathway
job_title=rna-xtl

#########################################################################
## 
cd $DIR
snakemake --cluster "qsub -V -N $job_title -l ncpus=12 -j oe " \
	  --use-conda --conda-base-path $conda_path \
	  --keep-going -s $smk -j 3 --cores 36 --configfile $config \
	  -d $DIR |tee -a >> $DIR/snakemake.logs.txt 2>&1

#########################################################################
## End
