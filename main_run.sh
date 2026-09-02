#!/bin/bash
#########################################################################
## Defaulted Settings
smk=
## user defined config file pathname 
config=/home/chengyu/workflows/snakemake/cuttag_chip-workflow/config.yaml
#########################################################################
## User Defined Settings
## work dir pathname
DIR=${1}  ## 
## 	job title for PBS jobs
job_title="chip-seq"
## conda enviroment pathname 
conda_path="/opt/anaconda3/"
#########################################################################
## rename 
cd $DIR
if [ -d "1.rawdata" ];then
	cd 1.rawdata	
	rename _R1.fastq.gz _1.fq.gz ./*gz && rename _R2.fastq.gz _2.fq.gz ./*gz	
	rename _1.fastq.gz _1.fq.gz ./*gz && rename _2.fastq.gz _2.fq.gz ./*gz
	rename _f1.fq.gz _1.fq.gz ./*gz && rename _r2.fq.gz _2.fq.gz ./*gz
	cd ../
fi
## main function
snakemake --cluster "qsub -V -N $job_title -l ncpus=6 -j oe " \
	  --use-conda --conda-base-path $conda_path \
	  --keep-going -s $smk -j 3 --cores 18 --configfile $config \
	  -d $DIR 
## clean up 
mv ${job_title}.o* ./logs/ 
rm -rf .snakemake/ 
#########################################################################
## The End