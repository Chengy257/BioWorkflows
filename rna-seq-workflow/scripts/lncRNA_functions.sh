#!/bin/bash
#########################################################################
# File Name: scripts/lncRNA_functions.sh
# Author: ChengYu
# Description: 
# Created Time: Tue 12 Sep 2023 04:17:50 PM CST
#########################################################################
## Description:
## Define the functions using in lncRNA denovo identification pipeline

CPC2=/home/chengyu/soft/CPC2_standalone-1.0.1/bin/CPC2.py
CNCI_dir=/home/chengyu/soft/github_source/CNCI/
PFAM_DB=/home/chengyu/data/database/pfam/

########################
# runStringtie() {
#     # cat ${ID}|while read id;    
#     # do  
#     # stringtie -p ${threads} --rf -o ${DIR}/4.Assembly/${id}.gtf -G ${GTF_REF} ${DIR}/3.align/${id}_Aligned.sortedByCoord.out.bam 
#     # echo "["`date +"%Y-%m-%d %H:%M.%S"`" Finished stringtie assembling of " ${id} " !]"
#     # done
#     stringtie --merge -p ${threads} -c 0 -F 0 -T 0 -o ${DIR}/4.Assembly/${NAME}_merged.gtf -G ${GTF_REF} `cat ${ID}|xargs -i ls ${DIR}/4.Assembly/{}.gtf|tr "\n" " "` 
#     echo "["`date +"%Y-%m-%d %H:%M.%S"`" Assembling: Stringite program finished!]"
# }
########################
getFasta() {
    ##  Usage: getFasta [input.gtf] [genome.fa]
    gffread -w ${1%.*}.fa -g ${2} $1
}
########################
runGFFcompare() {
# Usage: runGFFcompare [REF gtf file] [gtf file] [genome fasta file]
    gffcompare -T -r ${1} -o ${2}.compare ${2}
    cat ${2}.compare.tracking |awk '$4=="u"{split($5,a,"|");print a[2]}'| \
        sort -u > ${2}.compare.classcode_u.ID
    cat ${2}.compare.annotated.gtf |fgrep -w -f ${2}.compare.classcode_u.ID > ${2}.compare.classcode_u.gtf
    getFasta ${2}.compare.classcode_u.gtf ${3}
}
########################
runCPC2() { 
    ## 
    $CPC2 -i ${1} -o ${1}.CPC2.out
    cat ${1}.CPC2.out.txt|awk '$8=="noncoding"{print $1}'|sort -u > ${1}.CPC2.noncodingID
    }
########################
runCNCI() {
    ## Usage: runCNCI [input fasta,relative path] [threads]
    DIR_OLD=`pwd`
    cd ${CNCI_dir} 
    ./CNCI.py -f ${DIR_OLD}/${1} -o ${DIR_OLD}/${1}.CNCI.tmp -m pl -p ${2}
    # mv ${NAME}.CNCI.tmp/CNCI.index ${DIR_OLD}/${NAME}.CNCI.out.txt && rm -rf ${NAME}.CNCI.tmp 
    cat ${DIR_OLD}/${1}.CNCI.tmp/CNCI.index|awk '$2=="noncoding"{print $1}'|sort -u >${DIR_OLD}/${1}.CNCI.noncodingID 
    # rm ./*log 
    cd ${DIR_OLD}
    }
########################
runPfam() {
    pfam_scan.pl -translate -fasta $1 -dir ${PFAM_DB} -outfile ${NAME}.pfam_scan.out -as -cpu `echo 2*$threads|bc`
    cat ${NAME}.pfam_scan.out|grep -v "^#"|grep -v '^\s*$'|awk '($13<1e-5){print $1}'|sed -e 's/\.[1-9]*$//g'|sort -u > tmp.pfam.hitID
    }
########################

