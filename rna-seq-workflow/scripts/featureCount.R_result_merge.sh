#!/bin/bash
#########################################################################
# File Name: /home/chengyu/myscripts/featureCount.R_result_merge.sh
# Author: ChengYu
# Description: 
# Created Time: Fri Jul  7 17:06:08 2023
#########################################################################

DIR=${1}

cd $DIR
files=`ls *.count`
ls ${files}|cut -d"." -f1|sed 's/-/_/g'|tr '\n'  '\t' > samples 
sed -i '1 s/^/id\t/' samples
sed -i '1 s/\t$/\n/' samples 
/home/chengyu/myscripts/njoin.sh `ls ${files}` > all.feacount
/home/chengyu/myscripts/njoin.sh `ls ${files}` > all.log

num=`cat all.feacount|awk '{print NF}'|head -1`
cat all.feacount|cut -f1,`seq -s "," 5 4 ${num}`|sed '1d' > tmp.tpm 
cat all.feacount|cut -f1,`seq -s "," 4 4 ${num}`|sed '1d' > tmp.fpkm
cat all.feacount|cut -f1,`seq -s "," 3 4 ${num}`|sed '1d' > tmp.rawcounts
cat samples all.log > tmp.log
/home/chengyu/myscripts/transposition.sh tmp.log > GeneCount_Assigned_logs.xls
cat samples tmp.tpm > GeneExpression_TPM.xls && rm tmp.tpm 
cat samples tmp.fpkm > GeneExpression_FPKM.xls && rm tmp.fpkm
cat samples tmp.rawcounts > count.matrix.tsv && rm tmp.rawcounts

