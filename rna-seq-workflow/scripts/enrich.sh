#!/bin/bash
#########################################################################
# File Name: /home/chengyu/workflows/snakemake/rna-seq-workflow/scripts/enrich.sh
# Author: ChengYu
# Description: 
# Created Time: Sat Jul  8 17:25:47 2023
#########################################################################


function enrich() {
	mkdir -p $DIR/GO_KEGG_enrich $DIR/DEGs
	ls $DIR/Diff_Expr_Analysis_Reults/*_DESeq2.output.tsv|while read id; 
	do
	## get DEGs
		b_name=`basename $id`
		prefix=${b_name%_*}	
		cat $id |awk '$NF=="Up"{print $1}' > $DIR/DEGs/${prefix}_UP.DEGs.txt
		cat $id |awk '$NF=="Down"{print $1}' > $DIR/DEGs/${prefix}_DOWN.DEGs.txt
		cat $DIR/DEGs/${prefix}_UP.DEGs.txt $DIR/DEGs/${prefix}_DOWN.DEGs.txt > $DIR/DEGs/${prefix}_ALL.DEGs.txt

	## get FoldChange
		cat $id|cut -f1,3 > $DIR/DEGs/${prefix}_FoldChange.xls

	## GO & KEGG enrich	
		for i in UP DOWN ALL ;
		do
			mkdir -p $DIR/GO_KEGG_enrich/${prefix}/${i}
			/usr/bin/Rscript /home/chengyu/workflows/snakemake/rna-seq-workflow/scripts/enrich_GO_KEGG_clusterProfiler_gProfilerGO.R $DIR/DEGs/${prefix}_${i}.DEGs.txt $DIR/GO_KEGG_enrich/$prefix/${i} 
		done	
	## svg 2 pdf
		# ls ./enrich/$prefix/*/*pdf |xargs -i /usr/bin/rsvg-convert -f png -o {}.png {} 
		# ls ./*pdf|xargs -i /usr/bin/rsvg-convert -f png -o {}.png {}
	done
	# DEGs annotation
	cd $DIR/DEGs 
	ls *.DEGs.txt|while read id ; 
	do  
		cat $annotation |fgrep -w -f $id > ${id}.annotation.tsv ;
	done 	
	cd -
	ls $DIR/DEGs/*_FoldChange.xls > $DIR/DEGs/all_fc_filespath
	## GSEA 
	/usr/bin/Rscript /home/chengyu/workflows/snakemake/rna-seq-workflow/scripts/multiGSEA_gProfilerGO_231216.R $DIR/DEGs/all_fc_filespath $DIR/GSEA_enrich_GO osa
}

## main 
DIR=$1
spe=$2
annotation=/home/chengyu/references/osa/rap/20230315/annotation_full.tsv
# cd $DIR
# rm -rf $DIR/enrich $DIR/DEGs && mkdir -p $DIR/enrich $DIR/DEGs

enrich
echo `date` ": All Done!" > $DIR/GO_KEGG_enrich/flag.log

