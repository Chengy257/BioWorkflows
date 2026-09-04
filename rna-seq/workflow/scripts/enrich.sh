#!/bin/bash
#########################################################################
# DEG 功能富集：拆分 Up/Down/All DEGs → GO/KEGG 富集 + GSEA → DEGs 注释表
# Usage: enrich.sh <DEG_dir> <species> <annotation_tsv> <orgdb_tarball> <kegg_organism>
#   DEG_dir        runDESeq2 输出目录（含 Diff_Expr_Analysis_Results/*_DESeq2.output.tsv）
#   species        osa | hsa
#   annotation_tsv 基因 id 功能注释表
#   orgdb_tarball  osa 本地 OrgDb tarball（species=hsa 时忽略）
#   kegg_organism  KEGG 物种代码（osa→dosa，hsa→hsa）
#########################################################################

set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DIR=$1
spe=$2
annotation=$3
orgdb_tar=$4
kegg_org=$5

function enrich() {
	mkdir -p "$DIR/GO_KEGG_enrich" "$DIR/DEGs"
	ls "$DIR"/Diff_Expr_Analysis_Results/*_DESeq2.output.tsv | while read id;
	do
	## 拆分 DEGs
		b_name=$(basename "$id")
		prefix=${b_name%_*}
		awk '$NF=="Up"{print $1}' "$id" > "$DIR/DEGs/${prefix}_UP.DEGs.txt"
		awk '$NF=="Down"{print $1}' "$id" > "$DIR/DEGs/${prefix}_DOWN.DEGs.txt"
		cat "$DIR/DEGs/${prefix}_UP.DEGs.txt" "$DIR/DEGs/${prefix}_DOWN.DEGs.txt" > "$DIR/DEGs/${prefix}_ALL.DEGs.txt"

	## 提取 FoldChange（供 GSEA）
		cut -f1,3 "$id" > "$DIR/DEGs/${prefix}_FoldChange.xls"

	## GO & KEGG 富集
		for i in UP DOWN ALL;
		do
			mkdir -p "$DIR/GO_KEGG_enrich/${prefix}/${i}"
			Rscript "$SCRIPT_DIR/run_enrichment.R" \
				"$DIR/DEGs/${prefix}_${i}.DEGs.txt" "$DIR/GO_KEGG_enrich/$prefix/${i}" "$spe" "$orgdb_tar" "$kegg_org"
		done
	done

	## DEGs 功能注释
	(
	cd "$DIR/DEGs"
	ls *.DEGs.txt | while read f;
	do
		cat "$annotation" | fgrep -w -f "$f" > "${f}.annotation.tsv"
	done
	)
	ls "$DIR"/DEGs/*_FoldChange.xls > "$DIR/DEGs/all_fc_filespath"

	## GSEA（GO）
	Rscript "$SCRIPT_DIR/run_gsea.R" "$DIR/DEGs/all_fc_filespath" "$DIR/GSEA_enrich_GO" "$spe" "$orgdb_tar"
}

## main
enrich
echo "$(date) : All Done!" > "$DIR/GO_KEGG_enrich/flag.log"
