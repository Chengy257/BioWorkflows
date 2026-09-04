#!/bin/bash
#########################################################################
# lncRNA de novo 鉴定辅助函数（由 rules/lncrna.smk、rules/as.smk source 使用）。
# External tools are resolved from config/software.yaml and exported by run.sh:
#   CPC2 / CNCI_dir / PFAM_DB / PFAM_SCAN are runtime values.
# gffread / gffcompare use RNASEQ_TOOL_GFFREAD / RNASEQ_TOOL_GFFCOMPARE when set.
#########################################################################

die() { echo "[lncRNA_functions] ERROR: $*" >&2; exit 1; }
GFFREAD="${RNASEQ_TOOL_GFFREAD:-gffread}"
GFFCOMPARE="${RNASEQ_TOOL_GFFCOMPARE:-gffcompare}"
CNCI_PYTHON="${RNASEQ_TOOL_CNCI_PYTHON:-python2}"

########################
getFasta() {
    ## Usage: getFasta [input.gtf] [genome.fa]  ->  生成 <input 去扩展名>.fa
    "$GFFREAD" -w "${1%.*}.fa" -g "$2" "$1"
}

########################
runGFFcompare() {
    ## Usage: runGFFcompare [REF gtf] [gtf] [genome fasta]
    ## 取 classcode "u"（参考注释中不存在的转录本）并提取序列
    "$GFFCOMPARE" -T -r "$1" -o "$2".compare "$2"
    cat "$2".compare.tracking | awk '$4=="u"{split($5,a,"|");print a[2]}' | sort -u > "$2".compare.classcode_u.ID
    cat "$2".compare.annotated.gtf | fgrep -w -f "$2".compare.classcode_u.ID > "$2".compare.classcode_u.gtf
    getFasta "$2".compare.classcode_u.gtf "$3"
}

########################
runCPC2() {
    [ -n "${CPC2:-}" ] || die "未设置 CPC2 环境变量（source: software.yaml tools.cpc2）"
    ## Usage: runCPC2 [input fasta]
    "$CPC2" -i "$1" -o "$1".CPC2.out
    cat "$1".CPC2.out.txt | awk '$8=="noncoding"{print $1}' | sort -u > "$1".CPC2.noncodingID
}

########################
runCNCI() {
    [ -n "${CNCI_dir:-}" ] || die "未设置 CNCI_dir 环境变量（source: software.yaml paths.cnci_dir）"
    ## Usage: runCNCI [input fasta,相对路径] [threads]
    dir_old=$(pwd)
    cd "$CNCI_dir" || die "CNCI 目录不存在: $CNCI_dir"
    "$CNCI_PYTHON" ./CNCI.py -f "$dir_old/$1" -o "$dir_old/$1.CNCI.tmp" -m pl -p "$2"
    cat "$dir_old/$1".CNCI.tmp/CNCI.index | awk '$2=="noncoding"{print $1}' | sort -u > "$dir_old/$1".CNCI.noncodingID
    cd "$dir_old" || return 1
}

########################
runPfam() {
    [ -n "${PFAM_DB:-}" ] || die "未设置 PFAM_DB 环境变量（source: software.yaml databases.pfam）"
    ## Usage: runPfam [input fasta] [threads]
    local scan="${PFAM_SCAN:-pfam_scan.pl}"
    "$scan" -translate -fasta "$1" -dir "$PFAM_DB" -outfile "$1".pfam_scan.out -as -cpu $((2 * $2))
    cat "$1".pfam_scan.out | grep -v "^#" | grep -v '^[[:space:]]*$' | \
        awk '($13<1e-5){print $1}' | sed -e 's/\.[1-9]*$//g' | sort -u > tmp.pfam.hitID
}
