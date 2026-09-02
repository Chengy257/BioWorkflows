#!/bin/bash
#########################################################################
# lncRNA de novo 鉴定辅助函数（由 rules/RNA-seq_lncRNA_DenovoIdenti.smk、
# rules/RNA-seq_upstream_AS.smk source 使用）。
# 外部工具路径由调用方通过环境变量注入（配置来源见 config_lncRNA.yaml）:
#   CPC2      CPC2.py 可执行文件
#   CNCI_dir  CNCI 安装目录（含 CNCI.py，需 Python2 运行环境）
#   PFAM_DB   Pfam-A hmmer 数据库目录
# gffread / gffcompare / pfam_scan.pl 从 PATH 调用（见 envs/*.yaml）。
#########################################################################

die() { echo "[lncRNA_functions] ERROR: $*" >&2; exit 1; }

########################
getFasta() {
    ## Usage: getFasta [input.gtf] [genome.fa]  ->  生成 <input 去扩展名>.fa
    gffread -w "${1%.*}.fa" -g "$2" "$1"
}

########################
runGFFcompare() {
    ## Usage: runGFFcompare [REF gtf] [gtf] [genome fasta]
    ## 取 classcode "u"（参考注释中不存在的转录本）并提取序列
    gffcompare -T -r "$1" -o "$2".compare "$2"
    cat "$2".compare.tracking | awk '$4=="u"{split($5,a,"|");print a[2]}' | sort -u > "$2".compare.classcode_u.ID
    cat "$2".compare.annotated.gtf | fgrep -w -f "$2".compare.classcode_u.ID > "$2".compare.classcode_u.gtf
    getFasta "$2".compare.classcode_u.gtf "$3"
}

########################
runCPC2() {
    [ -n "${CPC2:-}" ] || die "未设置 CPC2 环境变量（来源: config_lncRNA.yaml 的 cpc2_bin）"
    ## Usage: runCPC2 [input fasta]
    "$CPC2" -i "$1" -o "$1".CPC2.out
    cat "$1".CPC2.out.txt | awk '$8=="noncoding"{print $1}' | sort -u > "$1".CPC2.noncodingID
}

########################
runCNCI() {
    [ -n "${CNCI_dir:-}" ] || die "未设置 CNCI_dir 环境变量（来源: config_lncRNA.yaml 的 cnci_dir）"
    ## Usage: runCNCI [input fasta,相对路径] [threads]
    dir_old=$(pwd)
    cd "$CNCI_dir" || die "CNCI 目录不存在: $CNCI_dir"
    ./CNCI.py -f "$dir_old/$1" -o "$dir_old/$1.CNCI.tmp" -m pl -p "$2"
    cat "$dir_old/$1".CNCI.tmp/CNCI.index | awk '$2=="noncoding"{print $1}' | sort -u > "$dir_old/$1".CNCI.noncodingID
    cd "$dir_old"
}

########################
runPfam() {
    [ -n "${PFAM_DB:-}" ] || die "未设置 PFAM_DB 环境变量（来源: config_lncRNA.yaml 的 pfam_DB）"
    ## Usage: runPfam [input fasta] [threads]
    pfam_scan.pl -translate -fasta "$1" -dir "$PFAM_DB" -outfile "$1".pfam_scan.out -as -cpu $((2 * $2))
    cat "$1".pfam_scan.out | grep -v "^#" | grep -v '^[[:space:]]*$' | \
        awk '($13<1e-5){print $1}' | sed -e 's/\.[1-9]*$//g' | sort -u > tmp.pfam.hitID
}
