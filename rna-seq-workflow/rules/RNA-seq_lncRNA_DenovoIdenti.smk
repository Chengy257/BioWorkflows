###############################################
## lncRNA de novo 鉴定规则：stringtie 组装 → gffcompare(classcode u) →
## CPC2/CNCI 编码潜能 → 长度过滤 → Pfam/NR 过滤 → 最终 lncRNA 及其表达定量
## 依赖 RNA-seq_upstream.smk 先行 include（使用其中 is_paired_end / SCRIPTS / ENVS）
###############################################

rule runStringtie:
    input:
        bam="3.align/{sample}_Aligned.sortedByCoord.out.bam",
        ref=config["gtf"],
        strandedness="3.align/{sample}.strandedness",
    output:
        "4.LncRNA/4.1.Assembly_stringtie/{sample}.gtf",
    params:
        extra=" --conservative -j 5 ",
    conda:
        os.path.join(ENVS, "lncrna.yaml")
    log:
        "logs/runStringtie/{sample}_stringtie.log.txt"
    shell:
        """
        strandedness=$(head -1 {input.strandedness} | awk '{{print $1}}')
        if [ "$strandedness" == "firststrand" ]; then
            stringtie -p {config[threads]} {params.extra} --rf -o {output} -G {input.ref} {input.bam} >> {log} 2>&1
        elif [ "$strandedness" == "secondstrand" ]; then
            stringtie -p {config[threads]} {params.extra} --fr -o {output} -G {input.ref} {input.bam} >> {log} 2>&1
        else
            stringtie -p {config[threads]} {params.extra} -o {output} -G {input.ref} {input.bam} >> {log} 2>&1
        fi
        """


rule gtf_merge_compare:
    input:
        gtfIN=expand("4.LncRNA/4.1.Assembly_stringtie/{sample}.gtf", sample=SAMPLES),
        ref=config["gtf_PcGs"],
        genome=config["genome"],
    output:
        "4.LncRNA/4.1.Assembly_stringtie/merged.gtf",
        "4.LncRNA/4.1.Assembly_stringtie/merged.gtf.compare.classcode_u.fa",
    conda:
        os.path.join(ENVS, "lncrna.yaml")
    log:
        "logs/gtf_merge_compare.log.txt"
    shell:
        """
        source {SCRIPTS}/lncRNA_functions.sh
        ## 合并各样本 gtf
        stringtie --merge -G {input.ref} -F 0.1 -T 0.1 -o 4.LncRNA/4.1.Assembly_stringtie/merged.gtf {input.gtfIN} >> {log} 2>&1
        ## 与参考蛋白编码注释比较，取 classcode "u" 的潜在 lincRNA
        runGFFcompare {input.ref} 4.LncRNA/4.1.Assembly_stringtie/merged.gtf {input.genome} >> {log} 2>&1
        """


rule Coding_predict:
    input:
        fa="4.LncRNA/4.1.Assembly_stringtie/merged.gtf.compare.classcode_u.fa",
        gtfIN="4.LncRNA/4.1.Assembly_stringtie/merged.gtf",
        genome=config["genome"],
    output:
        fa_nc="4.LncRNA/4.2.Coding_predict/noncoding.transciptIDs.fa",
        gtf_nc="4.LncRNA/4.2.Coding_predict/noncoding.transciptIDs.gtf",
    log:
        "logs/Coding_predict.log.txt"
    conda:
        os.path.join(ENVS, "lncrna.yaml")
    shell:
        """
        ## 外部工具路径由 config 注入（见 config_lncRNA.yaml）
        export CPC2={config[cpc2_bin]}
        export CNCI_dir={config[cnci_dir]}
        export PFAM_DB={config[pfam_DB]}
        source {SCRIPTS}/lncRNA_functions.sh
        ## CPC2 编码潜能预测
        if [ ! -f "merged.gtf.compare.classcode_u.fa.CPC2.noncodingID" ]; then
            runCPC2 {input.fa} >> {log} 2>&1
        fi
        ## CNCI 编码潜能预测
        if [ ! -f "merged.gtf.compare.classcode_u.fa.CNCI.noncodingID" ]; then
            runCNCI {input.fa} {config[threads]} >> {log} 2>&1
        fi
        ## 转录本长度过滤（>= 200 nt）
        cat {input.fa} | bioawk -c fastx 'length($seq)>=200{{print $name}}' > 4.LncRNA/4.2.Coding_predict/Tr_Len_Filtered.transciptIDs
        ## 合并 CPC2/CNCI/长度三重过滤
        cat 4.LncRNA/4.1.Assembly_stringtie/merged.gtf.compare.classcode_u.fa*noncodingID | fgrep -w -f 4.LncRNA/4.2.Coding_predict/Tr_Len_Filtered.transciptIDs | sort -u > 4.LncRNA/4.2.Coding_predict/noncoding.transciptIDs
        ## 回取序列（getFasta 生成 noncoding.transciptIDs.fa）
        cat {input.gtfIN} | fgrep -w -f 4.LncRNA/4.2.Coding_predict/noncoding.transciptIDs > {output.gtf_nc}
        getFasta {output.gtf_nc} {input.genome} >> {log} 2>&1
        """


rule Pfam_search:
    input:
        fa="4.LncRNA/4.2.Coding_predict/noncoding.transciptIDs.fa",
    output:
        "4.LncRNA/4.3.Pfam_search/pfam_scan.hitIDs",
    log:
        "logs/Pfam_search.log.txt"
    conda:
        os.path.join(ENVS, "lncrna.yaml")
    threads:
        int(config.get("lncrna_threads", 30))
    shell:
        """
        ## pfam_scan.pl 需在 PATH 中（见 envs/lncrna.yaml 说明）
        mkdir -p 4.LncRNA/4.3.Pfam_search
        pfam_scan.pl -translate -fasta {input.fa} -dir {config[pfam_DB]} \\
            -outfile 4.LncRNA/4.3.Pfam_search/pfam_scan.output -as -cpu {threads} >> {log} 2>&1
        ## 解析显著结构域命中（E-value < 1e-5）
        cat 4.LncRNA/4.3.Pfam_search/pfam_scan.output | grep -v "^#" | grep -v '^[[:space:]]*$' | \\
            awk '($13 < 1e-5){{print $1}}' | awk -F"." -v OFS="." '{{$4=null;print}}' | sed 's/\\.$//g' | sort -u > {output}
        """


rule Nr_search:
    input:
        fa="4.LncRNA/4.2.Coding_predict/noncoding.transciptIDs.fa",
    output:
        "4.LncRNA/4.4.Nr_search/nr.hitIDs",
    log:
        "logs/Nr_search/log.txt"
    conda:
        os.path.join(ENVS, "lncrna.yaml")
    threads:
        int(config.get("lncrna_threads", 30))
    shell:
        """
        mkdir -p 4.LncRNA/4.4.Nr_search
        diamond blastx --threads {threads} --db {config[nr_diamond_DB]} -q {input.fa} \\
            --out 4.LncRNA/4.4.Nr_search/diamond.output >> {log} 2>&1
        cat 4.LncRNA/4.4.Nr_search/diamond.output | awk '{{print $1}}' | sort -u > {output}
        """


rule Final_lncRNA:
    input:
        "4.LncRNA/4.3.Pfam_search/pfam_scan.hitIDs",
        "4.LncRNA/4.4.Nr_search/nr.hitIDs",
        gtf_noncode="4.LncRNA/4.2.Coding_predict/noncoding.transciptIDs.gtf",
        genome=config["genome"],
    output:
        fa="4.LncRNA/4.5.Final_lncRNA/final_lncRNA.fa",
        gtf="4.LncRNA/4.5.Final_lncRNA/final_lncRNA.gtf",
    log:
        "logs/Final_lnc/log.txt"
    conda:
        os.path.join(ENVS, "lncrna.yaml")
    shell:
        """
        mkdir -p 4.LncRNA/4.5.Final_lncRNA
        source {SCRIPTS}/lncRNA_functions.sh
        ## 去除有 Pfam/NR 命中的转录本
        cat 4.LncRNA/4.3.Pfam_search/pfam_scan.hitIDs 4.LncRNA/4.4.Nr_search/nr.hitIDs | sort -u > 4.LncRNA/4.5.Final_lncRNA/all.hitsIDs
        cat {input.gtf_noncode} | fgrep -w -v -f 4.LncRNA/4.5.Final_lncRNA/all.hitsIDs > {output.gtf}
        getFasta {output.gtf} {input.genome} >> {log} 2>&1
        ## 参考注释 + 新 lncRNA 合并 gtf（用于下游定量）
        cat {config[gtf]} {output.gtf} | sortBed -i - > 4.LncRNA/4.5.Final_lncRNA/ref_withLncRNA.gtf
        """


rule expression:
    input:
        bam="3.align/{sample}_Aligned.sortedByCoord.out.bam",
        gtf="4.LncRNA/4.5.Final_lncRNA/ref_withLncRNA.gtf",
        strand="3.align/{sample}.strandedness",
    output:
        count="5.expression/lncRNA/{sample}.count",
        stat="5.expression/lncRNA/{sample}.log",
    log:
        "logs/featureCount_R/{sample}.log.txt"
    params:
        is_pe=is_paired_end,
    conda:
        os.path.join(ENVS, "quant.yaml")
    shell:
        """
        strandedness=$(head -1 {input.strand} | awk '{{print $1}}')
        if [ "$strandedness" == "firststrand" ]; then
            strand="2"
        elif [ "$strandedness" == "secondstrand" ]; then
            strand="1"
        else
            strand="0"
        fi
        Rscript {SCRIPTS}/run-featurecounts.R -t {config[threads]} -b {input.bam} -g {input.gtf} \\
            -s $strand -i {params.is_pe} -o 5.expression/lncRNA/{wildcards.sample} >> {log} 2>&1
        """


rule count_merge2:
    input:
        counts=expand("5.expression/lncRNA/{sample}.count", sample=SAMPLES),
        logs=expand("5.expression/lncRNA/{sample}.log", sample=SAMPLES),
    output:
        "5.expression/lncRNA/count.matrix.tsv",
        "5.expression/lncRNA/GeneExpression_TPM.xls",
        "5.expression/lncRNA/GeneExpression_FPKM.xls",
        "5.expression/lncRNA/GeneCount_Assigned_logs.xls",
    log:
        "logs/count_merge/log.txt"
    shell:
        "python {SCRIPTS}/merge_featurecounts.py 5.expression/lncRNA >> {log} 2>&1"
