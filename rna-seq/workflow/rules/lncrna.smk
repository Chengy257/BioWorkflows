###############################################
## lncRNA de novo 鉴定规则（仅 pipeline=lncrna）：
## stringtie 组装（严格模式）→ gffcompare(classcode u) → CPC2/CNCI/长度编码势过滤
## → Pfam/NR 过滤 → 最终 lncRNA 及其表达定量
## 外部工具路径见 config.yaml 的 lncrna 段（cpc2_bin / cnci_dir / pfam_DB 等）
###############################################

rule runStringtie:
    input:
        bam=R("3.align/{sample}_Aligned.sortedByCoord.out.bam"),
        ref=lres("gtf"),
        strandedness=R("3.align/{sample}.strandedness"),
    output:
        R("4.LncRNA/4.1.Assembly_stringtie/{sample}.gtf"),
    params:
        extra=" --conservative -j 5 ",
    conda:
        os.path.join(ENVS, "lncrna.yaml")
    log:
        R("logs/runStringtie/{sample}_stringtie.log.txt"),
    shell:
        """
        mkdir -p {RD}4.LncRNA/4.1.Assembly_stringtie
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
        gtfIN=expand(R("4.LncRNA/4.1.Assembly_stringtie/{sample}.gtf"), sample=SAMPLES),
        ref=lres("gtf_PcGs"),
        genome=res("genome"),
    output:
        merged=R("4.LncRNA/4.1.Assembly_stringtie/merged.gtf"),
        classcode_fa=R("4.LncRNA/4.1.Assembly_stringtie/merged.gtf.compare.classcode_u.fa"),
    conda:
        os.path.join(ENVS, "lncrna.yaml")
    log:
        R("logs/gtf_merge_compare.log.txt"),
    shell:
        """
        source {SCRIPTS}/lncRNA_functions.sh
        ## 合并各样本 gtf
        stringtie --merge -G {input.ref} -F 0.1 -T 0.1 -o {output.merged} {input.gtfIN} >> {log} 2>&1
        ## 与参考蛋白编码注释比较，取 classcode "u" 的潜在 lincRNA
        runGFFcompare {input.ref} {output.merged} {input.genome} >> {log} 2>&1
        """


rule Coding_predict:
    input:
        fa=R("4.LncRNA/4.1.Assembly_stringtie/merged.gtf.compare.classcode_u.fa"),
        gtfIN=R("4.LncRNA/4.1.Assembly_stringtie/merged.gtf"),
        genome=res("genome"),
    output:
        ids=R("4.LncRNA/4.2.Coding_predict/noncoding.transciptIDs"),
        gtf_nc=R("4.LncRNA/4.2.Coding_predict/noncoding.transciptIDs.gtf"),
        fa_nc=R("4.LncRNA/4.2.Coding_predict/noncoding.transciptIDs.fa"),
    log:
        R("logs/Coding_predict.log.txt"),
    conda:
        os.path.join(ENVS, "lncrna.yaml")
    shell:
        """
        mkdir -p {RD}4.LncRNA/4.2.Coding_predict
        ## 外部工具路径由 config 注入（见 config.yaml 的 lncrna 段）
        export CPC2={config[lncrna][cpc2_bin]}
        export CNCI_dir={config[lncrna][cnci_dir]}
        export PFAM_DB={config[lncrna][pfam_DB]}
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
        cat {input.fa} | bioawk -c fastx 'length($seq)>=200{{print $name}}' > {RD}4.LncRNA/4.2.Coding_predict/Tr_Len_Filtered.transciptIDs
        ## 合并 CPC2/CNCI/长度三重过滤
        cat {RD}4.LncRNA/4.1.Assembly_stringtie/merged.gtf.compare.classcode_u.fa*noncodingID | fgrep -w -f {RD}4.LncRNA/4.2.Coding_predict/Tr_Len_Filtered.transciptIDs | sort -u > {output.ids}
        ## 回取序列（getFasta 生成 noncoding.transciptIDs.fa）
        cat {input.gtfIN} | fgrep -w -f {output.ids} > {output.gtf_nc}
        getFasta {output.gtf_nc} {input.genome} >> {log} 2>&1
        """


rule Pfam_search:
    input:
        fa=R("4.LncRNA/4.2.Coding_predict/noncoding.transciptIDs.fa"),
    output:
        hitIDs=R("4.LncRNA/4.3.Pfam_search/pfam_scan.hitIDs"),
        raw=R("4.LncRNA/4.3.Pfam_search/pfam_scan.output"),
    log:
        R("logs/Pfam_search.log.txt"),
    conda:
        os.path.join(ENVS, "lncrna.yaml")
    threads:
        int((config.get("lncrna") or {}).get("threads", config.get("threads", 30)))
    shell:
        """
        ## pfam_scan.pl 需在 PATH 中（见 envs/lncrna.yaml 说明）
        mkdir -p {RD}4.LncRNA/4.3.Pfam_search
        pfam_scan.pl -translate -fasta {input.fa} -dir {config[lncrna][pfam_DB]} \\
            -outfile {output.raw} -as -cpu {threads} >> {log} 2>&1
        ## 解析显著结构域命中（E-value < 1e-5）
        cat {output.raw} | grep -v "^#" | grep -v '^[[:space:]]*$' | \\
            awk '($13 < 1e-5){{print $1}}' | awk -F"." -v OFS="." '{{$4=null;print}}' | sed 's/\\.$//g' | sort -u > {output.hitIDs}
        """


rule Nr_search:
    input:
        fa=R("4.LncRNA/4.2.Coding_predict/noncoding.transciptIDs.fa"),
    output:
        hitIDs=R("4.LncRNA/4.4.Nr_search/nr.hitIDs"),
        raw=R("4.LncRNA/4.4.Nr_search/diamond.output"),
    log:
        R("logs/Nr_search/log.txt"),
    conda:
        os.path.join(ENVS, "lncrna.yaml")
    threads:
        int((config.get("lncrna") or {}).get("threads", config.get("threads", 30)))
    shell:
        """
        mkdir -p {RD}4.LncRNA/4.4.Nr_search
        diamond blastx --threads {threads} --db {config[lncrna][nr_diamond_DB]} -q {input.fa} \\
            --out {output.raw} >> {log} 2>&1
        cat {output.raw} | awk '{{print $1}}' | sort -u > {output.hitIDs}
        """


rule Final_lncRNA:
    input:
        pfam=R("4.LncRNA/4.3.Pfam_search/pfam_scan.hitIDs"),
        nr=R("4.LncRNA/4.4.Nr_search/nr.hitIDs"),
        gtf_noncode=R("4.LncRNA/4.2.Coding_predict/noncoding.transciptIDs.gtf"),
        genome=res("genome"),
    output:
        fa=R("4.LncRNA/4.5.Final_lncRNA/final_lncRNA.fa"),
        gtf=R("4.LncRNA/4.5.Final_lncRNA/final_lncRNA.gtf"),
        ref_with_lnc=R("4.LncRNA/4.5.Final_lncRNA/ref_withLncRNA.gtf"),
        hits=R("4.LncRNA/4.5.Final_lncRNA/all.hitsIDs"),
    log:
        R("logs/Final_lnc/log.txt"),
    params:
        ref_gtf=res("gtf"),
    conda:
        os.path.join(ENVS, "lncrna.yaml")
    shell:
        """
        mkdir -p {RD}4.LncRNA/4.5.Final_lncRNA
        source {SCRIPTS}/lncRNA_functions.sh
        ## 去除有 Pfam/NR 命中的转录本
        cat {input.pfam} {input.nr} | sort -u > {output.hits}
        cat {input.gtf_noncode} | fgrep -w -v -f {output.hits} > {output.gtf}
        getFasta {output.gtf} {input.genome} >> {log} 2>&1
        ## 参考注释 + 新 lncRNA 合并 gtf（用于下游定量）
        cat {params.ref_gtf} {output.gtf} | sortBed -i - > {output.ref_with_lnc}
        """


rule expression:
    input:
        bam=R("3.align/{sample}_Aligned.sortedByCoord.out.bam"),
        gtf=R("4.LncRNA/4.5.Final_lncRNA/ref_withLncRNA.gtf"),
        strand=R("3.align/{sample}.strandedness"),
    output:
        count=R("5.expression/lncRNA/{sample}.count"),
        stat=R("5.expression/lncRNA/{sample}.log"),
    log:
        R("logs/featureCount_R/{sample}.log.txt"),
    params:
        is_pe=is_paired_end,
        prefix=lambda wc: R(f"5.expression/lncRNA/{wc.sample}"),
    conda:
        os.path.join(ENVS, "quant.yaml")
    shell:
        """
        mkdir -p {RD}5.expression/lncRNA
        strandedness=$(head -1 {input.strand} | awk '{{print $1}}')
        if [ "$strandedness" == "firststrand" ]; then
            strand="2"
        elif [ "$strandedness" == "secondstrand" ]; then
            strand="1"
        else
            strand="0"
        fi
        Rscript {SCRIPTS}/run_featurecounts.R -t {config[threads]} -b {input.bam} -g {input.gtf} \\
            -s $strand -i {params.is_pe} -o {params.prefix} >> {log} 2>&1
        """


rule count_merge2:
    input:
        counts=expand(R("5.expression/lncRNA/{sample}.count"), sample=SAMPLES),
        logs=expand(R("5.expression/lncRNA/{sample}.log"), sample=SAMPLES),
    output:
        R("5.expression/lncRNA/count.matrix.tsv"),
        R("5.expression/lncRNA/GeneExpression_TPM.xls"),
        R("5.expression/lncRNA/GeneExpression_FPKM.xls"),
        R("5.expression/lncRNA/GeneCount_Assigned_logs.xls"),
    log:
        R("logs/count_merge/log.txt"),
    shell:
        "python {SCRIPTS}/merge_featurecounts.py {RD}5.expression/lncRNA >> {log} 2>&1"
