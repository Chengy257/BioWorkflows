###############################################
## lncRNA de novo 鉴定规则（仅 pipeline=lncrna）：
## stringtie 组装（严格模式）→ gffcompare(classcode u) → CPC2/CNCI/长度编码势过滤
## → Pfam/NR 过滤 → 最终 lncRNA 及其表达定量
## External executables and databases are resolved from config/software.yaml by run.sh.
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
        outdir=lambda wc, output: os.path.dirname(output[0]),
        stringtie=tool("stringtie", "stringtie"),
    log:
        R("logs/runStringtie/{sample}_stringtie.log.txt"),
    threads:
        rthreads("stringtie")
    resources:
        mem_mb=rmem("stringtie"),
        runtime_min=rruntime("stringtie"),
        runtime_sec=rruntime_sec("stringtie"),
    shell:
        """
        mkdir -p {params.outdir}
        strandedness=$(head -1 {input.strandedness} | awk '{{print $1}}')
        if [ "$strandedness" == "firststrand" ]; then
            {params.stringtie} -p {threads} {params.extra} --rf -o {output} -G {input.ref} {input.bam} >> {log} 2>&1
        elif [ "$strandedness" == "secondstrand" ]; then
            {params.stringtie} -p {threads} {params.extra} --fr -o {output} -G {input.ref} {input.bam} >> {log} 2>&1
        else
            {params.stringtie} -p {threads} {params.extra} -o {output} -G {input.ref} {input.bam} >> {log} 2>&1
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
    log:
        R("logs/gtf_merge_compare.log.txt"),
    params:
        funcs=os.path.join(SCRIPTS, "lncRNA_functions.sh"),
        stringtie=tool("stringtie", "stringtie"),
    threads:
        rthreads("gtf_merge")
    resources:
        mem_mb=rmem("gtf_merge"),
        runtime_min=rruntime("gtf_merge"),
        runtime_sec=rruntime_sec("gtf_merge"),
    shell:
        """
        source {params.funcs}
        ## 合并各样本 gtf
        {params.stringtie} --merge -G {input.ref} -F 0.1 -T 0.1 -o {output.merged} {input.gtfIN} >> {log} 2>&1
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
    params:
        outdir=lambda wc, output: os.path.dirname(output.ids),
        assembly_dir=lambda wc, input: os.path.dirname(input.fa),
        cpc2=runtime_env("RNASEQ_TOOL_CPC2", tool("cpc2", "CPC2.py")),
        cnci_dir=runtime_env("RNASEQ_CNCI_DIR", ""),
        pfam_db=runtime_env("RNASEQ_DB_PFAM", ""),
        bioawk=tool("bioawk", "bioawk"),
        funcs=os.path.join(SCRIPTS, "lncRNA_functions.sh"),
    threads:
        rthreads("lnc_coding")
    resources:
        mem_mb=rmem("lnc_coding"),
        runtime_min=rruntime("lnc_coding"),
        runtime_sec=rruntime_sec("lnc_coding"),
    shell:
        """
        mkdir -p {params.outdir}
        ## External tool/database paths are injected by the unified software runtime.
        export CPC2={params.cpc2}
        export CNCI_dir={params.cnci_dir}
        export PFAM_DB={params.pfam_db}
        source {params.funcs}
        ## CPC2 编码潜能预测
        if [ ! -f "merged.gtf.compare.classcode_u.fa.CPC2.noncodingID" ]; then
            runCPC2 {input.fa} >> {log} 2>&1
        fi
        ## CNCI 编码潜能预测
        if [ ! -f "merged.gtf.compare.classcode_u.fa.CNCI.noncodingID" ]; then
            runCNCI {input.fa} {threads} >> {log} 2>&1
        fi
        ## 转录本长度过滤（>= 200 nt）
        cat {input.fa} | {params.bioawk} -c fastx 'length($seq)>=200{{print $name}}' > {params.outdir}/Tr_Len_Filtered.transciptIDs
        ## 合并 CPC2/CNCI/长度三重过滤
        cat {params.assembly_dir}/merged.gtf.compare.classcode_u.fa*noncodingID | fgrep -w -f {params.outdir}/Tr_Len_Filtered.transciptIDs | sort -u > {output.ids}
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
    params:
        pfam_scan=tool("pfam_scan", "pfam_scan.pl"),
        outdir=lambda wc, output: os.path.dirname(output.raw),
        pfam_db=runtime_env("RNASEQ_DB_PFAM", ""),
    threads:
        rthreads("pfam")
    resources:
        mem_mb=rmem("pfam"),
        runtime_min=rruntime("pfam"),
        runtime_sec=rruntime_sec("pfam"),
    shell:
        """
        ## pfam_scan.pl is resolved from the main runtime/PATH or a software.yaml override.
        mkdir -p {params.outdir}
        {params.pfam_scan} -translate -fasta {input.fa} -dir {params.pfam_db} \\
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
    params:
        outdir=lambda wc, output: os.path.dirname(output.raw),
        nr_db=runtime_env("RNASEQ_DB_NR_DIAMOND", ""),
        diamond=tool("diamond", "diamond"),
    threads:
        rthreads("nr")
    resources:
        mem_mb=rmem("nr"),
        runtime_min=rruntime("nr"),
        runtime_sec=rruntime_sec("nr"),
    shell:
        """
        mkdir -p {params.outdir}
        {params.diamond} blastx --threads {threads} --db {params.nr_db} -q {input.fa} \\
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
        outdir=lambda wc, output: os.path.dirname(output.fa),
        funcs=os.path.join(SCRIPTS, "lncRNA_functions.sh"),
        bedtools=tool("bedtools", "bedtools"),
    threads:
        rthreads("lnc_final")
    resources:
        mem_mb=rmem("lnc_final"),
        runtime_min=rruntime("lnc_final"),
        runtime_sec=rruntime_sec("lnc_final"),
    shell:
        """
        mkdir -p {params.outdir}
        source {params.funcs}
        ## 去除有 Pfam/NR 命中的转录本
        cat {input.pfam} {input.nr} | sort -u > {output.hits}
        cat {input.gtf_noncode} | fgrep -w -v -f {output.hits} > {output.gtf}
        getFasta {output.gtf} {input.genome} >> {log} 2>&1
        ## 参考注释 + 新 lncRNA 合并 gtf（用于下游定量）
        cat {params.ref_gtf} {output.gtf} | {params.bedtools} sort -i - > {output.ref_with_lnc}
        """


rule expression:
    input:
        bam=R("3.align/{sample}_Aligned.sortedByCoord.out.bam"),
        gtf=R("4.LncRNA/4.5.Final_lncRNA/ref_withLncRNA.gtf"),
        strand=R("3.align/{sample}.strandedness"),
    output:
        count_file=R("5.expression/lncRNA/{sample}.count"),
        stat=R("5.expression/lncRNA/{sample}.log"),
    log:
        R("logs/featureCount_R/{sample}.log.txt"),
    params:
        is_pe=is_paired_end,
        prefix=lambda wc: R(f"5.expression/lncRNA/{wc.sample}"),
        outdir=lambda wc, output: os.path.dirname(output.count_file),
        script=os.path.join(SCRIPTS, "run_featurecounts.R"),
        rscript=RSCRIPT,
    threads:
        rthreads("lnc_featurecounts")
    resources:
        mem_mb=rmem("lnc_featurecounts"),
        runtime_min=rruntime("lnc_featurecounts"),
        runtime_sec=rruntime_sec("lnc_featurecounts"),
    shell:
        """
        mkdir -p {params.outdir}
        strandedness=$(head -1 {input.strand} | awk '{{print $1}}')
        if [ "$strandedness" == "firststrand" ]; then
            strand="2"
        elif [ "$strandedness" == "secondstrand" ]; then
            strand="1"
        else
            strand="0"
        fi
        {params.rscript} {params.script} -t {threads} -b {input.bam} -g {input.gtf} \\
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
    params:
        script=os.path.join(SCRIPTS, "merge_featurecounts.py"),
        python=PYTHON,
        expression_dir=lambda wc, output: os.path.dirname(output[0]),
    threads:
        rthreads("lnc_count_merge")
    resources:
        mem_mb=rmem("lnc_count_merge"),
        runtime_min=rruntime("lnc_count_merge"),
        runtime_sec=rruntime_sec("lnc_count_merge"),
    shell:
        "{params.python} {params.script} {params.expression_dir} >> {log} 2>&1"
