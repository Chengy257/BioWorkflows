###############################################
## 定量规则（upstream/deg/lncrna 管线使用）：featureCounts → 表达矩阵合并
###############################################

rule featureCount_R:
    input:
        bam=R("3.align/{sample}_Aligned.sortedByCoord.out.bam"),
        gtf=res("gtf"),
        strand=R("3.align/{sample}.strandedness"),
    output:
        count_file=R("4.expression/{sample}.count"),
        stat=R("4.expression/{sample}.log"),
    log:
        R("logs/featureCount_R/{sample}.log.txt"),
    params:
        is_pe=is_paired_end,
        prefix=lambda wc: R(f"4.expression/{wc.sample}"),
        outdir=lambda wc, output: os.path.dirname(output.count_file),
        script=os.path.join(SCRIPTS, "run_featurecounts.R"),
        rscript=RSCRIPT,
    threads:
        rthreads("featurecounts")
    resources:
        mem_mb=rmem("featurecounts"),
        runtime_min=rruntime("featurecounts"),
        runtime_sec=rruntime_sec("featurecounts"),
    shell:
        """
        mkdir -p {params.outdir}
        strandedness=$(head -1 {input.strand} | awk '{{print $1}}')
        ## 链特异性定量：featureCounts strandSpecific（1=整合链，2=反转链）
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


rule count_merge:
    input:
        counts=expand(R("4.expression/{sample}.count"), sample=SAMPLES),
        logs=expand(R("4.expression/{sample}.log"), sample=SAMPLES),
    output:
        R("4.expression/count.matrix.tsv"),
        R("4.expression/GeneExpression_TPM.xls"),
        R("4.expression/GeneExpression_FPKM.xls"),
        R("4.expression/GeneCount_Assigned_logs.xls"),
    log:
        R("logs/count_merge/log.txt"),
    params:
        script=os.path.join(SCRIPTS, "merge_featurecounts.py"),
        python=PYTHON,
        expression_dir=lambda wc, output: os.path.dirname(output[0]),
    threads:
        rthreads("count_merge")
    resources:
        mem_mb=rmem("count_merge"),
        runtime_min=rruntime("count_merge"),
        runtime_sec=rruntime_sec("count_merge"),
    shell:
        "{params.python} {params.script} {params.expression_dir} >> {log} 2>&1"
