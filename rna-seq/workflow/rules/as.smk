###############################################
## 转录本组装规则（仅 pipeline=as）：StringTie 组装 → 合并 gffcompare → isoform 定量
## 比对阶段由 align.smk 完成（pipeline=as 时 STAR 使用组装优化参数，见 common.smk）
###############################################

rule Assemble:
    input:
        bam=R("3.align/{sample}_Aligned.sortedByCoord.out.bam"),
        ref=res("gtf"),
        strandedness=R("3.align/{sample}.strandedness"),
    output:
        R("4.assembly/stringtie/{sample}.gtf"),
    log:
        R("logs/Assemble/{sample}_stringtie.log.txt"),
    params:
        outdir=lambda wc, output: os.path.dirname(output[0]),
        stringtie=tool("stringtie", "stringtie"),
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
        ## 链型选项：firststrand → --rf，secondstrand → --fr
        if [ "$strandedness" == "firststrand" ]; then
            {params.stringtie} -p {threads} --rf -o {output} -G {input.ref} {input.bam} >> {log} 2>&1
        elif [ "$strandedness" == "secondstrand" ]; then
            {params.stringtie} -p {threads} --fr -o {output} -G {input.ref} {input.bam} >> {log} 2>&1
        else
            {params.stringtie} -p {threads} -o {output} -G {input.ref} {input.bam} >> {log} 2>&1
        fi
        """


rule gtf_merge:
    input:
        gtf=expand(R("4.assembly/stringtie/{sample}.gtf"), sample=SAMPLES),
        ref=res("gtf"),
        genome=res("genome"),
    output:
        R("4.assembly/stringtie/merged.gtf"),
    log:
        R("logs/gtf_merge.log.txt"),
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
        ## 载入辅助函数（runGFFcompare / getFasta）
        source {params.funcs}
        ## 合并各样本 gtf 并与参考注释比较
        {params.stringtie} --merge -G {input.ref} -i -o {output} {input.gtf} >> {log} 2>&1
        runGFFcompare {input.ref} {output} {input.genome} >> {log} 2>&1
        """


rule isoform_expr:
    input:
        bam=R("3.align/{sample}_Aligned.sortedByCoord.out.bam"),
        ref=R("4.assembly/stringtie/merged.gtf"),
        strandedness=R("3.align/{sample}.strandedness"),
    output:
        gtf=R("4.assembly/isoform/{sample}.gtf"),
        tab=R("4.assembly/isoform/{sample}.tab"),
    log:
        R("logs/isoform_expr/{sample}_stringtie.log.txt"),
    params:
        outdir=lambda wc, output: os.path.dirname(output.gtf),
        stringtie=tool("stringtie", "stringtie"),
    threads:
        rthreads("isoform_expr")
    resources:
        mem_mb=rmem("isoform_expr"),
        runtime_min=rruntime("isoform_expr"),
        runtime_sec=rruntime_sec("isoform_expr"),
    shell:
        """
        mkdir -p {params.outdir}
        strandedness=$(head -1 {input.strandedness} | awk '{{print $1}}')
        if [ "$strandedness" == "firststrand" ]; then
            {params.stringtie} -p {threads} --rf -o {output.gtf} -e -G {input.ref} {input.bam} >> {log} 2>&1
        elif [ "$strandedness" == "secondstrand" ]; then
            {params.stringtie} -p {threads} --fr -o {output.gtf} -e -G {input.ref} {input.bam} >> {log} 2>&1
        else
            {params.stringtie} -p {threads} -o {output.gtf} -e -G {input.ref} {input.bam} >> {log} 2>&1
        fi
        ## 提取转录本表达
        cat {output.gtf} | grep -v "^#" | awk -v OFS="\\t" 'BEGIN{{print "transcript","FPKM","TPM"}} {{if($3=="transcript"){{print $12,$(NF-2),$NF}}}}' | sed 's/[;"]//g' > {output.tab}
        """
