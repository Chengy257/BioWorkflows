###############################################
## Transcript assembly rules (pipeline=as only): StringTie assembly -> merge + gffcompare -> isoform quantification
## The alignment stage is handled by align.smk (STAR uses assembly-optimized parameters when pipeline=as, see common.smk)
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
        ## Strand options: firststrand -> --rf, secondstrand -> --fr
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
        ## Load helper functions (runGFFcompare / getFasta)
        source {params.funcs}
        ## Merge per-sample GTFs and compare against the reference annotation
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
        ## Extract transcript-level expression
        cat {output.gtf} | grep -v "^#" | awk -v OFS="\\t" 'BEGIN{{print "transcript","FPKM","TPM"}} {{if($3=="transcript"){{print $12,$(NF-2),$NF}}}}' | sed 's/[;"]//g' > {output.tab}
        """
