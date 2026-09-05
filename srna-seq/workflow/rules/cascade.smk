# Sequential bowtie1 filtering cascade + genome alignment.
# Depends on workflow/rules/common.smk for cascade_input and helpers.

rule cascade_stage:
    input:
        fq=cascade_input,
        index=lambda wc: R(f"0.index/{wc.klass}/{wc.klass}.1.ebwt"),
    output:
        sam=R("3.align/filter/{klass}/{sample}.sam"),
        unmapped=R("3.align/filter/{klass}/{sample}_unmapped.fq"),
    params:
        index=lambda wc, input: str(input.index)[:-len(".1.ebwt")],
        extra=config["bowtie"]["extra"],
    log:
        R("logs/cascade/{klass}/{sample}.log"),
    threads: rthreads("cascade_stage")
    resources:
        mem_mb=rmem("cascade_stage"),
        runtime_min=rruntime("cascade_stage"),
        runtime_sec=rruntime_sec("cascade_stage"),
    shell:
        """
        bowtie -p {threads} -q -x {params.index} {input.fq} \
            -S {output.sam} --un {output.unmapped} {params.extra} > {log} 2>&1
        """


rule genome_align:
    # Faithful to the legacy script: the genome is aligned with the TRIMMED
    # reads directly, not with the cascade remainder.
    input:
        fq=R("2.cleandata/{sample}_trimmed.fq.gz"),
        index=R("0.index/genome/genome.1.ebwt"),
    output:
        sam=R("3.align/genome/{sample}.sam"),
        unmapped=R("3.align/genome/{sample}_unmapped.fq"),
    params:
        index=lambda wc, input: str(input.index)[:-len(".1.ebwt")],
        extra=config["bowtie"]["extra"],
    log:
        R("logs/genome_align/{sample}.log"),
    threads: rthreads("genome_align")
    resources:
        mem_mb=rmem("genome_align"),
        runtime_min=rruntime("genome_align"),
        runtime_sec=rruntime_sec("genome_align"),
    shell:
        """
        bowtie -p {threads} -q -x {params.index} {input.fq} \
            -S {output.sam} --un {output.unmapped} {params.extra} > {log} 2>&1
        """
