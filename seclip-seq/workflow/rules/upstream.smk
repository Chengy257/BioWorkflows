# Upstream steps: UMI extraction -> two-pass adapter trimming -> read sort
# by id -> FastQC. Depends on workflow/rules/common.smk for helpers.

rule umi_extract:
    input:
        fq=lambda wc: raw_fastq(wc.sample),
    output:
        fq=R("2.cleandata/{sample}_umi.fq.gz"),
        metrics=R("2.cleandata/logs/{sample}_umi_extract.metrics"),
    params:
        pattern=config["umi"]["pattern"],
    log:
        R("logs/umi_extract/{sample}.log"),
    threads: rthreads("umi_extract")
    resources:
        mem_mb=rmem("umi_extract"),
        runtime_min=rruntime("umi_extract"),
        runtime_sec=rruntime_sec("umi_extract"),
    shell:
        """
        umi_tools extract --random-seed 1 --bc-pattern {params.pattern} \
            --stdin {input.fq} --stdout {output.fq} --log {output.metrics} > {log} 2>&1
        """


rule cutadapt_trim:
    # Two-pass 3' trimming (faithful to the legacy pipeline): a relaxed -O 1
    # sweep against the 20 shifted adapter variants, then a strict -O 5 pass.
    input:
        R("2.cleandata/{sample}_umi.fq.gz"),
    output:
        fq=R("2.cleandata/{sample}_clean.fqTrTr.fq.gz"),
        metrics=R("2.cleandata/logs/{sample}_cutadapt.metrics"),
    params:
        min_len=config["cutadapt"]["min_len"],
        quality=config["cutadapt"]["quality_cutoff"],
        error_rate=config["cutadapt"]["error_rate"],
        adapters=config["cutadapt"]["adapters"],
        tmp_fq=lambda wc: f"{RD}2.cleandata/{wc.sample}_clean.fqTr.fq.gz",
    log:
        R("logs/cutadapt/{sample}.log"),
    threads: rthreads("cutadapt_trim")
    resources:
        mem_mb=rmem("cutadapt_trim"),
        runtime_min=rruntime("cutadapt_trim"),
        runtime_sec=rruntime_sec("cutadapt_trim"),
    shell:
        """
        adapter_args=$(printf -- '-a %s ' {params.adapters})
        cutadapt -j {threads} -O 1 --match-read-wildcards --times 1 \
            -e {params.error_rate} -q {params.quality} -m {params.min_len} \
            $adapter_args -o {params.tmp_fq} {input} > {log} 2>&1
        cutadapt -j {threads} -O 5 --match-read-wildcards --times 1 \
            -e {params.error_rate} -q {params.quality} -m {params.min_len} \
            $adapter_args -o {output.fq} {params.tmp_fq} > {output.metrics} 2>&1
        rm -f {params.tmp_fq}
        """


rule fastq_sort:
    # seqkit name sort (deterministic read order before STAR). Replaced the
    # legacy ea-utils `fastq-sort --id`: the current bioconda ea-utils build
    # (1.1.2.779) ships no fastq-sort binary.
    input:
        R("2.cleandata/{sample}_clean.fqTrTr.fq.gz"),
    output:
        R("2.cleandata/{sample}_clean.fqTrTr.sorted.fq.gz"),
    log:
        R("logs/fastq_sort/{sample}.log"),
    threads: rthreads("fastq_sort")
    resources:
        mem_mb=rmem("fastq_sort"),
        runtime_min=rruntime("fastq_sort"),
        runtime_sec=rruntime_sec("fastq_sort"),
    shell:
        """
        zcat {input} | seqkit sort -n -j {threads} | bgzip -@ {threads} > {output} 2> {log}
        """


rule fastqc:
    input:
        R("2.cleandata/{sample}_clean.fqTrTr.sorted.fq.gz"),
    output:
        html=R("2.cleandata/fastqc/{sample}_fastqc.html"),
        zip=R("2.cleandata/fastqc/{sample}_fastqc.zip"),
    params:
        outdir=lambda wc, output: os.path.dirname(str(output.html)),
    log:
        R("logs/fastqc/{sample}.log"),
    threads: rthreads("fastqc")
    resources:
        mem_mb=rmem("fastqc"),
        runtime_min=rruntime("fastqc"),
        runtime_sec=rruntime_sec("fastqc"),
    shell:
        """
        ## fastqc derives output names from the input file basename
        ## (sample_clean.fqTrTr.sorted_fastqc.*); rename to the
        ## sample-keyed contract paths. NOTE: snakemake formats the whole
        ## shell string, comments included — never write bare braces here.
        fastqc -f fastq -t {threads} -o {params.outdir} {input} > {log} 2>&1
        mv -f {params.outdir}/{wildcards.sample}_clean.fqTrTr.sorted_fastqc.html {output.html}
        mv -f {params.outdir}/{wildcards.sample}_clean.fqTrTr.sorted_fastqc.zip {output.zip}
        """
