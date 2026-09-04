# picard deduplication. Whether it runs is decided per assay by sample_bam()
# in workflow/rules/common.smk: config["dedup"] defaults chip/atac/faire to
# true and cuttag to false (PCR duplicates are kept).

rule dedup:
    input:
        R("3.align/bowtie2/{sample}_sorted.bam"),
    output:
        bam=R("3.align/bowtie2/{sample}_rmdup.bam"),
        bai=R("3.align/bowtie2/{sample}_rmdup.bam.bai"),
        metrics=R("3.align/bowtie2/{sample}_dup_metrics.txt"),
    log:
        R("logs/dedup/{sample}.log"),
    threads: rthreads("dedup")
    resources:
        mem_mb=rmem("dedup"),
        runtime_min=rruntime("dedup"),
        runtime_sec=rruntime_sec("dedup"),
    shell:
        """
        picard MarkDuplicates \
            --REMOVE_DUPLICATES true \
            --INPUT {input} \
            --OUTPUT {output.bam} \
            --METRICS_FILE {output.metrics} \
            --VALIDATION_STRINGENCY LENIENT \
            --CREATE_INDEX false > {log} 2>&1
        samtools index -@ {threads} {output.bam} {output.bai} >> {log} 2>&1
        """
