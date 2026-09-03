# picard 去重。是否执行由 workflow.smk 的 sample_bam() 按 assay 决定：
# config["dedup"] 中 chip/atac/faire 默认 true，cuttag 默认 false（保留 PCR 重复）。

rule dedup:
    input:
        "3.align/bowtie2/{sample}_sorted.bam",
    output:
        bam="3.align/bowtie2/{sample}_rmdup.bam",
        bai="3.align/bowtie2/{sample}_rmdup.bam.bai",
        metrics="3.align/bowtie2/{sample}_dup_metrics.txt",
    log:
        "logs/dedup/{sample}.log",
    threads: config["threads"]
    shell:
        """
        mkdir -p logs/dedup
        picard MarkDuplicates \
            --REMOVE_DUPLICATES true \
            --INPUT {input} \
            --OUTPUT {output.bam} \
            --METRICS_FILE {output.metrics} \
            --VALIDATION_STRINGENCY LENIENT \
            --CREATE_INDEX false > {log} 2>&1
        samtools index -@ {threads} {output.bam} {output.bai} >> {log} 2>&1
        """
