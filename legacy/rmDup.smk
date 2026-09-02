rule dedup_picard:
    input:
        bowtie2_mapping.output
    output:
        "3.align/bowtie2/${id}_rmdup.bam"
    threads:
        {config["threads"]}
    log:
        "logs/dedup_picard/{sample}.logs.txt"
    conda:
        """
            picard MarkDuplicates REMOVE_DUPLICATES=true I=3.align/{wildcards.sample}_sorted.bam O=3.align/{wildcards.sample}_rmdup.bam M=3.align/{wildcards.sample}_metrics >> {log} 2>&1
            samtools index -@ {threads} 3.align/bowtie2/{wildcards.sample}_rmdup.bam >> {log} 2>&1
        """