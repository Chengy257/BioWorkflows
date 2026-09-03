# 上游公共步骤：trim_galore 修剪 → fastqc → multiqc → bowtie2 建索引 → 比对
# 依赖 workflow.smk 提供的 SAMPLES / config

rule trim_adapter:
    input:
        fq1="1.rawdata/{sample}_1.fq.gz",
        fq2="1.rawdata/{sample}_2.fq.gz",
    output:
        fq1="2.cleandata/{sample}_1_val_1.fq.gz",
        fq2="2.cleandata/{sample}_2_val_2.fq.gz",
    params:
        quality=config["trim"]["quality"],
        stringency=config["trim"]["stringency"],
        error=config["trim"]["error_rate"],
        extra=config["trim"]["extra"],
    log:
        "logs/trim_galore/{sample}.log",
    threads: config["threads"]
    shell:
        """
        mkdir -p 2.cleandata logs/trim_galore
        trim_galore -q {params.quality} --stringency {params.stringency} \
            -e {params.error} {params.extra} --gzip -j {threads} \
            -o 2.cleandata/ --paired {input.fq1} {input.fq2} > {log} 2>&1
        """


rule fastqc:
    input:
        fq1="2.cleandata/{sample}_1_val_1.fq.gz",
        fq2="2.cleandata/{sample}_2_val_2.fq.gz",
    output:
        html1="2.cleandata/fastqc/{sample}_1_val_1_fastqc.html",
        html2="2.cleandata/fastqc/{sample}_2_val_2_fastqc.html",
        zip1="2.cleandata/fastqc/{sample}_1_val_1_fastqc.zip",
        zip2="2.cleandata/fastqc/{sample}_2_val_2_fastqc.zip",
    log:
        "logs/fastqc/{sample}.log",
    threads: config["threads"]
    shell:
        """
        mkdir -p 2.cleandata/fastqc logs/fastqc
        fastqc -f fastq -t {threads} -o 2.cleandata/fastqc/ {input} > {log} 2>&1
        """


rule multiqc:
    input:
        fastqc=expand("2.cleandata/fastqc/{sample}_{r}_val_{r}_fastqc.zip",
                      sample=SAMPLES, r=["1", "2"]),
        bowtie2_logs=expand("logs/bowtie2_mapping/{sample}.log", sample=SAMPLES),
        dup_metrics=expand("3.align/bowtie2/{sample}_dup_metrics.txt",
                           sample=[s for s in SAMPLES if assay_needs_dedup(SEQTYPE_OF[s])]),
        frip_mqc=lambda wc: (["5.QC/frip/FRiP_mqc.tsv"]
                              if config["qc"]["frip"] else []),
        spp_mqc=lambda wc: (["5.QC/spp/NSC_RSC_mqc.tsv"]
                             if config["qc"]["nsc_rsc"] else []),
    output:
        "2.cleandata/fastqc/multiqc/multiqc_report.html",
    log:
        "logs/multiqc.log",
    threads: 1
    shell:
        """
        mkdir -p 2.cleandata/fastqc/multiqc
        multiqc --force -o 2.cleandata/fastqc/multiqc {input} > {log} 2>&1
        """


rule bowtie2_index:
    input:
        config["genome_fa"],
    output:
        "0.index/bowtie2.1.bt2",
    log:
        "logs/bowtie2_index.log",
    threads: config["threads"]
    shell:
        """
        mkdir -p 0.index logs
        bowtie2-build --threads {threads} {input} 0.index/bowtie2 > {log} 2>&1
        """


rule bowtie2_mapping:
    input:
        index="0.index/bowtie2.1.bt2",
        fq1="2.cleandata/{sample}_1_val_1.fq.gz",
        fq2="2.cleandata/{sample}_2_val_2.fq.gz",
    output:
        bam="3.align/bowtie2/{sample}_sorted.bam",
        bai="3.align/bowtie2/{sample}_sorted.bam.bai",
    params:
        extra=config["bowtie2_extra"],
        min_mapq=config["min_mapq"],
    log:
        "logs/bowtie2_mapping/{sample}.log",
    threads: config["threads"]
    shell:
        """
        mkdir -p 3.align/bowtie2 logs/bowtie2_mapping
        # 比对，未比对上的双端读对单独存放；比对摘要写入日志
        bowtie2 -x 0.index/bowtie2 -p {threads} {params.extra} \
            -1 {input.fq1} -2 {input.fq2} \
            --un-conc-gz 3.align/bowtie2/{wildcards.sample}_unmapped.fq.gz \
            2> {log} \
        | samtools view -bS -q {params.min_mapq} - \
        | samtools sort -@ {threads} -O BAM \
            -o 3.align/bowtie2/{wildcards.sample}_sorted.bam -
        samtools index -@ {threads} 3.align/bowtie2/{wildcards.sample}_sorted.bam >> {log} 2>&1
        """
