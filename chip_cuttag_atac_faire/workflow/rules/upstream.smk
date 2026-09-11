# Upstream shared steps: trim_galore trimming -> fastqc -> multiqc ->
# bowtie2 index building -> alignment.
# Depends on workflow/rules/common.smk for SAMPLES / config / helpers.

rule trim_adapter:
    input:
        # Raw FASTQ naming variants (_1/_2 or _R1/_R2, .fq.gz or .fastq.gz)
        # are resolved by raw_fastq_pair() in common.smk (first complete pair
        # wins; clear error when nothing matches).
        fq1=lambda wc: raw_fastq_pair(wc.sample)[0],
        fq2=lambda wc: raw_fastq_pair(wc.sample)[1],
    output:
        fq1=R("2.cleandata/{sample}_1_val_1.fq.gz"),
        fq2=R("2.cleandata/{sample}_2_val_2.fq.gz"),
    params:
        quality=config["trim"]["quality"],
        stringency=config["trim"]["stringency"],
        error=config["trim"]["error_rate"],
        extra=config["trim"]["extra"],
        outdir=lambda wc, output: os.path.dirname(output.fq1),
    log:
        R("logs/trim_galore/{sample}.log"),
    threads: rthreads("trim_adapter")
    resources:
        mem_mb=rmem("trim_adapter"),
        runtime_min=rruntime("trim_adapter"),
        runtime_sec=rruntime_sec("trim_adapter"),
    shell:
        """
        trim_galore -q {params.quality} --stringency {params.stringency} \
            -e {params.error} {params.extra} --gzip -j {threads} \
            -o {params.outdir}/ --paired {input.fq1} {input.fq2} > {log} 2>&1
        """


rule fastqc:
    input:
        fq1=R("2.cleandata/{sample}_1_val_1.fq.gz"),
        fq2=R("2.cleandata/{sample}_2_val_2.fq.gz"),
    output:
        html1=R("2.cleandata/fastqc/{sample}_1_val_1_fastqc.html"),
        html2=R("2.cleandata/fastqc/{sample}_2_val_2_fastqc.html"),
        zip1=R("2.cleandata/fastqc/{sample}_1_val_1_fastqc.zip"),
        zip2=R("2.cleandata/fastqc/{sample}_2_val_2_fastqc.zip"),
    params:
        fastqc_dir=lambda wc, output: os.path.dirname(output.html1),
    log:
        R("logs/fastqc/{sample}.log"),
    threads: rthreads("fastqc")
    resources:
        mem_mb=rmem("fastqc"),
        runtime_min=rruntime("fastqc"),
        runtime_sec=rruntime_sec("fastqc"),
    shell:
        """
        fastqc -f fastq -t {threads} -o {params.fastqc_dir} {input} > {log} 2>&1
        """


rule multiqc:
    input:
        fastqc=expand(R("2.cleandata/fastqc/{sample}_{r}_val_{r}_fastqc.zip"),
                      sample=SAMPLES, r=["1", "2"]),
        bowtie2_logs=expand(R("logs/bowtie2_mapping/{sample}.log"), sample=SAMPLES),
        dup_metrics=expand(R("3.align/bowtie2/{sample}_dup_metrics.txt"),
                           sample=[s for s in SAMPLES if assay_needs_dedup(SEQTYPE_OF[s])]),
        frip_mqc=lambda wc: ([R("5.QC/frip/FRiP_mqc.tsv")]
                             if config["qc"]["frip"] else []),
        spp_mqc=lambda wc: ([R("5.QC/spp/NSC_RSC_mqc.tsv")]
                            if config["qc"]["nsc_rsc"] else []),
        replicate_mqc=lambda wc: ([R("5.QC/replicate_peaks/Replicate_summary_mqc.tsv")]
                                  if REPLICATE["enabled"] else []),
        tss_mqc=lambda wc: ([R("5.QC/tss/TSSE_summary_mqc.tsv")]
                            if (QC_TSS and TSS_SAMPLES) else []),
        organelle_mqc=lambda wc: ([R("5.QC/organelle/Organelle_summary_mqc.tsv")]
                                  if QC_ORGANELLE else []),
        gates_mqc=lambda wc: ([R("5.QC/gates/gate_summary_mqc.tsv")]
                              if GATES["enabled"] else []),
        spikein_mqc=lambda wc: ([R("5.QC/spike_in/Spikein_summary_mqc.tsv")]
                                if SPIKE_IN["enabled"] else []),
    output:
        R("2.cleandata/fastqc/multiqc/multiqc_report.html"),
    params:
        mqc_config=os.path.join(WORKFLOW_DIR, "multiqc_config.yaml"),
        outdir=lambda wc, output: os.path.dirname(str(output)),
    log:
        R("logs/multiqc.log"),
    threads: rthreads("multiqc")
    resources:
        mem_mb=rmem("multiqc"),
        runtime_min=rruntime("multiqc"),
        runtime_sec=rruntime_sec("multiqc"),
    shell:
        """
        multiqc --force -o {params.outdir} -c {params.mqc_config} {input} > {log} 2>&1
        """


rule bowtie2_index:
    input:
        config["genome_fa"],
    output:
        R("0.index/bowtie2.1.bt2"),
    params:
        # bowtie2-build writes <prefix>.1.bt2 .. <prefix>.rev.2.bt2; only the
        # first file is declared as the DAG sentinel.
        prefix=lambda wc, output: str(output)[:-len(".1.bt2")],
    log:
        R("logs/bowtie2_index.log"),
    threads: rthreads("bowtie2_index")
    resources:
        mem_mb=rmem("bowtie2_index"),
        runtime_min=rruntime("bowtie2_index"),
        runtime_sec=rruntime_sec("bowtie2_index"),
    shell:
        """
        bowtie2-build --threads {threads} {input} {params.prefix} > {log} 2>&1
        """


rule bowtie2_mapping:
    input:
        index=R("0.index/bowtie2.1.bt2"),
        fq1=R("2.cleandata/{sample}_1_val_1.fq.gz"),
        fq2=R("2.cleandata/{sample}_2_val_2.fq.gz"),
    output:
        bam=R("3.align/bowtie2/{sample}_sorted.bam"),
        bai=R("3.align/bowtie2/{sample}_sorted.bam.bai"),
        # The --un-conc-gz pair is declared ALWAYS (not gated by spike_in) so
        # the files carry a DAG edge: the spike_in stage re-aligns them
        # against the spike-in genome. bowtie2 appends .1/.2 to the pattern
        # with the trailing .gz replaced (verified on bowtie2 2.5.1), and it
        # creates both files even when nothing failed to align.
        unmapped1=R("3.align/bowtie2/{sample}_unmapped.fq.1.gz"),
        unmapped2=R("3.align/bowtie2/{sample}_unmapped.fq.2.gz"),
    params:
        extra=config["bowtie2_extra"],
        min_mapq=config["min_mapq"],
        idx_prefix=lambda wc, input: str(input.index)[:-len(".1.bt2")],
        # --un-conc-gz PATTERN (the produced files are the unmapped1/unmapped2
        # outputs declared above)
        unmapped=lambda wc, output: os.path.join(os.path.dirname(str(output.bam)),
                                                 f"{wc.sample}_unmapped.fq.gz"),
    log:
        R("logs/bowtie2_mapping/{sample}.log"),
    threads: rthreads("bowtie2_mapping")
    resources:
        mem_mb=rmem("bowtie2_mapping"),
        runtime_min=rruntime("bowtie2_mapping"),
        runtime_sec=rruntime_sec("bowtie2_mapping"),
    shell:
        """
        # Alignment; unaligned pairs are kept separately and the summary goes to the log.
        bowtie2 -x {params.idx_prefix} -p {threads} {params.extra} \
            -1 {input.fq1} -2 {input.fq2} \
            --un-conc-gz {params.unmapped} \
            2> {log} \
        | samtools view -bS -q {params.min_mapq} - \
        | samtools sort -@ {threads} -O BAM \
            -o {output.bam} -
        samtools index -@ {threads} {output.bam} {output.bai} >> {log} 2>&1
        """
