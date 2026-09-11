# Spike-in normalization (spike_in; included only when the switch is on):
# second-pass alignment of the read pairs that failed the primary alignment
# (the declared unmapped outputs of the main bowtie2 step) against the
# spike-in genome, per-sample idxstats/flagstat QC, and a project-wide
# summary with the scale factors (1e6 / spike-in mapped reads). The
# spike_in.scale_bigwigs wiring lives in callpeak.smk (bigwig rules).

rule spike_bowtie2_index:
    input:
        SPIKE_IN["fasta"],
    output:
        R("0.index/spike_bt2.1.bt2"),
    params:
        # bowtie2-build writes <prefix>.1.bt2 .. <prefix>.rev.2.bt2; only the
        # first file is declared as the DAG sentinel.
        prefix=lambda wc, output: str(output)[:-len(".1.bt2")],
    log:
        R("logs/spike_in/bowtie2_index.log"),
    threads: rthreads("spike_bowtie2_index")
    resources:
        mem_mb=rmem("spike_bowtie2_index"),
        runtime_min=rruntime("spike_bowtie2_index"),
        runtime_sec=rruntime_sec("spike_bowtie2_index"),
    shell:
        """
        bowtie2-build --threads {threads} {input} {params.prefix} > {log} 2>&1
        """


rule spike_align:
    # Re-align one sample's unmapped pair against the spike-in genome. The
    # alignment itself runs with plain bowtie2 defaults: the primary
    # recipe's bowtie2_extra / min_mapq tuning stays untouched (spike-in
    # quantification only needs the mapped count).
    input:
        index=R("0.index/spike_bt2.1.bt2"),
        unmapped1=R("3.align/bowtie2/{sample}_unmapped.fq.1.gz"),
        unmapped2=R("3.align/bowtie2/{sample}_unmapped.fq.2.gz"),
    output:
        bam=R("3.align/spike_in/{sample}_sorted.bam"),
        bai=R("3.align/spike_in/{sample}_sorted.bam.bai"),
        idxstats=R("5.QC/spike_in/{sample}_idxstats.txt"),
        flagstat=R("5.QC/spike_in/{sample}_flagstat.txt"),
    wildcard_constraints:
        sample=_group_regex(SAMPLES),
    params:
        idx_prefix=lambda wc, input: str(input.index)[:-len(".1.bt2")],
    log:
        R("logs/spike_in/{sample}.log"),
    threads: rthreads("spike_align")
    resources:
        mem_mb=rmem("spike_align"),
        runtime_min=rruntime("spike_align"),
        runtime_sec=rruntime_sec("spike_align"),
    shell:
        """
        bowtie2 -x {params.idx_prefix} -p {threads} \
            -1 {input.unmapped1} -2 {input.unmapped2} \
            2> {log} \
        | samtools view -bS - \
        | samtools sort -@ {threads} -O BAM \
            -o {output.bam} -
        samtools index -@ {threads} {output.bam} {output.bai} >> {log} 2>&1
        samtools idxstats {output.bam} > {output.idxstats} 2>> {log}
        samtools flagstat -@ {threads} {output.bam} > {output.flagstat} 2>> {log}
        """


rule spike_summary:
    input:
        expand(R("5.QC/spike_in/{sample}_flagstat.txt"), sample=SAMPLES),
    output:
        tsv=R("5.QC/spike_in/Spikein_summary.tsv"),
        mqc=R("5.QC/spike_in/Spikein_summary_mqc.tsv"),
    params:
        script=os.path.join(WORKFLOW_DIR, "scripts", "spikein_summary.py"),
        python=PYTHON_BIN,
        samples=" ".join(SAMPLES),
        spike_dir=R("5.QC/spike_in"),
        name=SPIKE_IN["name"],
    log:
        R("logs/spike_in/summary.log"),
    threads: rthreads("spike_summary")
    resources:
        mem_mb=rmem("spike_summary"),
        runtime_min=rruntime("spike_summary"),
        runtime_sec=rruntime_sec("spike_summary"),
    shell:
        """
        {params.python} {params.script} \
            --samples {params.samples} --spike-dir {params.spike_dir} \
            --name {params.name} \
            --out {output.tsv} --mqc {output.mqc} > {log} 2>&1
        """
