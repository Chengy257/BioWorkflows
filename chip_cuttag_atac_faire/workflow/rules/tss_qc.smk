# TSS enrichment for ATAC/FAIRE treat samples (qc.tss; included only when the
# switch is on and treat samples exist). ENCODE-flavored recipe: per-sample
# CPM coverage (binSize 10) over a strand-aware 1-bp TSS set derived from the
# configured gene-model BED, deeptools reference-point matrix +/- 2kb in 10bp
# bins, then the enrichment score (max of the baseline-normalized profile)
# computed by scripts/tss_score.py.

rule tss_bed:
    input:
        config["bed"],
    output:
        R("5.QC/tss/tss.bed"),
    params:
        script=os.path.join(WORKFLOW_DIR, "scripts", "tss_from_bed.py"),
        python=PYTHON_BIN,
    log:
        R("logs/tss/tss_bed.log"),
    threads: rthreads("tss_bed")
    resources:
        mem_mb=rmem("tss_bed"),
        runtime_min=rruntime("tss_bed"),
        runtime_sec=rruntime_sec("tss_bed"),
    shell:
        """
        {params.python} {params.script} {input} {output} > {log} 2>&1
        """


rule tss_matrix:
    input:
        bam=lambda wc: sample_bam(wc.sample),
        tss=R("5.QC/tss/tss.bed"),
    output:
        matrix=R("5.QC/tss/{sample}_matrix.gz"),
        png=R("5.QC/tss/{sample}_tss_profile.png"),
        score=R("5.QC/tss/{sample}_TSSE.txt"),
    wildcard_constraints:
        sample=_group_regex(TSS_SAMPLES),
    params:
        # the temporary CPM bigWig is an unmanaged intermediate, removed at the
        # end of the job (per-sample normalized tracks are a separate feature)
        bw=lambda wc, output: str(output.matrix)[:-len("_matrix.gz")] + "_cpm.bw",
        score_script=os.path.join(WORKFLOW_DIR, "scripts", "tss_score.py"),
        python=PYTHON_BIN,
    log:
        R("logs/tss/{sample}.log"),
    threads: rthreads("tss_matrix")
    resources:
        mem_mb=rmem("tss_matrix"),
        runtime_min=rruntime("tss_matrix"),
        runtime_sec=rruntime_sec("tss_matrix"),
    shell:
        """
        # ENCODE window: +/-2kb around each TSS in 10bp bins (400 bins);
        # CPM coverage at binSize 10 approximates the cut-site pileup
        # definition well enough for a QC curve
        bamCoverage -b {input.bam} --normalizeUsing CPM --binSize 10 \
            -p {threads} -o {params.bw} > {log} 2>&1
        computeMatrix reference-point -S {params.bw} -R {input.tss} \
            --referencePoint TSS -b 2000 -a 2000 --binSize 10 \
            --missingDataAsZero -p {threads} -o {output.matrix} >> {log} 2>&1
        plotProfile -m {output.matrix} -o {output.png} \
            --plotTitle "TSS profile {wildcards.sample}" >> {log} 2>&1
        {params.python} {params.score_script} {output.matrix} {wildcards.sample} {output.score} >> {log} 2>&1
        rm -f {params.bw}
        """


rule tss_summary:
    input:
        expand(R("5.QC/tss/{sample}_TSSE.txt"), sample=TSS_SAMPLES),
    output:
        tsv=R("5.QC/tss/TSSE_summary.tsv"),
        mqc=R("5.QC/tss/TSSE_summary_mqc.tsv"),
    log:
        R("logs/tss/summary.log"),
    threads: rthreads("tss_summary")
    resources:
        mem_mb=rmem("tss_summary"),
        runtime_min=rruntime("tss_summary"),
        runtime_sec=rruntime_sec("tss_summary"),
    shell:
        """
        printf 'sample\\tTSSE\\tcenter_norm\\tbaseline\\ttss_rows\\n' > {output.tsv}
        for f in {input}; do
            cat "$f" >> {output.tsv}
        done
        # MultiQC custom content (the _mqc.tsv convention)
        {{
            echo "# id: 'tss_table'"
            echo "# section_name: 'TSS enrichment (ATAC/FAIRE)'"
            echo "# format: 'tsv'"
            echo "# plot_type: 'table'"
            echo "# pconfig: {{'id': 'tss_table', 'title': 'TSS enrichment'}}"
            cat {output.tsv}
        }} > {output.mqc} 2> {log}
        """
