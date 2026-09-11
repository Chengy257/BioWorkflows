# Organelle read fraction (qc.organelle; included only when the switch is on):
# samtools idxstats per sample, then a project-wide summary of the
# chloroplast/mitochondrial mapped-read fraction (the plant analog of the
# mitochondrial-fraction check in ATAC-seq pipelines). Contig-name patterns
# come from qc.organelle_patterns (empty matches report a zero fraction).

rule organelle_idxstats:
    input:
        lambda wc: sample_bam(wc.sample),
    output:
        R("5.QC/organelle/{sample}_idxstats.tsv"),
    wildcard_constraints:
        sample=_group_regex(SAMPLES),
    log:
        R("logs/organelle/{sample}.log"),
    threads: rthreads("organelle_idxstats")
    resources:
        mem_mb=rmem("organelle_idxstats"),
        runtime_min=rruntime("organelle_idxstats"),
        runtime_sec=rruntime_sec("organelle_idxstats"),
    shell:
        """
        samtools idxstats {input} > {output} 2> {log}
        """


rule organelle_summary:
    input:
        expand(R("5.QC/organelle/{sample}_idxstats.tsv"), sample=SAMPLES),
    output:
        tsv=R("5.QC/organelle/Organelle_summary.tsv"),
        mqc=R("5.QC/organelle/Organelle_summary_mqc.tsv"),
    params:
        script=os.path.join(WORKFLOW_DIR, "scripts", "organelle_summary.py"),
        python=PYTHON_BIN,
        patterns=",".join(ORGANELLE_PATTERNS),
    log:
        R("logs/organelle/summary.log"),
    threads: rthreads("organelle_summary")
    resources:
        mem_mb=rmem("organelle_summary"),
        runtime_min=rruntime("organelle_summary"),
        runtime_sec=rruntime_sec("organelle_summary"),
    shell:
        """
        {params.python} {params.script} \
            --patterns {params.patterns} --out {output.tsv} --mqc {output.mqc} \
            {input} > {log} 2>&1
        """
