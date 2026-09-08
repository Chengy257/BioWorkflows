# Optional differential-expression stage (docs/TODO.md item 1, ported from
# the rna-seq workflow): DESeq2 over the per-class count matrices
# (4.expression/{klass}/{klass}_counts.tsv), driven by the optional
# group/batch columns of the sample table and the deg: config section.
# Parse-time guarded by _DEG_CLASSES from rules/common.smk: with the default
# deg.enabled=false this module declares nothing and the DAG is unchanged.
# Depends on workflow/rules/common.smk for R(), SCRIPTS, RSCRIPT, _SAMPLE_TABLE,
# _DEG_CONFIG/_DEG_CLASSES, and the rthreads/rmem/rruntime helpers.

if _DEG_CLASSES:
    rule deg_deseq2:
        input:
            counts=R("4.expression/{klass}/{klass}_counts.tsv"),
            sampleinfo=_SAMPLE_TABLE,
        output:
            R("6.DEG/{klass}/flag.log"),
        log:
            R("logs/deg_deseq2/{klass}.log"),
        params:
            rscript=RSCRIPT,
            script=os.path.join(SCRIPTS, "run_deseq2.R"),
            batch=str(_DEG_CONFIG.get("batch_correction", "F")),
            foldchange=_DEG_CONFIG.get("foldchange", 2),
            padj=_DEG_CONFIG.get("padj", 0.05),
            control_group=_DEG_CONFIG.get("control_group", "control"),
            pca_ntop=_DEG_CONFIG.get("pca_ntop", 2000),
            outdir=lambda wc, output: os.path.dirname(output[0]),
        threads:
            rthreads("deg_deseq2")
        resources:
            mem_mb=rmem("deg_deseq2"),
            runtime_min=rruntime("deg_deseq2"),
            runtime_sec=rruntime_sec("deg_deseq2"),
        shell:
            """
            {params.rscript} {params.script} \\
                -c {input.counts} -s {input.sampleinfo} \\
                -f {params.foldchange} -p {params.padj} -b {params.batch} -r {params.control_group} \\
                -n {params.pca_ntop} -k {wildcards.klass} \\
                -o {params.outdir} -t {threads} >> {log} 2>&1
            echo `date` " : All done!" > {output}
            """
