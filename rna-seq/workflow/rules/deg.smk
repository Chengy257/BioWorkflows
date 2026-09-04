###############################################
## 差异分析规则（仅 pipeline=deg）：DESeq2 → GO/KEGG 富集 + GSEA → 组间比较
###############################################

rule runDESeq2:
    input:
        count_matrix=R("4.expression/count.matrix.tsv"),
        sampleinfo=config["SampleListFile"],
    output:
        R("5.DEG/flag.log"),
    log:
        R("logs/DESeq2/runDESeq2.log.txt"),
    params:
        batch=lambda wc: config.get("batch_correction", "F"),
        script=os.path.join(SCRIPTS, "run_deseq2.R"),
        foldchange=config["FoldChange"],
        padj=config["padj"],
        control_group=config["control_group"],
        pca_ntop=config["pca_ntop"],
        outdir=lambda wc, output: os.path.dirname(output[0]),
        rscript=RSCRIPT,
    threads:
        rthreads("deseq2")
    resources:
        mem_mb=rmem("deseq2"),
        runtime_min=rruntime("deseq2"),
        runtime_sec=rruntime_sec("deseq2"),
    shell:
        """
        {params.rscript} {params.script} \\
            -c {input.count_matrix} -s {input.sampleinfo} \\
            -f {params.foldchange} -p {params.padj} -b {params.batch} -r {params.control_group} \\
            -n {params.pca_ntop} \\
            -o {params.outdir} -t {threads} >> {log} 2>&1
        echo `date` " : All done!" > {output}
        """


rule GO_KEGG_enrich:
    input:
        R("5.DEG/flag.log"),
    output:
        R("5.DEG/GO_KEGG_enrich/flag.log"),
    log:
        R("logs/GO_KEGG_enrich/GO_KEGG_enrich.log.txt"),
    params:
        script=os.path.join(SCRIPTS, "enrich.sh"),
        deg_dir=lambda wc, input: os.path.dirname(input[0]),
        species=config["species"],
        annotation_tsv=config["annotation_tsv"],
        kegg_organism=config["kegg_organism"],
    threads:
        rthreads("enrichment")
    resources:
        mem_mb=rmem("enrichment"),
        runtime_min=rruntime("enrichment"),
        runtime_sec=rruntime_sec("enrichment"),
    shell:
        """
        bash {params.script} {params.deg_dir} {params.species} {params.annotation_tsv} {params.kegg_organism} >> {log} 2>&1
        """


rule DEGgroupCompare:
    input:
        R("5.DEG/GO_KEGG_enrich/flag.log"),
        sampleinfo=config["SampleListFile"],
    output:
        R("6.DEGcompare/flag.log"),
    log:
        R("logs/DEGgroupCompare/DEGgroupCompare.log.txt"),
    params:
        script=os.path.join(SCRIPTS, "DEGgroupCompare.sh"),
        results_root=lambda wc, output: os.path.dirname(os.path.dirname(output[0])) + os.sep,
        control_group=config["control_group"],
        annotation_tsv=config["annotation_tsv"],
        species=config["species"],
    threads:
        rthreads("deg_compare")
    resources:
        mem_mb=rmem("deg_compare"),
        runtime_min=rruntime("deg_compare"),
        runtime_sec=rruntime_sec("deg_compare"),
    shell:
        """
        ## Species is forwarded to the comparison R jobs; OrgDb must already be installed.
        bash {params.script} {params.results_root} {params.control_group} {params.annotation_tsv} {input.sampleinfo} {threads} \\
            '{params.species}' >> {log} 2>&1
        """
