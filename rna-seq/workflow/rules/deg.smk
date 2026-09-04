###############################################
## 差异分析规则（仅 pipeline=deg）：DESeq2 → GO/KEGG 富集 + GSEA → 组间比较
###############################################

rule runDESeq2:
    input:
        count=R("4.expression/count.matrix.tsv"),
        sampleinfo=config["SampleListFile"],
    output:
        R("5.DEG/flag.log"),
    log:
        R("logs/DESeq2/runDESeq2.log.txt"),
    conda:
        os.path.join(ENVS, "deseq2.yaml")
    shell:
        """
        Rscript {SCRIPTS}/run_deseq2.R \\
            -c {input.count} -s {input.sampleinfo} \\
            -f {config[FoldChange]} -p {config[padj]} -b F -r {config[control_group]} \\
            -n {config[pca_ntop]} \\
            -o {RD}5.DEG -t {config[threads]} >> {log} 2>&1
        echo `date` " : All done!" > {output}
        """


rule GO_KEGG_enrich:
    input:
        R("5.DEG/flag.log"),
    output:
        R("5.DEG/GO_KEGG_enrich/flag.log"),
    log:
        R("logs/GO_KEGG_enrich/GO_KEGG_enrich.log.txt"),
    conda:
        os.path.join(ENVS, "enrich.yaml")
    shell:
        """
        bash {SCRIPTS}/enrich.sh {RD}5.DEG {config[species]} {config[annotation_tsv]} {config[orgdb_tarball]} {config[kegg_organism]} >> {log} 2>&1
        """


rule DEGgroupCompare:
    input:
        R("5.DEG/GO_KEGG_enrich/flag.log"),
        sampleinfo=config["SampleListFile"],
    output:
        R("6.DEGcompare/flag.log"),
    log:
        R("logs/DEGgroupCompare/DEGgroupCompare.log.txt"),
    conda:
        os.path.join(ENVS, "enrich.yaml")
    shell:
        """
        bash {SCRIPTS}/DEGgroupCompare.sh {RD} {config[control_group]} {config[annotation_tsv]} {input.sampleinfo} {config[threads]} >> {log} 2>&1
        """
