###############################################
## 下游分析规则：DESeq2 差异 → GO/KEGG 富集 → 组间 DEG 组合比较
###############################################
import os

SCRIPTS = os.path.join(workflow.basedir, "scripts")
ENVS = os.path.join(workflow.basedir, "envs")

rule runDESeq2:
    input:
        count="4.expression/count.matrix.tsv",
        sampleinfo=config["SampleListFile"],
    output:
        "5.DEG/flag.log",
    log:
        "logs/DESeq2/runDESeq2.log.txt",
    conda:
        os.path.join(ENVS, "deseq2.yaml")
    shell:
        """
        Rscript {SCRIPTS}/runDESeq2_multiGroup_230705.R \\
            -c {input.count} -s {input.sampleinfo} \\
            -f {config[FoldChange]} -p {config[padj]} -b F -r {config[control_group]} \\
            -o 5.DEG -t {config[threads]} >> {log} 2>&1
        echo `date` " : All done!" > 5.DEG/flag.log
        """


rule GO_KEGG_enrich:
    input:
        "5.DEG/flag.log",
    output:
        "5.DEG/GO_KEGG_enrich/flag.log",
    log:
        "logs/GO_KEGG_enrich/GO_KEGG_enrich.log.txt",
    conda:
        os.path.join(ENVS, "enrich.yaml")
    shell:
        """
        bash {SCRIPTS}/enrich.sh 5.DEG {config[species]} {config[annotation_tsv]} {config[orgdb_tarball]} >> {log} 2>&1
        """


rule DEGgroupCompare:
    input:
        "5.DEG/GO_KEGG_enrich/flag.log",
        sampleinfo=config["SampleListFile"],
    output:
        "6.DEGcompare/flag.log"
    log:
        "logs/DEGgroupCompare/DEGgroupCompare.log.txt"
    conda:
        os.path.join(ENVS, "enrich.yaml")
    shell:
        """
        bash {SCRIPTS}/DEGgroupCompare.sh {config[control_group]} {config[annotation_tsv]} {input.sampleinfo} {config[threads]} >> {log} 2>&1
        """
