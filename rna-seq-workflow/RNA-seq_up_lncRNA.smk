configfile: "/home/chengyu/workflows/snakemake/rna-seq-workflow/config_basic_defaulted.yaml"
configfile: "/home/chengyu/workflows/snakemake/rna-seq-workflow/config_lncRNA.yaml"

## function define: get sample IDs
def get_samples():
    ids = []
    with open(config["SampleListFile"], "r") as samples_list:
        next(samples_list)
        for line in samples_list:
            line = line.strip().split(",")
            ids.append(line[0])
    samples_list.close()
    return ids
SAMPLES = get_samples()

rule results:
    input:
        # "2.cleandata/fastqc/multiqc/multiqc_report.html",
        "4.expression/GeneExpression_TPM.xls",
        # "2.cleandata/trim/fastqc/multiqc_report.html",
        "3.align/mapping_stat.xls",
        expand("3.align/{sample}_Aligned.sortedByCoord.out.bam",sample=SAMPLES),
        "4.LncRNA/4.5.Final_lncRNA/final_lncRNA.fa",
        "5.expression/lncRNA/GeneExpression_TPM.xls",

include: "/home/chengyu/workflows/snakemake/rna-seq-workflow/rules/RNA-seq_upstream.smk"
# include: "/home/chengyu/workflows/snakemake/rna-seq-workflow/RNA-seq_downstream.smk"
include: "/home/chengyu/workflows/snakemake/rna-seq-workflow/rules/RNA-seq_lncRNA_DenovoIdenti.smk"
