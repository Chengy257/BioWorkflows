
configfile: "/home/chengyu/workflows/snakemake/rna-seq-workflow/config_basic_defaulted.yaml"
#configfile: "/home/chengyu/workflows/snakemake/rna-seq-workflow/config_user_defined.yaml"

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
        "2.cleandata/trim/fastqc/multiqc_report.html",
        "3.align/mapping_stat.xls",
        expand("3.align/{sample}_Aligned.sortedByCoord.out.bam",sample=SAMPLES),
        expand("4.assembly/4.2.IsoformExpr/{sample}.tab",sample=SAMPLES),
        expand("4.assembly/4.1.Assembly_stringtie/{sample}.gtf",sample=SAMPLES),

include: "/home/chengyu/workflows/snakemake/rna-seq-workflow/rules/RNA-seq_upstream_AS.smk"
include: "/home/chengyu/workflows/snakemake/rna-seq-workflow/rules/RNA-seq_downstream.smk"
#include: "/home/chengyu/workflows/snakemake/rna-seq-workflow/rules/RNA-seq_alterSplice.smk"
