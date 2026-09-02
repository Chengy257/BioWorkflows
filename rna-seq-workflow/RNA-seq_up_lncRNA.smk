###############################################
## 上游定量 + lncRNA de novo 鉴定 + lncRNA 表达定量
## 运行：bash run.sh lncrna <project_dir> [config.yaml] [jobs]
###############################################
import os

configfile: os.path.join(workflow.basedir, "config_basic_defaulted.yaml")
configfile: os.path.join(workflow.basedir, "config_lncRNA.yaml")   ## 同名键覆盖基础配置


def get_samples():
    ids = []
    with open(config["SampleListFile"], "r") as samples_list:
        next(samples_list)
        for line in samples_list:
            line = line.strip().split(",")
            ids.append(line[0])
    return ids


SAMPLES = get_samples()

rule results:
    input:
        "4.expression/GeneExpression_TPM.xls",
        "3.align/mapping_stat.xls",
        expand("3.align/{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        "4.LncRNA/4.5.Final_lncRNA/final_lncRNA.fa",
        "5.expression/lncRNA/GeneExpression_TPM.xls",


include: os.path.join(workflow.basedir, "rules", "RNA-seq_upstream.smk")
include: os.path.join(workflow.basedir, "rules", "RNA-seq_lncRNA_DenovoIdenti.smk")
