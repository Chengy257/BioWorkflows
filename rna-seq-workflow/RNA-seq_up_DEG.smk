###############################################
## 上游定量 + 差异分析 + 功能富集 + 组合比较
## 运行：bash run.sh deg <project_dir> [config.yaml] [jobs]
###############################################
import os

configfile: os.path.join(workflow.basedir, "config_basic_defaulted.yaml")
## 如需覆盖配置：复制 config_user_defined.yaml 模板修改后，以 --configfile 传入（命令行优先级更高）


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
        "2.cleandata/trim/fastqc/multiqc_report.html",
        "3.align/mapping_stat.xls",
        expand("3.align/{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        "4.expression/count.matrix.tsv",
        "4.expression/GeneExpression_TPM.xls",
        "5.DEG/GO_KEGG_enrich/flag.log",
        "5.DEG/flag.log",
        "6.DEGcompare/flag.log",


include: os.path.join(workflow.basedir, "rules", "RNA-seq_upstream.smk")
include: os.path.join(workflow.basedir, "rules", "RNA-seq_downstream.smk")
