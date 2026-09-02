###############################################
## 上游定量流程入口（trim → STAR → featureCounts → 表达矩阵）
## 运行：bash run.sh upstream <project_dir> [config.yaml] [jobs]
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
        expand("4.expression/{sample}.count", sample=SAMPLES),
        "4.expression/count.matrix.tsv",
        "4.expression/GeneExpression_TPM.xls",


include: os.path.join(workflow.basedir, "rules", "RNA-seq_upstream.smk")
