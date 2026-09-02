###############################################
## 上游定量 + 转录本组装/异构体表达（可变剪接分析基础）
## 运行：bash run.sh as <project_dir> [config.yaml] [jobs]
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
        expand("4.assembly/4.2.IsoformExpr/{sample}.tab", sample=SAMPLES),
        expand("4.assembly/4.1.Assembly_stringtie/{sample}.gtf", sample=SAMPLES),


include: os.path.join(workflow.basedir, "rules", "RNA-seq_upstream_AS.smk")
include: os.path.join(workflow.basedir, "rules", "RNA-seq_downstream.smk")
# include: os.path.join(workflow.basedir, "rules", "RNA-seq_alterSplice.smk")   ## 尚未实现，见 docs/优化路线图.md
