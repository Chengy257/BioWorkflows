###############################################
## 上游定量规则：trim → FastQC/MultiQC → STAR 比对 → 链型判断 → featureCounts → 表达矩阵
## 由各入口 Snakefile include；SAMPLES 与 config 由入口文件先行定义
###############################################
import os

SCRIPTS = os.path.join(workflow.basedir, "scripts")
ENVS = os.path.join(workflow.basedir, "envs")


def _raw_reads(sample):
    """探测样本原始 fastq：先 PE（<id>_1/<id>_2），后 SE（<id>.fastq.gz / <id>.fq.gz）。
    使用精确文件名匹配，避免样本 id 互为前缀时 glob 误配（如 ctal_1 与 ctal_1_1）。"""
    for pat1, pat2 in (("1.rawdata/{0}_1.fastq.gz", "1.rawdata/{0}_2.fastq.gz"),
                       ("1.rawdata/{0}_1.fq.gz", "1.rawdata/{0}_2.fq.gz")):
        if os.path.exists(pat1.format(sample)) and os.path.exists(pat2.format(sample)):
            return "PE", [pat1.format(sample), pat2.format(sample)]
    for pat in ("1.rawdata/{0}.fastq.gz", "1.rawdata/{0}.fq.gz"):
        if os.path.exists(pat.format(sample)):
            return "SE", [pat.format(sample)]
    raise ValueError(
        "sample {0}: 1.rawdata/ 下未找到原始 fastq"
        "（支持 {0}_1.fastq.gz+{0}_2.fastq.gz、{0}_1.fq.gz+{0}_2.fq.gz PE，"
        "或 {0}.fastq.gz / {0}.fq.gz SE）".format(sample))


def get_fastq(wildcards):
    """runSTAR 的输入依赖（trim 后文件），随文库类型自动切换。"""
    layout, _ = _raw_reads(wildcards.sample)
    s = wildcards.sample
    if layout == "PE":
        return [f"2.cleandata/trim/{s}_1_val_1.fq.gz",
                f"2.cleandata/trim/{s}_2_val_2.fq.gz"]
    return [f"2.cleandata/trim/{s}_trimmed.fq.gz"]


def trimmed_reads(sample):
    """STAR --readFilesIn 参数：PE 两个文件，SE 一个文件。"""
    layout, _ = _raw_reads(sample)
    if layout == "PE":
        return f"2.cleandata/trim/{sample}_1_val_1.fq.gz 2.cleandata/trim/{sample}_2_val_2.fq.gz"
    return f"2.cleandata/trim/{sample}_trimmed.fq.gz"


def is_paired_end(wildcards):
    return "True" if _raw_reads(wildcards.sample)[0] == "PE" else "False"


rule trimAdapter_SE:
    input:
        fq="1.rawdata/{sample}.fastq.gz",
    output:
        "2.cleandata/trim/{sample}_trimmed.fq.gz",
        "2.cleandata/trim/fastqc/{sample}.fastqc.flag",
    log:
        "logs/trim/{sample}_log.txt",
    conda:
        os.path.join(ENVS, "qc.yaml")
    shell:
        """
        mkdir -p 2.cleandata/trim/fastqc
        trim_galore -q 30 --stringency 3 -e 0.1 --gzip -j {config[threads]} -o 2.cleandata/trim/ {input.fq} >> {log} 2>&1
        touch 2.cleandata/trim/fastqc/{wildcards.sample}.fastqc.flag
        """


rule trimAdapter_PE:
    input:
        fq1="1.rawdata/{sample}_1.fastq.gz",
        fq2="1.rawdata/{sample}_2.fastq.gz",
    output:
        "2.cleandata/trim/{sample}_1_val_1.fq.gz",
        "2.cleandata/trim/{sample}_2_val_2.fq.gz",
        "2.cleandata/trim/fastqc/{sample}.fastqc.flag",
    log:
        "logs/trim/{sample}_log.txt",
    conda:
        os.path.join(ENVS, "qc.yaml")
    shell:
        """
        mkdir -p 2.cleandata/trim/fastqc
        trim_galore -q 30 --stringency 3 -e 0.1 --gzip -j {config[threads]} -o 2.cleandata/trim/ --paired {input.fq1} {input.fq2} >> {log} 2>&1
        touch 2.cleandata/trim/fastqc/{wildcards.sample}.fastqc.flag
        """


rule multiQC_cleaned:
    input:
        expand("2.cleandata/trim/fastqc/{sample}.fastqc.flag", sample=SAMPLES)
    output:
        "2.cleandata/trim/fastqc/multiqc_report.html",
    log:
        "logs/fastqc/multiqc_log.txt"
    conda:
        os.path.join(ENVS, "qc.yaml")
    shell:
        "multiqc --force 2.cleandata/trim/fastqc/ -o 2.cleandata/trim/fastqc/ >> {log} 2>&1"


rule STAR_index:
    input:
        GENOME=config["genome"],
        GTF=config["gtf"],
    output:
        "0.index/star_genome/SAindex",
        star_index=directory("0.index/star_genome/"),
    log:
        "logs/index/star_index_log.txt",
    conda:
        os.path.join(ENVS, "align.yaml")
    shell:
        """
        mkdir -p 0.index/star_genome
        STAR --runThreadN {config[threads]} --runMode genomeGenerate --genomeDir {output.star_index} \\
             --genomeFastaFiles {input.GENOME} --sjdbGTFfile {input.GTF} >> {log} 2>&1
        """


rule runSTAR:
    input:
        unpack(get_fastq),
        index="0.index/star_genome/SAindex",
    output:
        bam="3.align/{sample}_Aligned.sortedByCoord.out.bam",
        bai="3.align/{sample}_Aligned.sortedByCoord.out.bam.bai",
    log:
        "logs/runSTAR/{sample}.log.txt",
    params:
        reads=lambda wc: trimmed_reads(wc.sample),
    conda:
        os.path.join(ENVS, "align.yaml")
    shell:
        """
        STAR --runThreadN {config[threads]} --twopassMode Basic \\
             --genomeLoad NoSharedMemory --genomeDir 0.index/star_genome/ \\
             --readFilesCommand zcat --outSAMtype BAM Unsorted --outSAMattributes All \\
             --quantMode GeneCounts --readFilesIn {params.reads} \\
             --outFileNamePrefix 3.align/{wildcards.sample}_ --outSAMattrIHstart 0 >> {log} 2>&1
        samtools sort -@ {config[threads]} -O BAM -o {output.bam} 3.align/{wildcards.sample}_Aligned.out.bam && rm 3.align/{wildcards.sample}_Aligned.out.bam >> {log} 2>&1
        samtools index -@ {config[threads]} {output.bam} >> {log} 2>&1
        """


rule Mapping_stat:
    input:
        expand("3.align/{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES),
    output:
        "3.align/mapping_stat.xls",
    log:
        "logs/Mapping_stat.log.txt"
    shell:
        """
        ## 逐样本解析：样本 id 取自日志文件名（保留下划线），trim 报告按精确文件名匹配，
        ## 修复旧版 cat 全部报告导致的行错位与 cut -d_ -f1 导致的 id 截断
        echo -e "ID\\tTotal_Reads\\tClean_Reads\\tUniquely_mapped\\tUniquely_mapped_ratio" > 3.align/mapping_stat.xls
        for log in 3.align/*_Log.final.out; do
            [ -e "$log" ] || continue
            sid=$(basename "$log" _Log.final.out)
            input_reads=$(fgrep "Number of input reads" "$log" | awk '{{print $NF}}')
            unique_num=$(fgrep "Uniquely mapped reads number" "$log" | awk '{{print $NF}}')
            unique_pct=$(fgrep "Uniquely mapped reads %" "$log" | awk '{{print $NF}}')
            total="NA"
            for tr in "2.cleandata/trim/${{sid}}_1.fastq.gz_trimming_report.txt" \\
                      "2.cleandata/trim/${{sid}}_1.fq.gz_trimming_report.txt" \\
                      "2.cleandata/trim/${{sid}}.fastq.gz_trimming_report.txt" \\
                      "2.cleandata/trim/${{sid}}.fq.gz_trimming_report.txt"; do
                if [ -f "$tr" ]; then
                    total=$(fgrep "Total reads processed:" "$tr" | head -1 | awk '{{print $NF}}' | sed 's/,//g')
                    break
                fi
            done
            printf '%s\\t%s\\t%s\\t%s\\t%s\\n' "$sid" "$total" "$input_reads" "$unique_num" "$unique_pct" >> 3.align/mapping_stat.xls
        done
        """


rule check_strandedness:
    input:
        bam="3.align/{sample}_Aligned.sortedByCoord.out.bam",
        bed=config["bed"],
    output:
        strand="3.align/{sample}.strandedness",
        infer="3.align/{sample}_infer_experiment.out",
    log:
        "logs/check_strandedness/{sample}_check_strandedness.log.txt",
    conda:
        os.path.join(ENVS, "align.yaml")
    shell:
        """
        infer_experiment.py -r {input.bed} -i {input.bam} > {output.infer} 2> {log}
        tail -2 {output.infer} | awk '{{print $NF}}' | \\
            awk '{{f1=$0;getline;f2=$0; \\
                if(f1-f2 > 0.4) print "secondstrand"; \\
                else if(f2-f1 > 0.4) print "firststrand"; \\
                else print "unstrand"}}' > {output.strand}
        """


rule featureCount_R:
    input:
        bam="3.align/{sample}_Aligned.sortedByCoord.out.bam",
        gtf=config["gtf"],
        strand="3.align/{sample}.strandedness",
    output:
        count="4.expression/{sample}.count",
        stat="4.expression/{sample}.log",
    log:
        "logs/featureCount_R/{sample}.log.txt"
    params:
        is_pe=is_paired_end,
    conda:
        os.path.join(ENVS, "quant.yaml")
    shell:
        """
        strandedness=$(head -1 {input.strand} | awk '{{print $1}}')
        ## 链特异性定量：featureCounts strandSpecific（1=整合链，2=反转链）
        if [ "$strandedness" == "firststrand" ]; then
            strand="2"
        elif [ "$strandedness" == "secondstrand" ]; then
            strand="1"
        else
            strand="0"
        fi
        Rscript {SCRIPTS}/run-featurecounts.R -t {config[threads]} -b {input.bam} -g {input.gtf} \\
            -s $strand -i {params.is_pe} -o 4.expression/{wildcards.sample} >> {log} 2>&1
        """


rule count_merge:
    input:
        counts=expand("4.expression/{sample}.count", sample=SAMPLES),
        logs=expand("4.expression/{sample}.log", sample=SAMPLES),
    output:
        "4.expression/count.matrix.tsv",
        "4.expression/GeneExpression_TPM.xls",
        "4.expression/GeneExpression_FPKM.xls",
        "4.expression/GeneCount_Assigned_logs.xls",
    log:
        "logs/count_merge/log.txt"
    shell:
        "python {SCRIPTS}/merge_featurecounts.py 4.expression >> {log} 2>&1"
