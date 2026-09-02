###############################################
## 上游（.fq.gz 原始命名约定）+ 转录本组装/异构体表达
## 与 RNA-seq_upstream.smk 的差异：原始 fastq 命名为 {sample}.fq.gz / {sample}_1.fq.gz，
## trim_galore 附带 --fastqc，STAR 使用更严格的剪接参数
###############################################
import os

SCRIPTS = os.path.join(workflow.basedir, "scripts")
ENVS = os.path.join(workflow.basedir, "envs")


def _raw_reads(sample):
    """探测样本原始 fastq（精确文件名匹配，避免 id 互为前缀时 glob 误配）。"""
    for pat1, pat2 in (("1.rawdata/{0}_1.fq.gz", "1.rawdata/{0}_2.fq.gz"),
                       ("1.rawdata/{0}_1.fastq.gz", "1.rawdata/{0}_2.fastq.gz")):
        if os.path.exists(pat1.format(sample)) and os.path.exists(pat2.format(sample)):
            return "PE", [pat1.format(sample), pat2.format(sample)]
    for pat in ("1.rawdata/{0}.fq.gz", "1.rawdata/{0}.fastq.gz"):
        if os.path.exists(pat.format(sample)):
            return "SE", [pat.format(sample)]
    raise ValueError(
        "sample {0}: 1.rawdata/ 下未找到原始 fastq"
        "（支持 {0}_1.fq.gz+{0}_2.fq.gz、{0}_1.fastq.gz+{0}_2.fastq.gz PE，"
        "或 {0}.fq.gz / {0}.fastq.gz SE）".format(sample))


def get_fastq(wildcards):
    layout, _ = _raw_reads(wildcards.sample)
    s = wildcards.sample
    if layout == "PE":
        return [f"2.cleandata/trim/{s}_1_val_1.fq.gz",
                f"2.cleandata/trim/{s}_2_val_2.fq.gz"]
    return [f"2.cleandata/trim/{s}_trimmed.fq.gz"]


def trimmed_reads(sample):
    layout, _ = _raw_reads(sample)
    if layout == "PE":
        return f"2.cleandata/trim/{sample}_1_val_1.fq.gz 2.cleandata/trim/{sample}_2_val_2.fq.gz"
    return f"2.cleandata/trim/{sample}_trimmed.fq.gz"


rule trimAdapter_SE:
    input:
        fq="1.rawdata/{sample}.fq.gz",
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
        trim_galore -q 30 --stringency 3 -e 0.1 --gzip -j {config[threads]} -o 2.cleandata/trim/ {input.fq} \\
            --fastqc --fastqc_args "--outdir 2.cleandata/trim/fastqc " >> {log} 2>&1
        touch 2.cleandata/trim/fastqc/{wildcards.sample}.fastqc.flag
        """


rule trimAdapter_PE:
    input:
        fq1="1.rawdata/{sample}_1.fq.gz",
        fq2="1.rawdata/{sample}_2.fq.gz",
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
        trim_galore -q 30 --stringency 3 -e 0.1 --gzip -j {config[threads]} -o 2.cleandata/trim/ --paired {input.fq1} {input.fq2} \\
            --fastqc --fastqc_args "--outdir 2.cleandata/trim/fastqc " >> {log} 2>&1
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
        star=(" --twopassMode Basic --outFilterType BySJout --alignIntronMin 20 --alignIntronMax 5000 "
              "--alignMatesGapMax 5000 --outFilterMatchNminOverLread 0.66 --outFilterScoreMinOverLread 0.66 "
              "--winAnchorMultimapNmax 70 --seedSearchStartLmax 45 --outSAMattrIHstart 0 "
              "--outSAMstrandField intronMotif --genomeLoad NoSharedMemory --quantMode TranscriptomeSAM GeneCounts "),
        reads=lambda wc: trimmed_reads(wc.sample),
    conda:
        os.path.join(ENVS, "align.yaml")
    shell:
        """
        STAR {params.star} \\
             --runThreadN {config[threads]} --genomeDir 0.index/star_genome/ \\
             --readFilesCommand zcat --outSAMtype BAM Unsorted --outSAMattributes All \\
             --readFilesIn {params.reads} \\
             --outFileNamePrefix 3.align/{wildcards.sample}_ >> {log} 2>&1
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
        ## 逐样本解析（同 RNA-seq_upstream.smk），修复行错位与 id 截断
        echo -e "ID\\tTotal_Reads\\tClean_Reads\\tUniquely_mapped\\tUniquely_mapped_ratio" > 3.align/mapping_stat.xls
        for log in 3.align/*_Log.final.out; do
            [ -e "$log" ] || continue
            sid=$(basename "$log" _Log.final.out)
            input_reads=$(fgrep "Number of input reads" "$log" | awk '{{print $NF}}')
            unique_num=$(fgrep "Uniquely mapped reads number" "$log" | awk '{{print $NF}}')
            unique_pct=$(fgrep "Uniquely mapped reads %" "$log" | awk '{{print $NF}}')
            total="NA"
            for tr in "2.cleandata/trim/${{sid}}_1.fq.gz_trimming_report.txt" \\
                      "2.cleandata/trim/${{sid}}_1.fastq.gz_trimming_report.txt" \\
                      "2.cleandata/trim/${{sid}}.fq.gz_trimming_report.txt" \\
                      "2.cleandata/trim/${{sid}}.fastq.gz_trimming_report.txt"; do
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


rule Assemble:
    input:
        bam="3.align/{sample}_Aligned.sortedByCoord.out.bam",
        ref=config["gtf"],
        strandedness="3.align/{sample}.strandedness",
    output:
        "4.assembly/4.1.Assembly_stringtie/{sample}.gtf",
    conda:
        os.path.join(ENVS, "assembly.yaml")
    log:
        "logs/Assemble/{sample}_stringtie.log.txt"
    shell:
        """
        strandedness=$(head -1 {input.strandedness} | awk '{{print $1}}')
        ## 链型选项：firststrand → --rf，secondstrand → --fr
        if [ "$strandedness" == "firststrand" ]; then
            stringtie -p {config[threads]} --rf -o {output} -G {input.ref} {input.bam} >> {log} 2>&1
        elif [ "$strandedness" == "secondstrand" ]; then
            stringtie -p {config[threads]} --fr -o {output} -G {input.ref} {input.bam} >> {log} 2>&1
        else
            stringtie -p {config[threads]} -o {output} -G {input.ref} {input.bam} >> {log} 2>&1
        fi
        """


rule gtf_merge:
    input:
        gtf=expand("4.assembly/4.1.Assembly_stringtie/{sample}.gtf", sample=SAMPLES),
        ref=config["gtf"],
        genome=config["genome"],
    output:
        "4.assembly/4.1.Assembly_stringtie/merged.gtf",
    conda:
        os.path.join(ENVS, "assembly.yaml")
    log:
        "logs/gtf_merge.log.txt"
    shell:
        """
        ## 载入辅助函数（runGFFcompare / getFasta）
        source {SCRIPTS}/lncRNA_functions.sh
        ## 合并各样本 gtf 并与参考注释比较
        stringtie --merge -G {input.ref} -i -o 4.assembly/4.1.Assembly_stringtie/merged.gtf {input.gtf} >> {log} 2>&1
        runGFFcompare {input.ref} 4.assembly/4.1.Assembly_stringtie/merged.gtf {input.genome} >> {log} 2>&1
        """


rule isoform_expr:
    input:
        bam="3.align/{sample}_Aligned.sortedByCoord.out.bam",
        ref="4.assembly/4.1.Assembly_stringtie/merged.gtf",
        strandedness="3.align/{sample}.strandedness",
    output:
        gtf="4.assembly/4.2.IsoformExpr/{sample}.gtf",
        tab="4.assembly/4.2.IsoformExpr/{sample}.tab",
    conda:
        os.path.join(ENVS, "assembly.yaml")
    log:
        "logs/isoform_expr/{sample}_stringtie.log.txt"
    shell:
        """
        strandedness=$(head -1 {input.strandedness} | awk '{{print $1}}')
        if [ "$strandedness" == "firststrand" ]; then
            stringtie -p {config[threads]} --rf -o {output.gtf} -e -G {input.ref} {input.bam} >> {log} 2>&1
        elif [ "$strandedness" == "secondstrand" ]; then
            stringtie -p {config[threads]} --fr -o {output.gtf} -e -G {input.ref} {input.bam} >> {log} 2>&1
        else
            stringtie -p {config[threads]} -o {output.gtf} -e -G {input.ref} {input.bam} >> {log} 2>&1
        fi
        ## 提取转录本表达
        cat {output.gtf} | grep -v "^#" | awk -v OFS="\\t" 'BEGIN{{print "transcript","FPKM","TPM"}} {{if($3=="transcript"){{print $12,$(NF-2),$NF}}}}' | sed 's/[;"]//g' > {output.tab}
        """
