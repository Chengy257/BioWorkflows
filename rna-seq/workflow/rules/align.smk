###############################################
## 质控与比对规则（所有管线共用）：
##   trim（统一命名约定 + FastQC，修复 P1-1）→ MultiQC 全流程汇总（P2-9）
##   → STAR（SortedByCoordinate 直出，P1-4）→ 链特异性推断
## 依赖 common.smk 中的 SAMPLES / R() / _raw_reads / STAR_ARGS 等
###############################################

def raw_se(wildcards):
    """SE trim 规则的原始 fastq 输入（精确探测）。"""
    layout, files = _raw_reads(wildcards.sample)
    if layout != "SE":
        raise ValueError(f"样本 {wildcards.sample} 探测为 {layout} 文库，SE trim 规则被误触发")
    return files[0]


def raw_pe(wildcards):
    """PE trim 规则的原始 fastq 输入（精确探测）。"""
    layout, files = _raw_reads(wildcards.sample)
    if layout != "PE":
        raise ValueError(f"样本 {wildcards.sample} 探测为 {layout} 文库，PE trim 规则被误触发")
    return {"fq1": files[0], "fq2": files[1]}


rule trimAdapter_SE:
    input:
        fq=raw_se,
    output:
        fq_out=R("2.cleandata/trim/{sample}_trimmed.fq.gz"),
        report=R("2.cleandata/trim/{sample}_trimming_report.txt"),
        fastqc_html=R("2.cleandata/trim/fastqc/{sample}_trimmed_fastqc.html"),
        fastqc_zip=R("2.cleandata/trim/fastqc/{sample}_trimmed_fastqc.zip"),
    log:
        R("logs/trim/{sample}_log.txt"),
    conda:
        os.path.join(ENVS, "qc.yaml")
    shell:
        """
        mkdir -p {RD}2.cleandata/trim/fastqc
        trim_galore -q 30 --stringency 3 -e 0.1 --gzip -j {config[threads]} -o {RD}2.cleandata/trim/ {input.fq} \\
            --fastqc --fastqc_args "--outdir {RD}2.cleandata/trim/fastqc " >> {log} 2>&1
        ## trim 报告按原始文件名生成，统一重命名为以样本 id 为键的规范名
        mv "{RD}2.cleandata/trim/$(basename "{input.fq}")_trimming_report.txt" {output.report}
        """


rule trimAdapter_PE:
    input:
        unpack(raw_pe),
    output:
        fq1_out=R("2.cleandata/trim/{sample}_1_val_1.fq.gz"),
        fq2_out=R("2.cleandata/trim/{sample}_2_val_2.fq.gz"),
        report1=R("2.cleandata/trim/{sample}_1_trimming_report.txt"),
        report2=R("2.cleandata/trim/{sample}_2_trimming_report.txt"),
        fastqc1_html=R("2.cleandata/trim/fastqc/{sample}_1_val_1_fastqc.html"),
        fastqc1_zip=R("2.cleandata/trim/fastqc/{sample}_1_val_1_fastqc.zip"),
        fastqc2_html=R("2.cleandata/trim/fastqc/{sample}_2_val_2_fastqc.html"),
        fastqc2_zip=R("2.cleandata/trim/fastqc/{sample}_2_val_2_fastqc.zip"),
    log:
        R("logs/trim/{sample}_log.txt"),
    conda:
        os.path.join(ENVS, "qc.yaml")
    shell:
        """
        mkdir -p {RD}2.cleandata/trim/fastqc
        trim_galore -q 30 --stringency 3 -e 0.1 --gzip -j {config[threads]} -o {RD}2.cleandata/trim/ --paired {input.fq1} {input.fq2} \\
            --fastqc --fastqc_args "--outdir {RD}2.cleandata/trim/fastqc " >> {log} 2>&1
        mv "{RD}2.cleandata/trim/$(basename "{input.fq1}")_trimming_report.txt" {output.report1}
        mv "{RD}2.cleandata/trim/$(basename "{input.fq2}")_trimming_report.txt" {output.report2}
        """


def multiqc_inputs(wildcards=None):
    """全流程 MultiQC 输入：FastQC、trim 报告、STAR 日志、featureCounts 分配统计。"""
    files = []
    for s in SAMPLES:
        layout, _ = _raw_reads(s)
        if layout == "PE":
            files += [R(f"2.cleandata/trim/fastqc/{s}_1_val_1_fastqc.zip"),
                      R(f"2.cleandata/trim/fastqc/{s}_2_val_2_fastqc.zip")]
        else:
            files.append(R(f"2.cleandata/trim/fastqc/{s}_trimmed_fastqc.zip"))
        files += trim_reports(s)
        files.append(R(f"3.align/{s}_Log.final.out"))
    if PIPELINE in ("upstream", "deg", "lncrna"):
        files += [R(f"4.expression/{s}.log") for s in SAMPLES]
    if PIPELINE == "lncrna":
        files += [R(f"5.expression/lncRNA/{s}.log") for s in SAMPLES]
    return files


rule multiqc:
    input:
        multiqc_inputs,
    output:
        R("multiqc/multiqc_report.html"),
    log:
        R("logs/multiqc/multiqc_log.txt"),
    conda:
        os.path.join(ENVS, "qc.yaml")
    shell:
        """
        mkdir -p {RD}multiqc
        multiqc --force --config {WORKFLOW_DIR}/multiqc_config.yaml -o {RD}multiqc {RD} >> {log} 2>&1
        """


rule STAR_index:
    input:
        GENOME=res("genome"),
        GTF=res("gtf"),
    output:
        sa=R("0.index/star_genome/SAindex"),
        genome_sa=R("0.index/star_genome/genomeSA"),
        genome_sj=R("0.index/star_genome/genomeSJ"),
        chr_start=R("0.index/star_genome/chrStart.txt"),
    log:
        R("logs/index/star_index_log.txt"),
    conda:
        os.path.join(ENVS, "align.yaml")
    shell:
        """
        mkdir -p {RD}0.index/star_genome
        STAR --runThreadN {config[threads]} --runMode genomeGenerate --genomeDir {RD}0.index/star_genome/ \\
             --genomeFastaFiles {input.GENOME} --sjdbGTFfile {input.GTF} >> {log} 2>&1
        """


rule runSTAR:
    input:
        unpack(get_fastq),
        index=R("0.index/star_genome/SAindex"),
    output:
        bam=R("3.align/{sample}_Aligned.sortedByCoord.out.bam"),
        bai=R("3.align/{sample}_Aligned.sortedByCoord.out.bam.bai"),
        log_final=R("3.align/{sample}_Log.final.out"),
    log:
        R("logs/runSTAR/{sample}.log.txt"),
    params:
        star_args=STAR_ARGS,
        reads=lambda wc: trimmed_reads(wc.sample),
        prefix=lambda wc: R(f"3.align/{wc.sample}_"),
    conda:
        os.path.join(ENVS, "align.yaml")
    shell:
        """
        mkdir -p {RD}3.align
        STAR {params.star_args} \\
             --runThreadN {config[threads]} --genomeDir {RD}0.index/star_genome/ \\
             --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMattributes All \\
             --readFilesIn {params.reads} \\
             --outFileNamePrefix {params.prefix} >> {log} 2>&1
        samtools index -@ {config[threads]} {output.bam} >> {log} 2>&1
        """


rule Mapping_stat:
    input:
        bams=expand(R("3.align/{sample}_Aligned.sortedByCoord.out.bam"), sample=SAMPLES),
        reports=all_trim_reports,
    output:
        R("3.align/mapping_stat.xls"),
    log:
        R("logs/Mapping_stat.log.txt")
    shell:
        """
        ## 逐样本解析：样本 id 取自日志文件名（保留下划线）；trim 总数取规范报告的 R1 读数（PE 即读对数）
        echo -e "ID\\tTotal_Reads\\tClean_Reads\\tUniquely_mapped\\tUniquely_mapped_ratio" > {output}
        for log in {RD}3.align/*_Log.final.out; do
            [ -e "$log" ] || continue
            sid=$(basename "$log" _Log.final.out)
            input_reads=$(fgrep "Number of input reads" "$log" | awk '{{print $NF}}')
            unique_num=$(fgrep "Uniquely mapped reads number" "$log" | awk '{{print $NF}}')
            unique_pct=$(fgrep "Uniquely mapped reads %" "$log" | awk '{{print $NF}}')
            total="NA"
            for tr in "{RD}2.cleandata/trim/${{sid}}_1_trimming_report.txt" \\
                      "{RD}2.cleandata/trim/${{sid}}_trimming_report.txt"; do
                if [ -f "$tr" ]; then
                    total=$(fgrep "Total reads processed:" "$tr" | head -1 | awk '{{print $NF}}' | sed 's/,//g')
                    break
                fi
            done
            printf '%s\\t%s\\t%s\\t%s\\t%s\\n' "$sid" "$total" "$input_reads" "$unique_num" "$unique_pct" >> {output}
        done
        """


rule check_strandedness:
    input:
        bam=R("3.align/{sample}_Aligned.sortedByCoord.out.bam"),
        bed=res("bed"),
    output:
        strand=R("3.align/{sample}.strandedness"),
        infer=R("3.align/{sample}_infer_experiment.out"),
    log:
        R("logs/check_strandedness/{sample}_check_strandedness.log.txt"),
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
