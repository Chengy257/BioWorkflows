###############################################
## QC and alignment rules (shared by all pipelines):
##   trim (unified naming convention + FastQC, fix P1-1) -> MultiQC whole-pipeline summary (P2-9)
##   -> STAR (direct SortedByCoordinate output, P1-4) -> strandedness inference
## Depends on SAMPLES / R() / _raw_reads / STAR_ARGS etc. defined in common.smk
###############################################


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
    params:
        outdir=lambda wc, output: os.path.dirname(output.fq_out),
        fastqc_dir=lambda wc, output: os.path.dirname(output.fastqc_html),
        trim_galore=tool("trim_galore", "trim_galore"),
    threads:
        rthreads("trim")
    resources:
        mem_mb=rmem("trim"),
        runtime_min=rruntime("trim"),
        runtime_sec=rruntime_sec("trim"),
    shell:
        """
        mkdir -p {params.fastqc_dir}
        {params.trim_galore} -q 30 --stringency 3 -e 0.1 --gzip -j {threads} -o {params.outdir}/ {input.fq} \\
            --fastqc --fastqc_args "--outdir {params.fastqc_dir} " >> {log} 2>&1
        ## Trim reports are generated with the original file names; rename them uniformly to canonical names keyed by sample id
        mv "{params.outdir}/$(basename "{input.fq}")_trimming_report.txt" {output.report}
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
    params:
        outdir=lambda wc, output: os.path.dirname(output.fq1_out),
        fastqc_dir=lambda wc, output: os.path.dirname(output.fastqc1_html),
        trim_galore=tool("trim_galore", "trim_galore"),
    threads:
        rthreads("trim")
    resources:
        mem_mb=rmem("trim"),
        runtime_min=rruntime("trim"),
        runtime_sec=rruntime_sec("trim"),
    shell:
        """
        mkdir -p {params.fastqc_dir}
        {params.trim_galore} -q 30 --stringency 3 -e 0.1 --gzip -j {threads} -o {params.outdir}/ --paired {input.fq1} {input.fq2} \\
            --fastqc --fastqc_args "--outdir {params.fastqc_dir} " >> {log} 2>&1
        mv "{params.outdir}/$(basename "{input.fq1}")_trimming_report.txt" {output.report1}
        mv "{params.outdir}/$(basename "{input.fq2}")_trimming_report.txt" {output.report2}
        """



rule multiqc:
    input:
        multiqc_inputs,
    output:
        R("multiqc/multiqc_report.html"),
    log:
        R("logs/multiqc/multiqc_log.txt"),
    params:
        outdir=lambda wc, output: os.path.dirname(output[0]),
        scan_dir=lambda wc, output: os.path.dirname(os.path.dirname(output[0])),
        config_file=os.path.join(WORKFLOW_DIR, "multiqc_config.yaml"),
        multiqc=tool("multiqc", "multiqc"),
    threads:
        rthreads("multiqc")
    resources:
        mem_mb=rmem("multiqc"),
        runtime_min=rruntime("multiqc"),
        runtime_sec=rruntime_sec("multiqc"),
    shell:
        """
        mkdir -p {params.outdir}
        {params.multiqc} --force --config {params.config_file} -o {params.outdir} {params.scan_dir} >> {log} 2>&1
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
    params:
        index_dir=lambda wc, output: os.path.dirname(output.sa),
        star=tool("star", "STAR"),
    threads:
        rthreads("star_index")
    resources:
        mem_mb=rmem("star_index"),
        runtime_min=rruntime("star_index"),
        runtime_sec=rruntime_sec("star_index"),
    shell:
        """
        mkdir -p {params.index_dir}
        {params.star} --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.index_dir}/ \\
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
    threads:
        rthreads("star_align")
    resources:
        mem_mb=rmem("star_align"),
        runtime_min=rruntime("star_align"),
        runtime_sec=rruntime_sec("star_align"),
    params:
        star_args=STAR_ARGS,
        reads=lambda wc: trimmed_reads(wc.sample),
        prefix=lambda wc: R(f"3.align/{wc.sample}_"),
        align_dir=lambda wc, output: os.path.dirname(output.bam),
        index_dir=lambda wc, input: os.path.dirname(input.index),
        star=tool("star", "STAR"),
        samtools=tool("samtools", "samtools"),
    shell:
        """
        mkdir -p {params.align_dir}
        {params.star} {params.star_args} \\
             --runThreadN {threads} --genomeDir {params.index_dir}/ \\
             --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMattributes All \\
             --readFilesIn {params.reads} \\
             --outFileNamePrefix {params.prefix} >> {log} 2>&1
        {params.samtools} index -@ {threads} {output.bam} >> {log} 2>&1
        """


rule Mapping_stat:
    input:
        bams=expand(R("3.align/{sample}_Aligned.sortedByCoord.out.bam"), sample=SAMPLES),
        reports=all_trim_reports,
    output:
        R("3.align/mapping_stat.xls"),
    log:
        R("logs/Mapping_stat.log.txt")
    params:
        align_dir=lambda wc, output: os.path.dirname(output[0]),
        trim_dir=R("2.cleandata/trim"),
    threads:
        rthreads("mapping_stat")
    resources:
        mem_mb=rmem("mapping_stat"),
        runtime_min=rruntime("mapping_stat"),
        runtime_sec=rruntime_sec("mapping_stat"),
    shell:
        """
        ## Parse per sample: sample id comes from the log file name (underscores kept); the trim total is the R1 count
        ## from the canonical report (read pairs for PE data)
        echo -e "ID\\tTotal_Reads\\tClean_Reads\\tUniquely_mapped\\tUniquely_mapped_ratio" > {output}
        for log in {params.align_dir}/*_Log.final.out; do
            [ -e "$log" ] || continue
            sid=$(basename "$log" _Log.final.out)
            input_reads=$(fgrep "Number of input reads" "$log" | awk '{{print $NF}}')
            unique_num=$(fgrep "Uniquely mapped reads number" "$log" | awk '{{print $NF}}')
            unique_pct=$(fgrep "Uniquely mapped reads %" "$log" | awk '{{print $NF}}')
            total="NA"
            for tr in "{params.trim_dir}/${{sid}}_1_trimming_report.txt" \\
                      "{params.trim_dir}/${{sid}}_trimming_report.txt"; do
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
    params:
        infer_experiment=tool("infer_experiment", "infer_experiment.py"),
    threads:
        rthreads("strandedness")
    resources:
        mem_mb=rmem("strandedness"),
        runtime_min=rruntime("strandedness"),
        runtime_sec=rruntime_sec("strandedness"),
    shell:
        """
        {params.infer_experiment} -r {input.bed} -i {input.bam} > {output.infer} 2> {log}
        tail -2 {output.infer} | awk '{{print $NF}}' | \\
            awk '{{f1=$0;getline;f2=$0; \\
                if(f1-f2 > 0.4) print "secondstrand"; \\
                else if(f2-f1 > 0.4) print "firststrand"; \\
                else print "unstrand"}}' > {output.strand}
        """
