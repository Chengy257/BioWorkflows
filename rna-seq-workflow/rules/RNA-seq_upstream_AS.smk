


rule trimAdapter_SE: 
    input: 
        fq="1.rawdata/{sample}.fq.gz",

    output: 
        "2.cleandata/trim/{sample}_trimmed.fq.gz",
        "2.cleandata/trim/fastqc/{sample}.fastqc.flag",
        # "flag/trim/{sample}_SE_trimmed.flag",
    log:
        "logs/trim/{sample}_log.txt",
    conda:
        "rna-seq",
    shell: 
        """
            mkdir -p 2.cleandata/trim/fastqc
            trim_galore -q 30 --stringency 3 -e 0.1 --gzip -j {config[threads]} -o 2.cleandata/trim/ {input.fq}  --fastqc --fastqc_args "--outdir 2.cleandata/trim/fastqc "  >> {log} 2>&1
            touch 2.cleandata/trim/fastqc/{wildcards.sample}.fastqc.flag
        """

rule trimAdapter_PE: 
    input: 
        fq1="1.rawdata/{sample}_1.fq.gz",
        fq2="1.rawdata/{sample}_2.fq.gz",
        # csv=config["SampleListFile"]
    output: 
        "2.cleandata/trim/{sample}_1_val_1.fq.gz",
        "2.cleandata/trim/{sample}_2_val_2.fq.gz",
        "2.cleandata/trim/fastqc/{sample}.fastqc.flag",
    log:
        "logs/trim/{sample}_log.txt",
    conda:
        "rna-seq",
    shell: 
        """
            mkdir -p 2.cleandata/trim/fastqc
            trim_galore -q 30 --stringency 3 -e 0.1 --gzip -j {config[threads]} -o 2.cleandata/trim/ --paired {input.fq1} {input.fq2}  --fastqc --fastqc_args "--outdir 2.cleandata/trim/fastqc " >> {log} 2>&1 
            touch 2.cleandata/trim/fastqc/{wildcards.sample}.fastqc.flag
        """
# ruleorder: trimAdapter_PE > trimAdapter_SE

rule multiQC_cleaned:
    input: 
        expand("2.cleandata/trim/fastqc/{sample}.fastqc.flag",sample=SAMPLES)  
    log: 
        "logs/fastqc/multiqc_log.txt" 
    output:
        "2.cleandata/trim/fastqc/multiqc_report.html",
    conda:
        "common"
    shell:
        "multiqc --force 2.cleandata/trim/fastqc/ -o 2.cleandata/trim/fastqc/ >> {log} 2>&1 "

rule STAR_index:
    input: 
        GENOME={config["genome"]},
        GTF={config["gtf"]},
    output: 
        "0.index/star_genome/SAindex",
        star_index=directory("0.index/star_genome/"),
    log:
        "logs/index/star_index_log.txt",
    conda:
        "rna-seq",
    shell:
        """
            STAR --runThreadN {config[threads]} --runMode genomeGenerate --genomeDir {output.star_index} \
                 --genomeFastaFiles {input.GENOME} --sjdbGTFfile {input.GTF} ; 2>{log}
        """


import subprocess
def get_fastq(wildcards):
    command = "ls 1.rawdata/" + wildcards.sample + "*.gz |wc -l" 
    result = subprocess.Popen(command, shell=True,stdout=subprocess.PIPE,stderr = subprocess.STDOUT,encoding='utf-8',text=True)
    num = result.communicate()[0].strip() 
    if num == "2":
        return expand("2.cleandata/trim/{sample}_{R}_val_{R}.fq.gz",sample=wildcards.sample,R=["1","2"])
    elif num == "1":
        return expand("2.cleandata/trim/{sample}_trimmed.fq.gz",sample=wildcards.sample)  
    else:
        pass

print(unpack(get_fastq))

rule runSTAR:
    input:
        unpack(get_fastq),
        index = "0.index/star_genome/SAindex",
        # fq1 = "2.cleandata/{sample}_1_val_1.fq.gz",
        # fq2 = "2.cleandata/{sample}_2_val_2.fq.gz",        
    output:
        bam = "3.align/{sample}_Aligned.sortedByCoord.out.bam",
        bai = "3.align/{sample}_Aligned.sortedByCoord.out.bam.bai",
    log: 
        "logs/runSTAR/{sample}.log.txt",
    params:
        star = "--twopassMode Basic --outFilterType BySJout --alignIntronMin 20 --alignIntronMax 5000 --alignMatesGapMax 5000 --outFilterMatchNminOverLread 0.66 --outFilterScoreMinOverLread 0.66 --winAnchorMultimapNmax 70 --seedSearchStartLmax 45 --outSAMattrIHstart 0 --outSAMstrandField intronMotif --genomeLoad NoSharedMemory --quantMode TranscriptomeSAM GeneCounts "
    conda:
        "rna-seq"
    shell:
        """
            STAR {params.star} \
                 --runThreadN {config[threads]} --genomeDir 0.index/star_genome/ \
                 --readFilesCommand zcat --outSAMtype BAM Unsorted --outSAMattributes All \
                 --readFilesIn  2.cleandata/trim/{wildcards.sample}*fq.gz \
                 --outFileNamePrefix 3.align/{wildcards.sample}_ \
                 --outReadsUnmapped Fastx  >> {log} 2>&1
            samtools sort -@ {config[threads]} -O BAM -o {output.bam} 3.align/{wildcards.sample}_Aligned.out.bam && rm 3.align/{wildcards.sample}_Aligned.out.bam >> {log} 2>&1
            samtools index -@ {config[threads]} 3.align/{wildcards.sample}_Aligned.sortedByCoord.out.bam >> {log} 2>&1
        """

#                 --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
#                 --alignIntronMin 20 --alignIntronMax ${max_intron_size} \
#                 --alignMatesGapMax ${max_intron_size} \
#                 --outFilterMatchNminOverLread 0.66 --outFilterScoreMinOverLread 0.66 \
#                 --winAnchorMultimapNmax 70 --seedSearchStartLmax 45 \
#                 --outSAMattrIHstart 0 --outSAMstrandField intronMotif \
#                 --genomeLoad LoadAndKeep --outReadsUnmapped Fastx \
#                 --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM GeneCounts"

rule Mapping_stat:
    input:
        expand("3.align/{sample}_Aligned.sortedByCoord.out.bam",sample=SAMPLES),
    output:
        "3.align/mapping_stat.xls",
    log:
        "logs/Mapping_stat.log.txt"
    shell:
        """
            echo -e "ID\tTotal_Reads\tClean_Reads\tUniquely_mapped\tUniquely_mapped_ratio" > 3.align/mapping_stat.xls
            ls 3.align/*_Log.final.out|xargs -i basename {{}}|cut -d_ -f1 > temp.mapping_stat.id
            cat 2.cleandata/trim/*_trimming_report.txt|fgrep "Total reads processed:"|awk '{{print $NF}}'|sed 's/,//g'|head -1 > temp.mapping_stat.total
            cat 3.align/*_Log.final.out |fgrep "Number of input reads"|awk '{{print $NF}}' > temp.mapping_stat.input
            cat 3.align/*_Log.final.out |fgrep "Uniquely mapped reads number"|awk '{{print $NF}}' > temp.mapping_stat.unique
            cat 3.align/*_Log.final.out |fgrep "Uniquely mapped reads %"|awk '{{print $NF}}' > temp.mapping_stat.unique_ratio
            paste -d "\t" temp.mapping_stat.id temp.mapping_stat.total temp.mapping_stat.input temp.mapping_stat.unique temp.mapping_stat.unique_ratio >> 3.align/mapping_stat.xls && rm temp.mapping_stat*
        """

rule check_strandedness:
    input:
        bam="3.align/{sample}_Aligned.sortedByCoord.out.bam",
        bed={config["bed"]},
    output:
        "3.align/{sample}.strandedness",
    conda:
        "rna-seq",
    # log:
    #     "logs/check_strandedness/{sample}_check_strandedness.log.txt",
    shell:
        """
            infer_experiment.py -r {input.bed} -i {input.bam} > 3.align/{wildcards.sample}_infer_experiment.out  
            cat 3.align/{wildcards.sample}_infer_experiment.out|tail -2|awk '{{print $NF}}'| \
                awk '{{f1=$0;getline;f2=$0; \
                    if(f1-f2 > 0.4) {{print "secondstrand"}}  \
                    else if(f2-f1 > 0.4) {{print "firststrand"}} \
                    else {{print "unstrand"}} }}' > 3.align/{wildcards.sample}.strandedness 
            # strandedness=`cat 3.align/{wildcards.sample}.strandedness`
        """

rule Assemble:
    input:
        bam = "3.align/{sample}_Aligned.sortedByCoord.out.bam",
        ref = {config["gtf"]},
        strandedness = "3.align/{sample}.strandedness",
    output:
        "4.assembly/4.1.Assembly_stringtie/{sample}.gtf",
    conda:
        "lncrna"
    log:
        "logs/Assemble/{sample}_stringtie.log.txt"
    shell:
        """
            strandedness=`cat 3.align/{wildcards.sample}.strandedness`
            ## run stringtie depends on  different library strandedness type 
            if [ ${{strandedness}} == "firststrand" ];then  ## 
                stringtie -p {config[threads]} --rf -o {output} -G {input.ref} {input.bam}  >> {log} 2>&1
            elif [ ${{strandedness}} == "secondstrand" ];then ## 
                stringtie -p {config[threads]} --fr -o {output} -G {input.ref} {input.bam}  >> {log} 2>&1
            elif [ ${{strandedness}} == "unstrand" ];then
                stringtie -p {config[threads]} -o {output} -G {input.ref} {input.bam}  >> {log} 2>&1
            fi
        """

rule gtf_merge:
    input:
        gtf = expand("4.assembly/4.1.Assembly_stringtie/{sample}.gtf",sample=SAMPLES),
        ref = {config["gtf"]},
        genome = {config["genome"]},
    output:
        "4.assembly/4.1.Assembly_stringtie/merged.gtf",
    conda:
        "lncrna"
    log:
        "logs/gtf_merge.log.txt"
    shell:
        """
            ## loading the defined functions
            source /home/chengyu/workflows/snakemake/rna-seq-workflow/scripts/lncRNA_functions.sh             
            ## merge all gtf files
            stringtie --merge -G {input.ref} -i -o 4.assembly/4.1.Assembly_stringtie/merged.gtf {input.gtf} >> {log} 2>&1
            ## compare with reference gene annotaion
            runGFFcompare {input.ref} 4.assembly/4.1.Assembly_stringtie/merged.gtf {input.genome} >> {log} 2>&1
        """

rule isoform_expr:
    input:
        bam = "3.align/{sample}_Aligned.sortedByCoord.out.bam",
        ref = "4.assembly/4.1.Assembly_stringtie/merged.gtf",
        strandedness = "3.align/{sample}.strandedness",
    output:
        gtf="4.assembly/4.2.IsoformExpr/{sample}.gtf",
        tab="4.assembly/4.2.IsoformExpr/{sample}.tab",
    conda:
        "lncrna"
    log:
        "logs/isoform_expr/{sample}_stringtie.log.txt"
    shell:
        """
            strandedness=`cat 3.align/{wildcards.sample}.strandedness`
            ## run stringtie depends on  different library strandedness type 
            if [ ${{strandedness}} == "firststrand" ];then  ## 
                stringtie -p {config[threads]} --rf -o {output.gtf} -e -G {input.ref} {input.bam} >> {log} 2>&1
            elif [ ${{strandedness}} == "secondstrand" ];then ## 
                stringtie -p {config[threads]} --fr -o {output.gtf} -e -G {input.ref} {input.bam} >> {log} 2>&1
            elif [ ${{strandedness}} == "unstrand" ];then
                stringtie -p {config[threads]} -o {output.gtf} -e -G {input.ref} {input.bam} >> {log} 2>&1
            fi
            ## extract expression
            cat {output.gtf}|grep -v "^#"|awk -v OFS="\t" 'BEGIN{{print "transcript","FPKM","TPM"}} {{if($3=="transcript"){{print $12,$(NF-2),$NF}}}}'|sed 's/[;"]//g' > {output.tab}

        """

# rule run_rMATs:

