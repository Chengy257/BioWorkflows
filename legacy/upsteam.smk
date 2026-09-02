rule trimAdapter: 
    input: 
        fq1 = "1.rawdata/{sample}_1.fq.gz",
        fq2 = "1.rawdata/{sample}_2.fq.gz",   
    output: 
        "2.cleandata/{sample}_1_val_1.fq.gz",
        "2.cleandata/{sample}_2_val_2.fq.gz",
    # threads:
    #     int({config[threads]})
    log:
        "logs/trim/{sample}_log.txt",
    conda:
        "chip"
    shell: 
        """
            trim_galore -q 25 --stringency 3 -e 0.1 --gzip \
                -j {config[threads]} -o 2.cleandata/ \
                --paired {input.fq1} {input.fq2} >> {log} 2>&1  
        """

rule fastQC:    
    input: 
        expand("2.cleandata/{sample}_{R}_val_{R}.fq.gz",sample=SAMPLES,R=["1","2"]),  
    output:
        expand("2.cleandata/fastqc/{sample}_{R}_val_{R}_fastqc.html",sample=SAMPLES,R=["1","2"]),
        expand("2.cleandata/fastqc/{sample}_{R}_val_{R}_fastqc.zip",sample=SAMPLES,R=["1","2"]),        
    log: 
        "logs/fastqc/fastqc_log.txt"    
    conda:
        "chip"
    shell: 
        """
            fastqc -f fastq -t {config[threads]} -o 2.cleandata/fastqc/ {input} >> {log} 2>&1
        """       

rule multiQC:
    input: 
        expand("2.cleandata/fastqc/{sample}_{R}_val_{R}_fastqc.zip", sample=SAMPLES,R=["1","2"]),
        expand("2.cleandata/fastqc/{sample}_{R}_val_{R}_fastqc.html", sample=SAMPLES,R=["1","2"]),        
    output:
        "2.cleandata/fastqc/multiqc/multiqc_report.html",
    log:
        "logs/multiQC.log.txt"
    conda:
        "chip"
    shell:
        """
            multiqc --force 2.cleandata/fastqc/ -o 2.cleandata/fastqc/multiqc
        """

rule bowtie2_index:
    input:
        GENOME = {config["GENOME"]}
    output:
        "0.index/bowtie2.1.bt2",
    conda:
        "chip"
    log:
        "logs/bowtie2_index.log.txt"
    shell:
        """
            bowtie2-build --threads {config[threads]} {input.GENOME} 0.index/bowtie2 >> {log} 2>&1 
        """

rule bowtie2_mapping:
    input:
        index="0.index/bowtie2.1.bt2",
        fq1="2.cleandata/{sample}_1_val_1.fq.gz",
        fq2="2.cleandata/{sample}_2_val_2.fq.gz",        
    output:
        "3.align/bowtie2/{sample}_sorted.bam",
    params:
        " --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 "
    conda:
        "chip"
    log:
        "logs/bowtie2_mapping/{sample}.logs.txt"
    shell:
        """
            bowtie2 -x 0.index/bowtie2 -p {config[threads]} {params} -1 {input.fq1} -2 {input.fq2} --un-conc-gz 3.align/bowtie2/{wildcards.sample}_unmapped.fq.gz -S |/opt/anaconda3/bin/samtools view -@ {config[threads]} -q 5 -bS -o 3.align/bowtie2/{wildcards.sample}.bam - >> {log} 2>&1  ## Simply filter the bam with MAPQ (mapping quality of the reads), 5 or 10 is usually reasonable
            /opt/anaconda3/bin/samtools sort -@ {config[threads]} -O BAM -o 3.align/bowtie2/{wildcards.sample}_sorted.bam 3.align/bowtie2/{wildcards.sample}.bam  && rm 3.align/bowtie2/{wildcards.sample}.bam >> {log} 2>&1 
            /opt/anaconda3/bin/samtools index -@ {config[threads]} 3.align/bowtie2/{wildcards.sample}_sorted.bam >> {log} 2>&1
        """
