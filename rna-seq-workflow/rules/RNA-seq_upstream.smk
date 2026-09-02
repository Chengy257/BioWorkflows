

rule trimAdapter_SE: 
    input: 
        fq="1.rawdata/{sample}.fastq.gz",

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
            trim_galore -q 30 --stringency 3 -e 0.1 --gzip -j {config[threads]} -o 2.cleandata/trim/ {input.fq}  >> {log} 2>&1
            touch 2.cleandata/trim/fastqc/{wildcards.sample}.fastqc.flag
        """

rule trimAdapter_PE: 
    input: 
        fq1="1.rawdata/{sample}_1.fastq.gz",
        fq2="1.rawdata/{sample}_2.fastq.gz",
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
            trim_galore -q 30 --stringency 3 -e 0.1 --gzip -j {config[threads]} -o 2.cleandata/trim/ --paired {input.fq1} {input.fq2} >> {log} 2>&1 
            touch 2.cleandata/trim/fastqc/{wildcards.sample}.fastqc.flag
        """

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
    conda:
        "rna-seq"
    shell:
        """
            STAR --runThreadN {config[threads]} --twopassMode Basic \
                 --genomeLoad NoSharedMemory --genomeDir 0.index/star_genome/ \
                 --readFilesCommand zcat --outSAMtype BAM Unsorted --outSAMattributes All \
                 --quantMode GeneCounts --readFilesIn  2.cleandata/trim/{wildcards.sample}*fq.gz \
                 --outFileNamePrefix 3.align/{wildcards.sample}_ --outSAMattrIHstart 0 \
                 --outReadsUnmapped Fastx >> {log} 2>&1
            samtools sort -@ {config[threads]} -O BAM -o {output.bam} 3.align/{wildcards.sample}_Aligned.out.bam && rm 3.align/{wildcards.sample}_Aligned.out.bam >> {log} 2>&1
            samtools index -@ {config[threads]} 3.align/{wildcards.sample}_Aligned.sortedByCoord.out.bam >> {log} 2>&1
        """

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
            cat 2.cleandata/trim/*_trimming_report.txt|fgrep "Total reads processed:"|awk '{{print $NF}}'|sed 's/,//g' > temp.mapping_stat.total
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

rule featureCount_R:
    input:
        bam = "3.align/{sample}_Aligned.sortedByCoord.out.bam",
        gtf = {config["gtf"]},
        strand = "3.align/{sample}.strandedness",
        # layout = (isPAIRED),
    output:
        "4.expression/{sample}.count",
    conda:
        "rna-seq",
    log:
        "logs/featureCount_R/{sample}.log.txt"
    shell:
        """
            ## 
            strandedness=`cat {input.strand}|head -1|awk '{{print $1}}'`            
            ## run featureCount_R depends on  different library strandedness type 
            if [ ${{strandedness}} == "firststrand" ];then  ## 
                strand="2"
            elif [ ${{strandedness}} == "secondstrand" ];then ## 
                strand="1"
            elif [ ${{strandedness}} == "unstrand" ];then
                strand="0"
            fi
            ## 
            ends=`ls 1.rawdata/{wildcards.sample}*.gz|wc -l`
            if [ "${{ends}}" == 2 ];then
                isPairedEnd="True"
            elif [ "${{ends}}" == 1 ];then
                isPairedEnd="False"
            fi   
            ##
            Rscript /share/workflows/rna-seq/scripts/run-featurecounts.R -t {config[threads]} -b {input.bam} -g {input.gtf} -s ${{strand}} -i ${{isPairedEnd}} -o 4.expression/{wildcards.sample}      
            touch 4.expression/{wildcards.sample}.flag     
        """

rule count_merge:
    input:
        expand("4.expression/{sample}.count",sample = SAMPLES)
    output:
        "4.expression/count.matrix.tsv",
        "4.expression/GeneExpression_TPM.xls",
    log:
        "logs/count_merge/log.txt"
    shell:
        """
            bash /home/chengyu/workflows/snakemake/rna-seq-workflow/scripts/featureCount.R_result_merge.sh 4.expression
        """
