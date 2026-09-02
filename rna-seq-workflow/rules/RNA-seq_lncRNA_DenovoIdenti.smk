
rule runStringtie:
    input:
        bam = "3.align/{sample}_Aligned.sortedByCoord.out.bam",
        ref = {config["gtf"]},
        strandedness = "3.align/{sample}.strandedness",
    output:
        "4.LncRNA/4.1.Assembly_stringtie/{sample}.gtf",
    params:
        " --conservative -j 5 "
    conda:
        "lncrna"
    log:
        "logs/runStringtie/{sample}_stringtie.log.txt"
    shell:
        """
            strandedness=`cat 3.align/{wildcards.sample}.strandedness`
            ## run stringtie depends on  different library strandedness type 
            if [ ${{strandedness}} == "firststrand" ];then  ## 
                stringtie -p {config[threads]} {params} --rf -o {output} -G {input.ref} {input.bam}  >> {log} 2>&1
            elif [ ${{strandedness}} == "secondstrand" ];then ## 
                stringtie -p {config[threads]} {params} --fr -o {output} -G {input.ref} {input.bam}  >> {log} 2>&1
            elif [ ${{strandedness}} == "unstrand" ];then
                stringtie -p {config[threads]} {params} -o {output} -G {input.ref} {input.bam}  >> {log} 2>&1
            fi
        """

rule gtf_merge_compare:
    input:
        gtfIN = expand("4.LncRNA/4.1.Assembly_stringtie/{sample}.gtf",sample=SAMPLES),
        ref = {config["gtf_PcGs"]},
        genome = {config["genome"]},
    output:
        "4.LncRNA/4.1.Assembly_stringtie/merged.gtf",
        "4.LncRNA/4.1.Assembly_stringtie/merged.gtf.compare.classcode_u.fa"
    conda:
        "lncrna"
    log:
        "logs/gtf_merge_compare.log.txt"
    shell:
        """
            ## loading the defined functions
            source /home/chengyu/workflows/snakemake/rna-seq-workflow/scripts/lncRNA_functions.sh 
            ## merge all gtf files
            stringtie --merge -G {input.ref} -F 0.1 -T 0.1 -o 4.LncRNA/4.1.Assembly_stringtie/merged.gtf {input.gtfIN} >> {log} 2>&1
            ## compare with reference gene annotaion, get potential lincRNA transcripts
            runGFFcompare {input.ref} 4.LncRNA/4.1.Assembly_stringtie/merged.gtf {input.genome} >> {log} 2>&1
        """

rule Coding_predict:
    input:
        fa = "4.LncRNA/4.1.Assembly_stringtie/merged.gtf.compare.classcode_u.fa",
        gtfIN = "4.LncRNA/4.1.Assembly_stringtie/merged.gtf",
        genome = {config["genome"]},
    output:
        "4.LncRNA/4.2.Coding_predict/noncoding.transciptIDs.fa",
        "4.LncRNA/4.2.Coding_predict/noncoding.transciptIDs.gtf",
    log:
        "logs/Coding_predict.log.txt"

    conda:
        "lncrna"
    shell:
        """
            source /home/chengyu/workflows/snakemake/rna-seq-workflow/scripts/lncRNA_functions.sh
            ## CPC2 prediction 
            if [ ! -f "merged.gtf.compare.classcode_u.fa.CPC2.noncodingID" ];then
                runCPC2 {input.fa} >> {log} 2>&1
            fi
            ## CNCI prediction 
            if [ ! -f "merged.gtf.compare.classcode_u.fa.CNCI.noncodingID" ];then
                runCNCI {input.fa} {config[threads]} >> {log} 2>&1
            fi
            ## Tr length filter 
            cat {input.fa} | /home/chengyu/soft/bin/bioawk -c fastx 'length($seq)>=200{{print $name}}' > 4.LncRNA/4.2.Coding_predict/Tr_Len_Filtered.transciptIDs 
            
            ## combined filters
            cat 4.LncRNA/4.1.Assembly_stringtie/merged.gtf.compare.classcode_u.fa*noncodingID|fgrep -w -f 4.LncRNA/4.2.Coding_predict/Tr_Len_Filtered.transciptIDs|sort -u > 4.LncRNA/4.2.Coding_predict/noncoding.transciptIDs  

            ## get seqs
            cat {input.gtfIN}|fgrep -w -f 4.LncRNA/4.2.Coding_predict/noncoding.transciptIDs > 4.LncRNA/4.2.Coding_predict/noncoding.transciptIDs.gtf 
            getFasta 4.LncRNA/4.2.Coding_predict/noncoding.transciptIDs.gtf {input.genome} >> {log} 2>&1

        """

rule Pfam_search:
    input:
        fa = "4.LncRNA/4.2.Coding_predict/noncoding.transciptIDs.fa",
        gtfIN = "4.LncRNA/4.2.Coding_predict/noncoding.transciptIDs.gtf",
        genome = {config["genome"]},
    output:
        # protected("4.LncRNA/4.3.Pfam_search/pfam_scan.output"),
        # "4.LncRNA/4.3.Pfam_search/pfam_scan_filtered.fa",
        # "4.LncRNA/4.3.Pfam_search/pfam_scan_filtered.gtf",
        "4.LncRNA/4.3.Pfam_search/pfam_scan.hitIDs",
    log:
        "logs/Pfam_search.log.txt"
    conda:
        "lncrna"
    threads:
        30
    shell:
        """
            source /home/chengyu/workflows/snakemake/rna-seq-workflow/scripts/lncRNA_functions.sh
            # echo {threads} >> {log} 2>&1
            pfam_scan.pl -translate -fasta {input.fa} -dir {config[pfam_DB]} -outfile 4.LncRNA/4.3.Pfam_search/pfam_scan.output -as -cpu {threads} >> {log} 2>&1
            ## pfam_scan results parse
            cat 4.LncRNA/4.3.Pfam_search/pfam_scan.output |grep -v "^#"|grep -v '^\s*$'|awk '($13 < 1e-5){{print $1}}'|awk -F"." -v OFS="." '{{$4=null;print}}'|sed 's/\.$//g'|sort -u > 4.LncRNA/4.3.Pfam_search/pfam_scan.hitIDs 
            ## get seqs
            # cat {input.gtfIN} |fgrep -w -v -f 4.LncRNA/4.3.Pfam_search/pfam_scan.hitIDs > 4.LncRNA/4.3.Pfam_search/pfam_scan_filtered.gtf 
            # getFasta 4.LncRNA/4.3.Pfam_search/pfam_scan_filtered.gtf {input.genome} >> {log} 2>&1
        """

rule Nr_search:
    input:
        fa = "4.LncRNA/4.2.Coding_predict/noncoding.transciptIDs.fa",
        gtfIN = "4.LncRNA/4.2.Coding_predict/noncoding.transciptIDs.gtf",
        genome = {config["genome"]},
    output:
        # protected("4.LncRNA/4.4.Nr_search/diamond.output"),
        # "4.LncRNA/4.4.Nr_search/nr_filtered.gtf",
        # "4.LncRNA/4.4.Nr_search/nr_filtered.fa",
        "4.LncRNA/4.4.Nr_search/nr.hitIDs",
    conda:
        "lncrna"
    log:
        "logs/Nr_search/log.txt"
    threads:
        30
    shell:
        """

            source /home/chengyu/workflows/snakemake/rna-seq-workflow/scripts/lncRNA_functions.sh
            echo {threads} >> {log} 2>&1
            diamond blastx --threads {threads} --db {config[nr_diamond_DB]} -q {input.fa} --out 4.LncRNA/4.4.Nr_search/diamond.output >> {log} 2>&1
            cat 4.LncRNA/4.4.Nr_search/diamond.output| awk '{{print $1}}' |sort -u > 4.LncRNA/4.4.Nr_search/nr.hitIDs 
            # cat {input.gtfIN} |fgrep -w -v -f 4.LncRNA/4.4.Nr_search/nr.hitIDs >  4.LncRNA/4.4.Nr_search/nr_filtered.gtf 
            # getFasta 4.LncRNA/4.4.Nr_search/nr_filtered.gtf {input.genome}  >> {log} 2>&1

        """

rule Final_lncRNA:
    input:
        "4.LncRNA/4.3.Pfam_search/pfam_scan.hitIDs",
        "4.LncRNA/4.4.Nr_search/nr.hitIDs",
        gtf_noncode = "4.LncRNA/4.2.Coding_predict/noncoding.transciptIDs.gtf",
        genome = {config["genome"]},
    output:
        "4.LncRNA/4.5.Final_lncRNA/final_lncRNA.fa",
        "4.LncRNA/4.5.Final_lncRNA/final_lncRNA.gtf",
    conda:
        "lncrna"
    log:
        "logs/Final_lnc/log.txt"
    shell:
        """
            source /home/chengyu/workflows/snakemake/rna-seq-workflow/scripts/lncRNA_functions.sh
            cat 4.LncRNA/4.3.Pfam_search/pfam_scan.hitIDs 4.LncRNA/4.4.Nr_search/nr.hitIDs |sort -u > 4.LncRNA/4.5.Final_lncRNA/all.hitsIDs
            cat {input.gtf_noncode} |fgrep -w -v -f 4.LncRNA/4.5.Final_lncRNA/all.hitsIDs > 4.LncRNA/4.5.Final_lncRNA/final_lncRNA.gtf
            getFasta 4.LncRNA/4.5.Final_lncRNA/final_lncRNA.gtf {input.genome} >> {log} 2>&1
            cat {config[gtf]} 4.LncRNA/4.5.Final_lncRNA/final_lncRNA.gtf|sortBed -i - > 4.LncRNA/4.5.Final_lncRNA/ref_withLncRNA.gtf
        """

rule expression:
    input:
        bam = "3.align/{sample}_Aligned.sortedByCoord.out.bam",
        gtf = "4.LncRNA/4.5.Final_lncRNA/ref_withLncRNA.gtf",
        strand = "3.align/{sample}.strandedness",
        # layout = (isPAIRED),
    output:
        "5.expression/lncRNA/{sample}.count",
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
            /opt/anaconda3/envs/R/bin/Rscript /share/workflows/rna-seq/scripts/run-featurecounts.R -t {config[threads]} -b {input.bam} -g 4.LncRNA/4.5.Final_lncRNA/ref_withLncRNA.gtf -s ${{strand}} -i ${{isPairedEnd}} -o 5.expression/lncRNA/{wildcards.sample}      
            touch 5.expression/lncRNA/{wildcards.sample}.flag     
        """

rule count_merge2:
    input:
        expand("5.expression/lncRNA/{sample}.count",sample = SAMPLES)
    output:
        "5.expression/lncRNA/count.matrix.tsv",
        "5.expression/lncRNA/GeneExpression_TPM.xls",
    log:
        "logs/count_merge/log.txt"

    shell:
        """
			# bash /home/chengyu/workflows/snakemake/rna-seq-workflow/scripts/featureCount.R_result_merge.sh 5.expression/lncRNA
			python /home/chengyu/myscripts/featureCount.R_result_merge.py 5.expression/lncRNA
        """
