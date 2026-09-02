
# rule featureCount_R:
#     input:
#         bam = "3.align/{sample}_Aligned.sortedByCoord.out.bam",
#         gtf = {config["gtf"]},
#         strand = "3.align/{sample}.strandedness",
#         # layout = (isPAIRED),
#     output:
#         "4.expression/{sample}.count",
#     conda:
#         "rna-seq",
#     log:
#         "logs/featureCount_R/{sample}.log.txt"
#     shell:
#         """
#             ## 
#             strandedness=`cat {input.strand}|head -1|awk '{{print $1}}'`            
#             ## run featureCount_R depends on  different library strandedness type 
#             if [ ${{strandedness}} == "firststrand" ];then  ## 
#                 strand="2"
#             elif [ ${{strandedness}} == "secondstrand" ];then ## 
#                 strand="1"
#             elif [ ${{strandedness}} == "unstrand" ];then
#                 strand="0"
#             fi
#             ## 
#             ends=`ls 1.rawdata/{wildcards.sample}*.gz|wc -l`
#             if [ "${{ends}}" == 2 ];then
#                 isPairedEnd="True"
#             elif [ "${{ends}}" == 1 ];then
#                 isPairedEnd="False"
#             fi   
#             ##
#             Rscript /share/workflows/rna-seq/scripts/run-featurecounts.R -t {config[threads]} -b {input.bam} -g {input.gtf} -s ${{strand}} -i ${{isPairedEnd}} -o 4.expression/{wildcards.sample}      
#             touch 4.expression/{wildcards.sample}.flag     
#         """

# rule count_merge:
#     input:
#         expand("4.expression/{sample}.count",sample = SAMPLES)
#     output:
#         "4.expression/count.matrix.tsv",
#         "4.expression/GeneExpression_TPM.xls",
#     log:
#         "logs/count_merge/log.txt"

#     shell:
#         """
#             bash /home/chengyu/workflows/snakemake/rna-seq-workflow/scripts/featureCount.R_result_merge.sh 4.expression
#         """

rule runDESeq2:
    input: 
        "4.expression/count.matrix.tsv",
        "sample_info.csv",
    output:
        # "5.DEG/Rplots.pdf",
        "5.DEG/flag.log",
    log: 
        "logs/DESeq2/runDESeq2.log.txt",
    shell:
        """
            /usr/bin/Rscript /share/workflows/rna-seq/scripts/runDESeq2_multiGroup_230705.R -c 4.expression/count.matrix.tsv -s sample_info.csv -f {config[FoldChange]} -b F -o 5.DEG -t {config[threads]} >> {log} 2>&1
            echo `date` " : All done!" > 5.DEG/flag.log

        """

rule GO_KEGG_enrich:
    input:
        "5.DEG/flag.log",
    output:
        "5.DEG/GO_KEGG_enrich/flag.log",
    log: 
        "logs/GO_KEGG_enrich/GO_KEGG_enrich.log.txt",
    shell:
        """
            bash /home/chengyu/workflows/snakemake/rna-seq-workflow/scripts/enrich.sh 5.DEG >> {log} 2>&1
        """

rule DEGgroupCompare:
    input:
        "5.DEG/GO_KEGG_enrich/flag.log",
    output:
        "6.DEGcompare/flag.log"
    log:
        "logs/DEGgroupCompare/DEGgroupCompare.log.txt"
    shell:
        """
            bash /home/chengyu/workflows/snakemake/rna-seq-workflow/scripts/DEGgroupCompare.sh >> {log} 2>&1
        """
