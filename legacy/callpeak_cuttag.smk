## skip dedup_picard step 

rule callpeak_macs2:
    input:
        expand("3.align/bowtie2/{sample}_sorted.bam",sample=SAMPLES),
        grouplist = {config["grouplist"]},
    output:
        "4.peak/all_peakfiles.txt"
    log:
        "logs/callpeak_macs2.txt"
    params:
        macs2 = " --keep-dup -f BAMPE --nomodel -B -g 3.7e8 -q 0.5 "
    shell:
        """
            cat {input.grouplist}|sed '1d'|while read line;
                do
                    ID_control=`echo "${{line}}"|cut -d, -f2`
                    ID_treat=`echo "${{line}}"|cut -d, -f1`
                    GROUP_name="${{ID_treat}}_${{ID_control}}"
                # macs2 call peak
                    for j in $ID_control $ID_treat
                        do
                            macs2 callpeak {params.macs2} -c 3.align/bowtie2/${{ID_control}}_rmdup.bam -t 3.align/bowtie2/${{ID_treat}}_rmdup.bam --outdir 4.peak/ -n ${{GROUP_name}}  >> {log} 2>&1
                        done
                done
                ls 4.peak/*.narrowPeak > 4.peak/all_peakfiles.txt
        """

rule peakAnno:
    input:
        peakfiles = "4.peak/all_peakfiles.txt",
        gtf = {config["gtf"]},
    output:
        "4.peak/peakAnno.flag"
    log:
        "logs/peakAnno.log.txt"
    shell:
        """
            Rscript /home/chengyu/workflows/annoPeak_chipseeker.R {input.gtf} {input.peakfiles} >> {log} 2>&1
            touch 4.peak/peakAnno.flag
        """
