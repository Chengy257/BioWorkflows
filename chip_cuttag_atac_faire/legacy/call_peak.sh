#!/bin/bash
#########################################################################
# File Name: call_peak.sh
# Author: ChengYu
# Description: 
# Created Time: Sat 02 Mar 2024 09:12:31 PM CST
# Sun Mar  3 21:10:49 CST 2024: v1.0 by chengyu
#########################################################################
######################## DEFINE FUNCTIONS ########################
## Extract the fragment length estimate from column 3 of the cross-correlation scores file
frglen_spp(){
    ## Usage: frglen_spp [input treat bam]
    echo "["`date`"]: " "FUNCTION[frglen_spp]: processing " ${1}
    /opt/anaconda3/envs/phantompeakqualtools/bin/Rscript \
        /opt/anaconda3/envs/phantompeakqualtools/bin/run_spp.R \
        -p=${threads} -savp -c=${1} \
        -out=${1}_spp_cross_correlation.txt
    cat ${1}_spp_cross_correlation.txt| sed -r 's/,[^\t]+//g' |cut -f3|\
        cut -d, -f1 > ${1}_fragment_len
    cat ${1}_spp_cross_correlation.txt| cut -f9 >  ${1}_NSC
    cat ${1}_spp_cross_correlation.txt| cut -f10 >  ${1}_RSC
    echo "["`date`"]: " "Finished! "
}
    # NSC values range from a minimum of 1 to larger positive numbers. 1.1 is the critical threshold.
    # Datasets with NSC values much less than 1.1 (< 1.05) tend to have low signal to noise or few peaks 
    # (this could be biological eg.a factor that truly binds only a few sites in a particular tissue type 
    # OR it could be due to poor quality)
    # RSC values range from 0 to larger positive values. 1 is the critical threshold. RSC values significantly 
    # lower than 1 (< 0.8) tend to have low signal to noise. The low scores can be due to failed and poor 
    # quality ChIP, low read sequence quality and hence lots of mismappings, shallow sequencing depth 
    # (significantly below saturation) or a combination of these. Like the NSC, datasets with few binding 
    # sites (< 200) which is biologically justifiable also show low RSC scores.
    
MACS2_narrow_peak(){
    ## Usage: MACS2_narrow_peak [control bam file] [treatment bam file]
    echo "["`date`"]: " "FUNCTION[MACS2_narrow_peak]: processing " ${1}
    macs2 callpeak --keep-dup all -f BAMPE  \
          -B --SPMR -g 3.7e8 \
          --nomodel --extsize ${fragment_len} -q 0.5 \
          --call-summits --cutoff-analysis   \
          -t ${1} -c ${2} \
          --outdir ${output_dir} -n ${GROUP_name}_MACS2_narrow
    echo "["`date`"]: " "Finished! "
}
## H3K27ac, H3K4me1 and H3K4me3 are narrow. H3K27me3 is broad
MACS2_broad_peak(){
    ## For H3K27me3 PolII ATAC FAIRE 
    ## Usage: MACS2_broad_peak [control bam file] [treatment bam file]
    echo "["`date`"]: " "FUNCTION[MACS2_broad_peak]: processing " ${1}
    macs2 callpeak --keep-dup all -f BAMPE \
          -B --SPMR -g 3.7e8 \
          --nomodel --extsize ${fragment_len} -q 0.5 \
          --broad --broad-cutoff 0.5 \
          --cutoff-analysis  \
          -t ${1} -c ${2} \
          --outdir ${output_dir} -n ${GROUP_name}_MACS2_broad
    echo "["`date`"]: " "Finished! "
}

MACS2_ATAC_peak(){
    ## For H3K27me3 PolII ATAC FAIRE 
    ## Usage: MACS2_broad_peak [control bam file] [treatment bam file]
    echo "["`date`"]: " "FUNCTION[MACS2_ATAC_peak]: processing " ${1}
    macs2 callpeak --keep-dup all -f BAMPE \
          -B --SPMR -g 3.7e8 \
          --nomodel --shift 100 --extsize 200 -q 0.5 \
          --call-summits --cutoff-analysis  \
          -t ${1} \
          --outdir ${output_dir} -n ${GROUP_name}_MACS2_ATAC
    echo "["`date`"]: " "Finished! "
}

HMMRATAC_peak(){
    ## For ATAC (only for paired data)
    ## Usage: HMMRATAC_peak [input bam file] 
    echo "["`date`"]: " "FUNCTION[HMMRATAC_peak]: processing " ${1}
    java -jar /opt/biosoft/HMMRATAC_V1.2.10_exe.jar -b ${1} \
         -i ${1}.bai -g ${chrom_size} \
         --removeDuplicates false \
         -o ${output_dir}/${GROUP_name}_HMMRATAC_peak   
    echo "["`date`"]: " "Finished! "
}

# HistoneHMM_peak(){
# to be continued
# }

FSeq2_peak(){
    ## Usage: FSeq2_peak [treat_bam] [control_bam] 
    ## For FAIRE-seq, ATAC-seq, DNase-seq, broad peak
    echo "["`date`"]: " "FUNCTION[FSeq2_peak]: processing " ${1}
    fseq2 callpeak ${1} -control_file ${2} \
        -cpus ${threads} -pe -chrom_size_file ${chrom_size} \
        -o ${output_dir} -name ${GROUP_name}_fseq2_peak \
        -sig_format bigwig -f 0 -l 600 -t 8.0 -q_thr 0.5  
    echo "["`date`"]: " "Finished! "
}
###############################################################
######################## MAIN FUNCTION ########################
grouplist=${1}  ## sample info csv file
output_dir=${2} 
threads=${3}
chrom_size="/share/data/reference/osa/chrom.sizes"
##
cat ${grouplist}|sed '1d'|\
while IFS="," read {ID_treat,ID_control,Seq_type,Peak_type}
do
    ## BAM files pathname
    treat_bam="3.align/bowtie2/${ID_treat}_rmdup.bam"
    control_bam="3.align/bowtie2/${ID_control}_rmdup.bam"
    GROUP_name="${ID_treat}_${ID_control}"
    echo "["`date`"]: " "Running Peaks calling of: " ${GROUP_name}
    ## BAM index 
    for bams in ${treat_bam} ${control_bam}
    do
        if [ ! -f "${bams}.bai" ];then
            echo "["`date`"]: " "BAM index of " ${bams}
            samtools index -@ ${threads} ${bams}
        fi
    done
    ## Get fragment length
    if [ ! -f "${treat_bam}_fragment_len" ];then 
        frglen_spp ${treat_bam}
        fragment_len=`cat ${treat_bam}_fragment_len`
        echo "["`date`"]: " "Setting fragment length as: " ${fragment_len}
    else
        echo "["`date`"]: " "Skip running FUNCTION[frglen_spp] due to file existed for: " ${treat_bam}
        fragment_len=`cat ${treat_bam}_fragment_len`
        echo "["`date`"]: " "Setting fragment length as: " ${fragment_len}
    fi
    ## Call peak depends on seq types 
    case ${Seq_type} in
    "ChIP") ## 
        echo "["`date`"]: Choosing seq type as ChIP-seq"
        if [ "${Peak_type}" == "narrow" ];then
            MACS2_narrow_peak ${treat_bam} ${control_bam}
        elif [ "${Peak_type}" == "broad" ];then
            MACS2_broad_peak ${treat_bam} ${control_bam}
        else
            echo "["`date`"]: " "ChIP-seq Peak_type Error!" ${GROUP_name}
        fi
    ;;
    "CutTag") ## To improve: use _sorted.bam not _rmdup.bam 
        echo "["`date`"]: Choosing seq type as Cut&Tag"
        if [ "${Peak_type}" == "narrow" ];then
            MACS2_narrow_peak ${treat_bam} ${control_bam}
        elif [ "${Peak_type}" == "broad" ];then
            MACS2_broad_peak ${treat_bam} ${control_bam}
        fi
    ;;
    "ATAC") ## broad_peak
        echo "["`date`"]: Choosing seq type as ATAC-seq"
        MACS2_ATAC_peak ${treat_bam}
        # HMMRATAC_peak ${treat_bam} ## too slow 
        # FSeq2_peak ${treat_bam} ${control_bam}
    ;;
    "FAIRE")  ## broad_peak
        echo "["`date`"]: Choosing seq type as FAIRE-seq"
        MACS2_ATAC_peak ${treat_bam} ${control_bam}
        # HMMRATAC_peak ${treat_bam} ## too slow 
        FSeq2_peak ${treat_bam} ${control_bam}
    ;;
    esac 
    echo "["`date`"]: " "Finished sample group : " ${GROUP_name}
done
echo "["`date`"]: " "Finished All Samples! "
