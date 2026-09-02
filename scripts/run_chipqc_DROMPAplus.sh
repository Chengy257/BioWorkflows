#!/bin/bash
DIR=${1}
QC_samples=${2}
cd ${DIR} && mkdir -p ${DIR}/tempdir_runDROMPAplus
# cp ${BAM} ${DIR}/tempdir_runDROMPAplus/
# cp ${PEAK_BED} ${DIR}/tempdir_runDROMPAplus/
cp /share/data/reference/osa/chrom.sizes ${DIR}/tempdir_runDROMPAplus/
##
cat ${QC_samples}|while read {bam,peak} 
do
    name=`basename $bam|cut -d_ -f1`
    docker run --rm -v `pwd`:/mnt rnakato/ssp_drompa parse2wig+ \
        -i /mnt/${bam} \
        -o ${name} --odir /mnt/5.QC/DROMPAplus/ \
        --gt /mnt/tempdir_runDROMPAplus/chrom.sizes --pair \
        --bed /mnt/${peak}
done
# docker run  --rm -v ${DIR}:/mnt rnakato/ssp_drompa parse2wig+ -i /mnt/3.align/bowtie2/A-1_rmdup.bam -o A-1 --odir /mnt/5.QC/DROMPAplus/ --gt /mnt/tempdir_runDROMPAplus/chrom.sizes --pair --bed /mnt/4.peaks/A-1_in-1_MACS2_ATAC_peaks.narrowPeak
