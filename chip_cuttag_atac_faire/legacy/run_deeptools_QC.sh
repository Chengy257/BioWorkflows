#!/bin/bash
#########################################################################
# File Name: scripts/QC_deeptools.sh
# Author: ChengYu
# Description: 
# Created Time: Wed 28 Feb 2024 10:34:21 PM CST
#########################################################################


ids=`ls 3.align/bowtie2/*rmdup.bw |xargs -i basename {{}}|cut -d_ -f1|tr "\n" " "`

## computeMatrix 
computeMatrix scale-regions -R {input.bed} -S `ls 3.align/bowtie2/*_rmdup.bw` -b 3000 -a 3000 -o {output.matrix} --outFileNameMatrix 5.QC_deeptools/matrix_scaled.tab -p {config[threads]} >> {log} 2>&1

## plotProfile
plotProfile -m {output.matrix} -out 5.QC_deeptools/plotprofiler.png --plotTitle "plotProfile" >> {log} 2>&1

## multiBamSummary
multiBamSummary bins -p {config[threads]} --bamfiles 3.align/bowtie2/*_rmdup.bam --minMappingQuality 30 --labels ${{ids}} -out 5.QC_deeptools/multiBamSummary.npz --outRawCounts 5.QC_deeptools/readCounts.tab >> {log} 2>&1

## plotCorrelation
plotCorrelation -in 5.QC_deeptools/multiBamSummary.npz --corMethod spearman --skipZeros --plotTitle  "Spearman Correlation of Read Counts" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o 5.QC_deeptools/heatmap_SpearmanCorr_readCounts.png  --outFileCorMatrix 5.QC_deeptools/SpearmanCorr_readCounts.tab >> {log} 2>&1

## PCA
plotPCA -in 5.QC_deeptools/multiBamSummary.npz -o 5.QC_deeptools/PCA_readCounts.png -T "PCA of read counts" >> {log} 2>&1

## PEFragmentSize
bamPEFragmentSize -p {config[threads]} -hist 5.QC_deeptools/fragmentsize.png -T "Fragment size of PE data" --maxFragmentLength 1000 -b 3.align/bowtie2/*_rmdup.bam --samplesLabel ${{ids}}  >> {log} 2>&1

## plotFingerprint
plotFingerprint -b 3.align/bowtie2/*_rmdup.bam --labels ${{ids}} --minMappingQuality 30 --skipZeros -T "Fingerprints of different samples" --plotFile 5.QC_deeptools/fingerprints.png --outRawCounts 5.QC_deeptools/fingerprints.tab
