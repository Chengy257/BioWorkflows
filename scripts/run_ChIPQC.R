#!/usr/bin/env Rscript
#########################################################################
# File Name: run_ChIPQC.R
# Author: ChengYu
# Description: 
# Created Time: Sat 02 Mar 2024 08:44:17 PM CST
#########################################################################
## Load libraries
pkgs <- c('ChIPQC','GenomicFeatures')
lapply(pkgs, function(x){
   suppressMessages(library(x, character.only = T))})

## Load sample data
samples <- read.csv('~/chipseq/results/chip_qc/ChIPQC/samplesheet.csv')

# The sample sheet contains metadata information for our dataset. Each row represents a peak set (which in most cases is every ChIP sample) and several columns of required information, which allows us to easily load the associated data in one single command.
# NOTE: The column headers have specific names that are expected by ChIPQC!!.
# SampleID: Identifier string for sample
# Tissue, Factor, Condition: Identifier strings for up to three different factors (You will need to have all columns listed. If you don't have infomation, then set values to NA)
# Replicate: Replicate number of sample
# bamReads: file path for BAM file containing aligned reads for ChIP sample
# ControlID: an identifier string for the control sample
# bamControl: file path for bam file containing aligned reads for control sample
# Peaks: path for file containing peaks for sample
# PeakCaller: Identifier string for peak caller used. Possible values include “raw”, “bed”, “narrow”, “macs”

Dir <- basename()

gtf_file <- args[1]
# gtf_file <- "/home/chengyu/references/osa/Oryza_sativa.IRGSP-1.0.56.Chr.gtf"
TxDb <- makeTxDbFromGFF(gtf_file)

## Create ChIPQC object
chipObj <- ChIPQC(samples, annotation=TxDb)

## Save the chipObj to file
save(chipObj, file="~/chipseq/results/chip_qc/ChIPQC/chipObj.RData")

## Create ChIPQC report
ChIPQCreport(chipObj, reportName="Nanog_and_Pou5f1", reportFolder="~/chipseq/results/chip_qc/ChIPQC/ChIPQCreport")

