#!/usr/bin/env Rscript
# =====================================================================
# ChIPQC 质控报告（独立工具，不在主流程 DAG 中）
# 用法: Rscript run_ChIPQC.R <samplesheet.csv> <gtf_file> <outdir> [report_name]
#
# samplesheet 列（ChIPQC 约定，参见 Bioconductor ChIPQC 文档）:
#   SampleID, Tissue, Factor, Condition, Replicate,
#   bamReads, ControlID, bamControl, Peaks, PeakCaller("bed"/"narrow"/"macs")
# =====================================================================

pkgs <- c("ChIPQC", "GenomicFeatures")
lapply(pkgs, function(x) suppressMessages(library(x, character.only = TRUE)))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_ChIPQC.R <samplesheet.csv> <gtf_file> <outdir> [report_name]")
}
samplesheet <- args[1]
gtf_file    <- args[2]
outdir      <- args[3]
report_name <- ifelse(length(args) >= 4, args[4], "ChIPQC_report")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

samples <- read.csv(samplesheet)
message("样本表: ", nrow(samples), " 行")

# 峰路径与 BAM 相对工作目录解析（ChIPQC 以当前工作目录为基准）
TxDb <- makeTxDbFromGFF(gtf_file)

chipObj <- ChIPQC(samples, annotation = TxDb, chromosomes = NULL)
save(chipObj, file = file.path(outdir, "chipObj.RData"))

ChIPQCreport(chipObj, reportName = report_name,
             reportFolder = file.path(outdir, paste0(report_name, "_report")))

message("ChIPQC 报告完成: ", outdir)
