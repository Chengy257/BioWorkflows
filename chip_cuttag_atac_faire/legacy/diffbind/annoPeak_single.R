#!/usr/bin/env Rscript
# =====================================================================
# ChIPseeker 单峰文件注释（独立小工具，不在主流程 DAG 中）
# 用法: Rscript annoPeak_single.R <gtf_file> <peak_file>
# 输出: <peak_file>_peakAnnoTable.xls 与 <peak_file>_peakAnnoPlot.pdf
# =====================================================================

pkgs <- c("ChIPseeker", "GenomicFeatures")
lapply(pkgs, function(x) suppressMessages(library(x, character.only = TRUE)))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript annoPeak_single.R <gtf_file> <peak_file>")
}
gtf_file   <- args[1]
peak_files <- args[2]

my_txdb <- makeTxDbFromGFF(gtf_file)
peaks <- readPeakFile(peak_files)
peakAnno <- annotatePeak(peaks, TxDb = my_txdb, level = "gene", verbose = FALSE,
                         addFlankGeneInfo = TRUE, flankDistance = 3000)
write.table(as.data.frame(peakAnno), paste0(peak_files, "_peakAnnoTable.xls"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

TSSRegion <- getPromoters(TxDb = my_txdb, upstream = 3000, downstream = 3000)
tagMatrix <- getTagMatrix(peaks, windows = TSSRegion)
pdf(paste0(peak_files, "_peakAnnoPlot.pdf"), height = 4, width = 6)
  plotAnnoBar(peakAnno)
  plotAnnoPie(peakAnno)
  vennpie(peakAnno)
  plotDistToTSS(peakAnno, title = "Distribution of peak relative to TSS")
  plotAvgProf(tagMatrix, xlim = c(-3000, 3000), conf = 0.95, resample = 500, facet = "row")
dev.off()
