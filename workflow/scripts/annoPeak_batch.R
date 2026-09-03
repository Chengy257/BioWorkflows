#!/usr/bin/env Rscript
# =====================================================================
# ChIPseeker 批量峰注释与分布图
# 用法: Rscript annoPeak_batch.R <gtf_file> <peak_files逗号分隔> <outdir> [flank_bp]
# 输出: <outdir>/<sample>.Anno.xls 与 Peakanno_PeakDistributions.pdf
# =====================================================================

pkgs <- c("ChIPseeker", "GenomicFeatures")
lapply(pkgs, function(x) suppressMessages(library(x, character.only = TRUE)))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript annoPeak_batch.R <gtf_file> <peak_files逗号分隔> <outdir> [flank_bp]")
}
gtf_file   <- args[1]
peak_files <- strsplit(args[2], ",", fixed = TRUE)[[1]]
outdir     <- args[3]
flank      <- if (length(args) >= 4) as.integer(args[4]) else 3000L

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 过滤不存在或空文件（峰调用结果可能为空集）
peak_files <- peak_files[file.exists(peak_files)]
peak_files <- peak_files[file.info(peak_files)$size > 0]
if (length(peak_files) == 0) {
  stop("没有有效的峰文件可注释（全部缺失或为空）")
}

load_file <- function(paths) {
  peaks <- list()
  for (f in paths) {
    sample <- sub("_peaks\\.(narrow|broad)Peak.*", "", basename(f))
    message("Reading: ", f)
    peaks[[sample]] <- readPeakFile(f)
  }
  peaks
}

my_txdb <- makeTxDbFromGFF(gtf_file)
peaks <- load_file(peak_files)
if (length(peaks) == 0) stop("所有峰文件均为空集")

peakAnnoList <- lapply(peaks, annotatePeak, TxDb = my_txdb, level = "gene",
                       verbose = FALSE, addFlankGeneInfo = TRUE, flankDistance = flank)

lapply(names(peakAnnoList), function(x) {
  file_name <- file.path(outdir, paste0(x, ".Anno.xls"))
  write.table(as.data.frame(peakAnnoList[[x]]), file_name,
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
})

TSSRegion <- getPromoters(TxDb = my_txdb, upstream = flank, downstream = flank)
tagMatrixList <- lapply(peaks, getTagMatrix, windows = TSSRegion)

pdf(file.path(outdir, "Peakanno_PeakDistributions.pdf"), height = 8, width = 8)
  plotAnnoBar(peakAnnoList)
  plotDistToTSS(peakAnnoList, title = "Distribution of binding loci relative to TSS")
  plotAvgProf(tagMatrixList, xlim = c(-flank, flank), conf = 0.95, resample = 500, facet = "row")
  tagHeatmap(tagMatrixList)
dev.off()

message("Peak annotation done: ", outdir)
