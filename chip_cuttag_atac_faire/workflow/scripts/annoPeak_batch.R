#!/usr/bin/env Rscript
# =====================================================================
# ChIPseeker batch peak annotation and distribution plots
# Usage: Rscript annoPeak_batch.R <gtf_file> <peak_files comma-separated> <outdir> [flank_bp]
# Output: <outdir>/<sample>.Anno.xls and Peakanno_PeakDistributions.pdf
# =====================================================================

pkgs <- c("ChIPseeker", "GenomicFeatures")
lapply(pkgs, function(x) suppressMessages(library(x, character.only = TRUE)))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript annoPeak_batch.R <gtf_file> <peak_files comma-separated> <outdir> [flank_bp]")
}
gtf_file   <- args[1]
peak_files <- strsplit(args[2], ",", fixed = TRUE)[[1]]
outdir     <- args[3]
flank      <- if (length(args) >= 4) as.integer(args[4]) else 3000L

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Drop missing or empty files (peak calling may produce an empty set)
peak_files <- peak_files[file.exists(peak_files)]
peak_files <- peak_files[file.info(peak_files)$size > 0]
if (length(peak_files) == 0) {
  stop("No valid peak files to annotate (all missing or empty)")
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
if (length(peaks) == 0) stop("All peak files are empty")

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
