#!/usr/bin/env Rscript
# Differential binding analysis for chip_cuttag_atac_faire (DiffBind, v0.5).
# Reads one DiffBind sample sheet built by scripts/diffbind_sheet.py
# (SampleID/Condition/Replicate/bamReads/bamControl/Peaks/PeakCaller/Batch), counts
# reads over the consensus peak set (summit-recentered when summit_flank > 0),
# fits the contrast with DESeq2 or edgeR (optionally blocked on the Batch
# column), and writes into the sheet's directory:
#   DB_results.tsv            every consensus region with Fold / FDR / p-value
#   DB_significant.tsv        subset passing the FDR (and |Fold| when > 1)
#   MA_plot.png / Volcano_plot.png / PCA_plot.png
#   sessionInfo.txt
# R runtime and library paths are provided by config/software.yaml.
suppressMessages(library(DiffBind))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop(paste("Usage: run_diffbind.R <samplesheet.tsv> <DESeq2|edgeR>",
             "<summit_flank> <fdr> <foldchange> <block T|F>"))
}
sheet_file <- args[1]
analysis <- toupper(args[2])
summit_flank <- as.integer(args[3])
fdr_th <- as.numeric(args[4])
fc_th <- as.numeric(args[5])
use_block <- toupper(args[6]) == "T"
outdir <- dirname(sheet_file)

sheet <- read.delim(sheet_file, header = TRUE, colClasses = "character",
                    check.names = FALSE)
sheet$Replicate <- as.integer(sheet$Replicate)
sheet$Condition <- as.factor(sheet$Condition)
if (analysis == "DESEQ2") {
  method <- DBA_DESEQ2
} else {
  method <- DBA_EDGER
}
# The optional Batch column stays in the sheet (DiffBind turns extra factor
# columns into attributes); it becomes a blocking factor only when the user
# asked for batch correction and it carries >= 2 levels.
batch_levels <- character(0)
if ("Batch" %in% colnames(sheet)) {
  batch_levels <- unique(na.omit(sheet$Batch[sheet$Batch != ""]))
}
has_batch <- use_block && length(batch_levels) >= 2
if (use_block && !has_batch) {
  message("[diffbind] batch correction requested but the Batch column is ",
          "absent or has fewer than two levels; running without blocking")
}

# Peak files and their caller (narrowpeak|broadpeak, encoding the score
# column convention) arrive in the sheet's Peaks/PeakCaller columns —
# DiffBind's own sample-sheet names (WSL real-run finding 2026-09-10: the
# PeakFormat/ScoreCol overrides were vignette-era and a wrong column name
# silently yields an empty peak set).

dba <- dba(sampleSheet = sheet)
message("[diffbind] consensus peak set: ", dba$numTotal, " regions in ",
        nrow(sheet), " samples")

count_args <- list(DBA = dba, bUseSummarizeOverlaps = FALSE)
if (summit_flank > 0) count_args$summits <- summit_flank
dba <- do.call(dba.count, count_args)

conds <- as.character(unique(sheet$Condition))
# dba.contrast's block argument takes a DBA_ attribute constant or a logical
# vector aligned with the sample sheet (a character vector of batch labels is
# the old API); block on the first batch level = the classic paired design
block_mask <- if (has_batch) sheet$Batch == batch_levels[1] else NULL
dba <- dba.contrast(dba,
                    group1 = dba$masks[[conds[1]]],
                    group2 = dba$masks[[conds[2]]],
                    block = block_mask)
dba <- dba.analyze(dba, method = method)

res <- dba.report(dba, method = method, contrast = 1,
                  th = 1, bUsePval = FALSE, DataType = DBA_DATA_FRAME)
colnames(res)[1] <- "chr"
write.table(res, file.path(outdir, "DB_results.tsv"), sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

sig <- subset(res, FDR <= fdr_th & abs(Fold) >= fc_th)
write.table(sig, file.path(outdir, "DB_significant.tsv"), sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
message("[diffbind] ", nrow(sig), " / ", nrow(res), " regions pass FDR<=",
        fdr_th, " |Fold|>=", fc_th)

png(file.path(outdir, "MA_plot.png"), width = 1200, height = 900, res = 130)
dba.plotMA(dba, method = method, contrast = 1, th = fdr_th)
dev.off()
png(file.path(outdir, "Volcano_plot.png"), width = 1200, height = 900, res = 130)
dba.plotVolcano(dba, method = method, contrast = 1, th = fdr_th)
dev.off()
png(file.path(outdir, "PCA_plot.png"), width = 1200, height = 900, res = 130)
dba.plotPCA(dba, attributes = DBA_CONDITION, label = DBA_REPLICATE)
dev.off()

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
message("[diffbind] done: ", outdir)
