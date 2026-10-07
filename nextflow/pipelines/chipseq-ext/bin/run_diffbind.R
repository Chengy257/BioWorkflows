#!/usr/bin/env Rscript
# Differential binding analysis with DiffBind for the chipseq-ext route.
# Reads the DiffBind sample sheet built by diffbind_sheet.py (SampleID /
# Condition / Replicate / bamReads / bamControl / Peaks / PeakCaller / Batch),
# counts reads over the per-sample peak sets (summit-recentered when
# summit_flank > 0), fits the contrast with DESeq2 or edgeR (optionally
# blocked on the Batch column), and writes into the working directory:
#   DB_results.tsv            every region with Fold / FDR / p-value
#   DB_significant.tsv        subset passing the FDR (and |Fold| when > 1)
#   MA_plot.png / Volcano_plot.png / PCA_plot.png
#   sessionInfo.txt
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
outdir <- getwd()

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
# column convention) arrive in the sheet's Peaks/PeakCaller columns -
# DiffBind's own sample-sheet names.

dba <- dba(sampleSheet = sheet)
message("[diffbind] peak set: ", dba$numTotal, " regions in ",
        nrow(sheet), " samples")

count_args <- list(DBA = dba, bUseSummarizeOverlaps = FALSE)
if (summit_flank > 0) count_args$summits <- summit_flank
dba <- do.call(dba.count, count_args)

conds <- as.character(unique(sheet$Condition))
# dba.contrast's block argument must stay ABSENT when there is no batch
# factor: passing block=NULL trips the attribute validation ("attribute must
# be a DBA_ attribute, a logical vector, or a list of logical vectors") - a
# logical vector over the sheet (first batch level = classic paired design)
# is only attached when blocking is actually possible.
contrast_args <- list(DBA = dba,
                      group1 = dba$masks[[conds[1]]],
                      group2 = dba$masks[[conds[2]]])
if (has_batch) {
  contrast_args$block <- sheet$Batch == batch_levels[1]
}
dba <- do.call(dba.contrast, contrast_args)
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
