#!/usr/bin/env Rscript
## GSEA (gseGO) on a full ranked gene list (clusterProfiler).
##
## Ported from the legacy rna-seq workflow/scripts/run_gsea.R core with the
## species switch removed (decision record D3): OrgDb package name and
## keyType are arguments; the ranked table carries gene_id + the ranking
## statistic (stat by default, log2FoldChange fallback - see
## prepare_enrichment_inputs.py).
##
## Usage:
##   Rscript run_gsea.R <ranked tsv> <output dir> <contrast tag> \
##       <OrgDb package name> <keyType> [<OrgDb tarball path>]
args <- commandArgs(T)
if (length(args) < 5) {
    stop("usage: run_gsea.R <ranked tsv> <outdir> <contrast tag> <OrgDb pkg> <keyType> [<OrgDb tarball>]")
}

suppressMessages(library(clusterProfiler))
suppressMessages(library(ggplot2))
suppressMessages(library(enrichplot))
suppressMessages(library(dplyr))

ranked_file <- args[1]
output_dir <- args[2]
contrast <- args[3]
orgdb_name <- args[4]
keytype <- args[5]
orgdb_tarball <- if (length(args) >= 6 && args[6] != "") args[6] else NULL

if (!is.null(orgdb_tarball)) {
    lib_dir <- file.path(getwd(), "Rlib")
    dir.create(lib_dir, showWarnings = FALSE, recursive = TRUE)
    install.packages(orgdb_tarball, lib = lib_dir, repos = NULL, type = "source")
    .libPaths(c(lib_dir, .libPaths()))
}
if (!requireNamespace(orgdb_name, quietly = TRUE)) {
    stop(orgdb_name, " is not installed; pass the OrgDb tarball via params.enrichment.orgdb_file")
}
suppressMessages(library(orgdb_name, character.only = TRUE))
orgdb <- getExportedValue(orgdb_name, orgdb_name)

df <- read.table(ranked_file, header = TRUE)
ordered <- df[order(df$stat, decreasing = TRUE), ]
fc <- as.vector(ordered$stat)
names(fc) <- as.character(ordered$gene_id)
fc <- fc[!duplicated(names(fc)) & !is.na(fc)]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
for (ont in c("BP", "MF", "CC")) {
    tryCatch({
        res <- clusterProfiler::gseGO(geneList = fc, OrgDb = orgdb, ont = ont,
                                      keyType = keytype, pvalueCutoff = 1)
        outname <- file.path(output_dir, paste0(contrast, "_", ont, "_multiGSEA"))
        write.table(as.data.frame(res@result),
                    file = paste0(outname, "_ResultTable_all.tsv"),
                    col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
        sig <- as.data.frame(res@result %>% filter(abs(NES) >= 1 & pvalue <= 0.05 & qvalue <= 0.25))
        write.table(sig, file = paste0(outname, "_ResultTable_significant.tsv"),
                    col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
        if (nrow(as.data.frame(res@result)) > 0) {
            dotplot(res, showCategory = 20) + facet_grid(~.sign) +
                labs(title = paste0(contrast, " GO:", ont))
            ggsave(file = paste0(outname, "_dotplot.pdf"), device = "pdf")
        }
    }, error = function(e) {
        message("[WARN] GSEA ", ont, " skipped for ", contrast, ": ", conditionMessage(e))
    })
}
message("[", date(), "] GSEA finished: ", contrast)
