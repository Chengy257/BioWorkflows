#!/usr/bin/env Rscript
## GO/KEGG over-representation for one DEG gene list (clusterProfiler).
##
## Ported from the legacy rna-seq workflow/scripts/run_enrichment.R with the
## species switch removed (decision record D3: the pipeline never binds a
## species - OrgDb, keyType and KEGG code are all arguments). aPEAR stays
## optional and is skipped silently when absent (dead upstream - never a
## hard requirement).
##
## Usage:
##   Rscript run_enrichment.R <gene list txt> <output dir> <contrast direction tag> \
##       <OrgDb package name> <keyType> <KEGG organism code|''> [<OrgDb tarball path>]
##
## KEGG queries rely on the KEGG REST API (internet access); on failure KEGG
## is skipped and GO results are unaffected (documented offline behavior).
args <- commandArgs(T)
if (length(args) < 6) {
    stop("usage: run_enrichment.R <gene list> <outdir> <direction tag> <OrgDb pkg> <keyType> <KEGG code|''> [<OrgDb tarball>]")
}

suppressMessages(library(clusterProfiler))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

genelist_file <- args[1]
output_dir <- args[2]
direction <- args[3]
orgdb_name <- args[4]
keytype <- args[5]
kegg_org <- args[6]
orgdb_tarball <- if (length(args) >= 7 && args[7] != "") args[7] else NULL

## Container route: install the staged OrgDb tarball into a task-local
## library first; conda/dev route: the package must already be installed.
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

genelist <- read.table(genelist_file)[, 1]
genelist <- genelist[genelist != ""]
if (length(genelist) == 0) {
    message("[", direction, "] empty gene list - nothing to enrich")
    quit(save = "no", status = 0)
}

out <- paste0(output_dir, "/", sep = "")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

for (ont in c("MF", "BP", "CC", "KEGG")) {
    tryCatch({
        if (ont == "KEGG") {
            if (kegg_org == "") {
                message("[WARN] KEGG skipped: kegg_organism param is empty")
                return()
            }
            ego <- enrichKEGG(genelist, organism = kegg_org, keyType = "kegg",
                              pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1)
            outname_prefix <- "KEGG"
        } else {
            ego <- enrichGO(genelist, OrgDb = orgdb, ont = ont, keyType = keytype,
                            pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH")
            outname_prefix <- paste0("GO-", ont)
        }
        write.table(as.data.frame(ego@result),
                    file = paste0(out, direction, "_", outname_prefix, "_AllResults.tsv"),
                    quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
        write.table(as.data.frame(ego@result %>% filter(p.adjust <= 0.05)),
                    file = paste0(out, direction, "_", outname_prefix, "_PadjSignificant.tsv"),
                    quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
        if (nrow(as.data.frame(ego@result)) > 0) {
            dotplot(ego, showCategory = 10) + labs(title = paste0(outname_prefix, " top10 Terms"))
            ggsave(filename = paste0(out, direction, "_", outname_prefix, "_dotplot.pdf"),
                   device = "pdf", width = 5, height = 5)
        }
        if (requireNamespace("aPEAR", quietly = TRUE)) {
            tryCatch({
                aPEAR::enrichmentNetwork(ego@result %>% filter(pvalue <= 0.05),
                                         colorBy = 'p.adjust', colorType = 'pval',
                                         drawEllipses = FALSE, repelLabels = TRUE,
                                         verbose = FALSE, nodeSize = "Count")
                ggsave(filename = paste0(out, direction, "_", outname_prefix, "_sigNetwork.pdf"),
                       device = "pdf", width = 6, height = 6)
            }, error = function(e) {
                message("[WARN] aPEAR network skipped for ", outname_prefix, ": ", conditionMessage(e))
            })
        }
    }, error = function(e) {
        message("[WARN] ", ont, " enrichment skipped: ", conditionMessage(e))
    })
}
message("[", date(), "] enrichment finished: ", direction)
