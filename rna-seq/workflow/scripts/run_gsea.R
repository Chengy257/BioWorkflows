#!/usr/bin/env Rscript
## Usage:
##   Rscript run_gsea.R <FoldChange文件列表> <输出目录前缀> <species: osa|hsa>
## GSEA（gseGO）；物种与 OrgDb 路径由 enrich.sh 从 config 传入；R runtime and library paths are provided by config/software.yaml
args <- commandArgs(T)
pkgs <- c('clusterProfiler','ggplot2','enrichplot','dplyr')
lapply(pkgs, function(x){
   suppressMessages(library(x, character.only = T))})

spe <- "osa"
if (length(args) >= 3) spe <- args[3]

if (spe == "hsa") {
    suppressMessages(library(org.Hs.eg.db))
    orgdb <- org.Hs.eg.db
    keytype <- "ENSEMBL"
} else {
    if (!requireNamespace("org.Osativa.eg.db", quietly = TRUE)) {
        stop("org.Osativa.eg.db is not installed in the configured R libraries; install it before running the workflow")
    }
    suppressMessages(library(org.Osativa.eg.db))
    orgdb <- org.Osativa.eg.db
    keytype <- "GID"
}

filenamelist <- read.table(args[1])
print(filenamelist)
out_prefix <- as.character(args[2])

if (dir.exists(out_prefix) == FALSE) {
    system(paste0("mkdir -p ", out_prefix))
}

ONT <- c("BP", "MF", "CC")
for (i in 1:length(filenamelist[,1])) {
    filename <- filenamelist[i,1]
    print(paste0("parsing ", filename, " ..."))
    df <- read.table(filename, header = T)
    groupname <- strsplit(basename(filename), split = "_FoldChange")[[1]][1]

    fc <- as.vector(df$log2FoldChange[order(df$log2FoldChange, decreasing = T)])
    gids <- df$gene_id[order(df$log2FoldChange, decreasing = T)]
    if (spe == "hsa") gids <- sub("\\..*$", "", gids)   ## 去掉 ENSEMBL 版本号
    names(fc) <- gids

    for (j in 1:3) {
        tryCatch({
            res <- clusterProfiler::gseGO(geneList = fc, OrgDb = orgdb, ont = ONT[j], keyType = keytype, pvalueCutoff = 1)
            outname <- paste0(out_prefix, "/", groupname, "_", ONT[j], "_multiGSEA_", Sys.Date())
            write.table(as.data.frame(res@result), file = paste0(outname, "_ResultTable_all.xls"), col.name = T, row.name = F, sep = "\t", quote = F)
            write.table(as.data.frame(res@result %>% filter(abs(NES) >= 1 & pvalue <= 0.05 & qvalue <= 0.25)), file = paste0(outname, "_ResultTable_significant.xls"), col.name = T, row.name = F, sep = "\t", quote = F)
            dotplot(res, showCategory = 20) + facet_grid(~.sign) + labs(title = paste0(groupname, " GO:", ONT[j]))
            ggsave(file = paste0(outname, "_multiGSEA_dotplot_", Sys.Date(), ".pdf"), device = "pdf")
            print(paste0("Finished ", groupname, " GO:", ONT[j], "..."))
        }, error = function(e) {
            print(paste0("Error processing ", groupname, " GO:", ONT[j], ": ", conditionMessage(e)))
        })
    }
}

print("All done!")
