#!/usr/bin/env Rscript
## Usage:
##   Rscript run_gsea.R <FoldChange文件列表> <输出目录前缀> <species: osa|hsa> <osa OrgDb tarball>
## GSEA（gseGO）；物种与 OrgDb 路径由 enrich.sh 从 config 传入；依赖环境见 envs/enrich.yaml
args <- commandArgs(T)
pkgs <- c('clusterProfiler','ggplot2','enrichplot','dplyr')
lapply(pkgs, function(x){
   suppressMessages(library(x, character.only = T))})

spe <- "osa"
if (length(args) >= 3) spe <- args[3]
orgdb_tar <- NA
if (length(args) >= 4) orgdb_tar <- args[4]

if (spe == "hsa") {
    suppressMessages(library(org.Hs.eg.db))
    orgdb <- org.Hs.eg.db
    keytype <- "ENSEMBL"
} else {
    if (!require("org.Osativa.eg.db", quietly = TRUE)) {
        if (is.na(orgdb_tar)) stop("osa 需要 OrgDb：请在 config 的 orgdb_tarball 中提供本地 tarball 路径")
        install.packages(orgdb_tar, repos = NULL)
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
