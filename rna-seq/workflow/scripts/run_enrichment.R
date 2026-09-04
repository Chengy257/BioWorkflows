#!/usr/bin/env Rscript
## Usage:
##   Rscript run_enrichment.R <gene list> <output dir> <species: osa|hsa> <osa OrgDb tarball> <kegg organism code>
## 依赖环境见 envs/enrich.yaml；物种与 OrgDb 路径由 enrich.sh 从 config 传入
args <- commandArgs(T)
pkgs <- c('clusterProfiler','ggplot2','aPEAR','svglite','magrittr','dplyr')
lapply(pkgs, function(x){
   suppressMessages(library(x, character.only = T))})

spe <- "osa"
if (length(args) >= 3) spe <- args[3]
orgdb_tar <- NA
if (length(args) >= 4) orgdb_tar <- args[4]
kegg_org <- ifelse(spe == "hsa", "hsa", "dosa")   # 默认按物种推导，可被第 5 参数覆盖
if (length(args) >= 5 && !is.na(args[5]) && args[5] != "") kegg_org <- args[5]

if (spe == "hsa") {
    suppressMessages(library(org.Hs.eg.db))
    orgdb <- org.Hs.eg.db
    keytype <- "ENSEMBL"
} else {
    if (!require("org.Osativa.eg.db", quietly = TRUE)) {
        if (is.na(orgdb_tar) || orgdb_tar == "")
            stop("osa 需要 OrgDb：请在 config 的 orgdb_tarball 中提供本地 tarball 路径")
        install.packages(orgdb_tar, repos = NULL)
    }
    suppressMessages(library(org.Osativa.eg.db))
    orgdb <- org.Osativa.eg.db
    keytype <- "GID"
}

genelist <- read.table(args[1])[,1]
if (spe == "hsa") genelist <- sub("\\..*$", "", genelist)   ## 去掉 ENSEMBL 版本号

output_dir <- args[2]
out <- paste(output_dir, "/", sep = "")
prefix <- basename(args[1])

###############################################
## KEGG 查询依赖 KEGG REST API（需外网）；失败时跳过 KEGG，不影响 GO 结果
R.utils::setOption("clusterProfiler.download.method", "auto")
###############################################
for (ont in c("MF", "BP", "CC", "KEGG")) {
    tryCatch({
        if (ont == "KEGG") {
            if (spe == "hsa") {
                eg <- bitr(genelist, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = orgdb)
                ego <- enrichKEGG(eg$ENTREZID, organism = kegg_org, keyType = "kegg",
                                  pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1)
            } else {
                ego <- enrichKEGG(genelist, organism = kegg_org, keyType = "kegg",
                                  pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1)
            }
            outname_prefix <- "KEGG"
            if (kegg_org == "dosa") {
                newnames <- c()
                for (i in 1:length(ego@result$Description)) {
                    name <- unlist(strsplit(ego@result$Description[i], split = " - Oryza sativa japonica"))[[1]][1]
                    newnames <- c(newnames, name)
                }
                ego@result$Description <- newnames
            }
        } else {
            ego <- enrichGO(genelist, OrgDb = orgdb, ont = ont, keyType = keytype,
                            pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH")
            outname_prefix <- paste0("GO:", ont)
        }
        write.table(as.data.frame(ego@result), file = paste0(out, prefix, "_", outname_prefix, "_AllResults_", Sys.Date(), ".xls"), quote = F, sep = "\t", col.names = T, row.names = F)
        write.table(as.data.frame(ego@result %>% filter(p.adjust <= 0.05)), file = paste0(out, prefix, "_", outname_prefix, "_PadjSignificant_", Sys.Date(), ".xls"), quote = F, sep = "\t", col.names = T, row.names = F)

        dotplot(ego, showCategory = 10) + labs(title = paste0(outname_prefix, " top10 Terms"))
        ggsave(filename = paste0(out, prefix, "_", outname_prefix, "_", Sys.Date(), "_dotplot.pdf"), device = "pdf", width = 5, height = 5)

        tryCatch({
            aPEAR::enrichmentNetwork(ego@result %>% filter(pvalue <= 0.05), colorBy = 'p.adjust', colorType = 'pval', drawEllipses = FALSE, repelLabels = TRUE, verbose = F, nodeSize = "Count") + labs(title = paste0(outname_prefix, "top 100 Terms Network"))
            ggsave(filename = paste0(out, prefix, "_", outname_prefix, "_", Sys.Date(), "_sigNetwork.pdf"), device = "pdf", width = 6, height = 6)
            pathcl <- aPEAR::findPathClusters(ego@result %>% filter(pvalue <= 0.05), verbose = F)
            write.table(pathcl$clusters, file = paste0(out, prefix, "_", outname_prefix, "_", Sys.Date(), "_sigNetwork_PathClusters.xls"), quote = F, sep = "\t", col.names = T, row.names = F)
        }, error = function(e) {
            print(paste0("Error while running aPEAR!", outname_prefix, conditionMessage(e)))
        })
    }, error = function(e) {
        print(paste0("[WARN] ", ont, " enrichment skipped: ", conditionMessage(e)))
    })
}
print(paste0("[", date(), "] ", "All Finished!"))
