#!/usr/bin/Rscript
## Usage:
## /usr/bin/Rscript /home/chengyu/workflows/enrich_GO_KEGG_clusterProfiler_gProfilerGO.R [gene list] [output dir prefix]     
.libPaths("/home/chengyu/R/Rlib_4.2.3")
.libPaths()
args <- commandArgs(T)
# suppressMessages(library(clusterProfiler))
# suppressMessages(library(ggplot2))
# suppressMessages(library(aPEAR))
# suppressMessages(library(svglite))
# if (!require("BiocManager", quietly = TRUE)){
#     install.packages("BiocManager")
#     BiocManager::install("pacman")
# }
# suppressMessages(library(pacman))
# pacman::p_load('clusterProfiler','ggplot2','aPEAR','svglite')
pkgs <- c('clusterProfiler','ggplot2','aPEAR','svglite','magrittr')
lapply(pkgs, function(x){
   suppressMessages(library(x, character.only = T))})


genelist <- read.table(args[1])[,1]
head(genelist)
output_dir <- args[2]
out <- paste(output_dir,"/",sep="")
prefix <- gsub(".DEGs.txt","",basename(args[1]))

## KEGG
R.utils::setOption("clusterProfiler.download.method",'auto')
T2G <- read.table("/share/data/reference/osa/GO/KEGG/KEGG_TranscriptID2GeneIDs",header=F,sep="\t")
## translist <- T2G[(which( genelist %in% T2G[,2])),1] 

tmp <- merge(x=T2G,y=genelist,by.x=2,by.y=1)
translist <- unique(tmp[,2])
head(translist)
ekegg <- enrichKEGG(translist, organism = "dosa",keyType = "kegg",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1)
newnames <- c()
for(i in 1:length(ekegg@result$Description)){
  name <- unlist(strsplit(ekegg@result$Description[i],split = " - Oryza sativa japonica")[[1]][1])
  newnames <- c(newnames,name)
}
ekegg@result$Description <- newnames
dotplot(ekegg,showCategory =10)+labs(title = "KEGG")


print(paste0("[",date(),"] ","KEGG enrich done!"))
write.table(as.data.frame(ekegg@result),file = paste0(out,prefix,"_EnrichResult_KEGG_",Sys.Date(),".xls"),quote = F,sep = "\t",col.names = T,row.names = F)
write.table(as.data.frame(ekegg@result %>% p.adjust <= 0.05)),file = paste0(out,prefix,"_EnrichResult_KEGG_significant_",Sys.Date(),".xls"),quote = F,sep = "\t",col.names = T,row.names = F)




## GO enrichment redundancy removal
## filter(pvalue <= 0.05 & p.adjust <= 0.3)   head(100)
aPEAR::enrichmentNetwork(ego_BP@result %>%  head(100),colorBy = 'p.adjust', colorType = 'pval', drawEllipses = FALSE,repelLabels = TRUE,verbose = F,nodeSize = "Count")
ggsave(filename = paste0(out,prefix,"_GO_BP_",Sys.Date(),"_top100Network.pdf"),device = "pdf",width = 6,height = 6)
aPEAR::enrichmentNetwork(ego_MF@result %>%  head(100),colorBy = 'p.adjust', colorType = 'pval', drawEllipses = FALSE,repelLabels = TRUE,verbose = F,nodeSize = "Count")
ggsave(filename = paste0(out,prefix,"_GO_MF_",Sys.Date(),"_top100Network.pdf"),device = "pdf",width = 6,height = 6)
aPEAR::enrichmentNetwork(ego_CC@result %>%  head(100),colorBy = 'p.adjust', colorType = 'pval', drawEllipses = FALSE,repelLabels = TRUE,verbose = F,nodeSize = "Count")
ggsave(filename = paste0(out,prefix,"_GO_CC_",Sys.Date(),"_top100Network.pdf"),device = "pdf",width = 6,height = 6)
## output aPEAR::findPathClusters results 
pathcl_BP <- aPEAR::findPathClusters(ego_BP@result %>%  head(100),verbose = F)
pathcl_MF <- aPEAR::findPathClusters(ego_MF@result %>%  head(100),verbose = F)
pathcl_CC <- aPEAR::findPathClusters(ego_CC@result %>%  head(100),verbose = F)
write.table(pathcl_BP$clusters,file = paste0(out,prefix,"_GO_BP_",Sys.Date(),"_top100Network_PathClusters.xls"),quote = F,sep = "\t",col.names = T,row.names = F)
write.table(pathcl_MF$clusters,file = paste0(out,prefix,"_GO_MF_",Sys.Date(),"_top100Network_PathClusters.xls"),quote = F,sep = "\t",col.names = T,row.names = F)
write.table(pathcl_CC$clusters,file = paste0(out,prefix,"_GO_CC_",Sys.Date(),"_top100Network_PathClusters.xls"),quote = F,sep = "\t",col.names = T,row.names = F)

print(paste0("[",date(),"] ","plot done!"))
print(paste0("[",date(),"] ","All Finished!"))

