#!/usr/bin/Rscript
args <- commandArgs(T)
.libPaths("/home/chengyu/R/Rlib_4.2.3")
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("pacman")
# suppressMessages(library(pacman))
# pacman::p_load('clusterProfiler','ggplot2','aPEAR','svglite','VennDiagram','UpSetR')

pkgs <- c('clusterProfiler','ggplot2','aPEAR','svglite','VennDiagram','UpSetR','magrittr')
lapply(pkgs, function(x){
   suppressMessages(library(x, character.only = T))})


if (!require("org.Osativa.eg.db", quietly = TRUE))
    install.packages("/home/chengyu/references/osa/GO/gProfiler/org.Osativa.eg.db_0.1.tar.gz",repos=NULL)
suppressMessages(library(org.Osativa.eg.db))
T2G <- read.table("/share/data/reference/osa/GO/KEGG/KEGG_TranscriptID2GeneIDs",header=F,sep="\t")

# geneID2TrID <- function(genelist){
#   tmp <- merge(x=T2G,y=genelist,by.x=2,by.y=1)
#   translist <- unique(tmp[,2])
#   return(translist)
# }


## read in multi gene set as a list 
getGeneSetList <- function(file_list){
  geneSets <- list()
  for(i in 1:length(file_list[,1])){
    tmp <- read.table(file_list[i,1],header = F,sep=" ")

    name <- file_list[i,2]
    # tr <- geneID2TrID(tmp[,1])

    # geneSets_Tr[[name]] <- tr
    geneSets[[name]] <- tmp[,1]
    
  }
  # res["GO"] <- geneSets
  # res["KEGG"] <- geneSets_Tr

  return(geneSets)
}

## get intersections 
getIntersect <-  function(list,name){
  newlist <- list()
  inter <- get.venn.partitions(list);
  venn.diagram(list,filename = paste0(name,"_venn_",Sys.Date(),".png"),imagetype = 'png');  ## plot venn diagram  
  png(paste0(name,"_upset_",Sys.Date(),".png"))
    upset(fromList(list),order.by = "freq")
  dev.off()
  # ggsave(filename = paste0(name,"_upset_",Sys.Date(),".pdf"),device = "pdf",width = 8,height = 8)
  for (i in 1:nrow(inter)) {
    newlist[[inter[i,'..set..']]] <- (inter[i,'..values..'])[[1]];
    inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ');
  } 
  write.table(inter[,!colnames(inter) %in% c('..values..')],paste(name,"_intersetion.xls",sep=""),quote = F,col.names = T,row.names = F,sep = "\t")
  ## only for TWO sets compare, rename the 3 parts! 
  names(newlist)[1] <- "intersections"
  names(newlist)[2] <- paste0(gsub("\\(","",strsplit(names(newlist)[2],split = ")")[[1]][1]),"_complementarySet")
  names(newlist)[3] <- paste0(gsub("\\(","",strsplit(names(newlist)[3],split = ")")[[1]][1]),"_complementarySet")  
  return(newlist)
}

## compare GO enrich 
compare <- function(geneSets,name){
  print(paste0("[",date(),"] ","running GO enrichment compare..."))
  ## GO 
  mego_MF <- compareCluster(geneSets,fun = "enrichGO",OrgDb = org.Osativa.eg.db, keyType= "GID",ont="MF",pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1)
  mego_BP <- compareCluster(geneSets,fun = "enrichGO",OrgDb = org.Osativa.eg.db, keyType= "GID",ont="BP",pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1)
  mego_CC <- compareCluster(geneSets,fun = "enrichGO",OrgDb = org.Osativa.eg.db, keyType= "GID",ont="CC",pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1)
  # mego_KEGG <- compareCluster(geneSets,fun = "enrichKEGG",organism = "dosa",keyType = "kegg",pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1)

  write.table(as.data.frame(mego_MF@compareClusterResult),file = paste0(name,"_compareCluster_GO_MF_",Sys.Date(),".xls"),quote = F,sep = "\t",col.names = T,row.names = F)
  write.table(as.data.frame(mego_BP@compareClusterResult),file = paste0(name,"_compareCluster_GO_BP_",Sys.Date(),".xls"),quote = F,sep = "\t",col.names = T,row.names = F)
  write.table(as.data.frame(mego_CC@compareClusterResult),file = paste0(name,"_compareCluster_GO_CC_",Sys.Date(),".xls"),quote = F,sep = "\t",col.names = T,row.names = F)
  # write.table(as.data.frame(mego_KEGG@compareClusterResult),file = paste0(name,"_compareCluster_KEGG_",Sys.Date(),".xls"),quote = F,sep = "\t",col.names = T,row.names = F)

  dotplot(mego_MF)+theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+labs(title = "GO:MF")
  ggsave(filename = paste0(name,"_GO_MF_",Sys.Date(),"_dotplot.pdf"),device = "pdf",width = 8,height = 10)
  dotplot(mego_BP)+theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+labs(title = "GO:BP")
  ggsave(filename = paste0(name,"_GO_BP_",Sys.Date(),"_dotplot.pdf"),device = "pdf",width = 8,height = 10)
  dotplot(mego_CC)+theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+labs(title = "GO:CC")
  ggsave(filename = paste0(name,"_GO_CC_",Sys.Date(),"_dotplot.pdf"),device = "pdf",width = 8,height = 10)
  # dotplot(mego_KEGG)+theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))+labs(title = "KEGG")
  # ggsave(filename = paste0(name,"_KEGG_",Sys.Date(),"_dotplot.pdf"),device = "pdf",width = 8,height = 10)


  for(i in 1:length(names(geneSets))){
    gs_MF <- mego_MF@compareClusterResult %>% filter(Cluster == names(geneSets)[i]) %>% filter(pvalue <= 0.05 & p.adjust <= 0.3)
    gs_BP <- mego_BP@compareClusterResult %>% filter(Cluster == names(geneSets)[i]) %>% filter(pvalue <= 0.05 & p.adjust <= 0.3)
    gs_CC <- mego_CC@compareClusterResult %>% filter(Cluster == names(geneSets)[i]) %>% filter(pvalue <= 0.05 & p.adjust <= 0.3)
    # gs_KEGG <- mego_KEGG@compareClusterResult %>% filter(Cluster == names(geneSets)[i]) %>% filter(pvalue <= 0.05 & p.adjust <= 0.3)
    
    aPEAR::enrichmentNetwork( gs_MF ,colorBy = 'p.adjust', colorType = 'pval',verbose = F,nodeSize = "Count") + labs(title = names(geneSets)[i]);
    ggsave(filename = paste0(name,"_GO_MF_",names(geneSets)[i],"_",Sys.Date(),"_top100Network.pdf"),device = "pdf",width = 6,height = 6)
    aPEAR::enrichmentNetwork( gs_BP ,colorBy = 'p.adjust', colorType = 'pval',verbose = F,nodeSize = "Count") + labs(title = names(geneSets)[i]);
    ggsave(filename = paste0(name,"_GO_BP_",names(geneSets)[i],"_",Sys.Date(),"_top100Network.pdf"),device = "pdf",width = 6,height = 6)
    aPEAR::enrichmentNetwork( gs_CC ,colorBy = 'p.adjust', colorType = 'pval',verbose = F,nodeSize = "Count") + labs(title = names(geneSets)[i]);
    ggsave(filename = paste0(name,"_GO_CC_",names(geneSets)[i],"_",Sys.Date(),"_top100Network.pdf"),device = "pdf",width = 6,height = 6)
    # aPEAR::enrichmentNetwork( gs_KEGG ,colorBy = 'p.adjust', colorType = 'pval',verbose = F,nodeSize = "Count") + labs(title = names(geneSets)[i]);
    # ggsave(filename = paste0(name,"_KEGG_",names(geneSets)[i],"_",Sys.Date(),"_top100Network.pdf"),device = "pdf",width = 6,height = 6)

    pathcl_BP <- aPEAR::findPathClusters(gs_BP,verbose = F)
    pathcl_MF <- aPEAR::findPathClusters(gs_MF,verbose = F)
    pathcl_CC <- aPEAR::findPathClusters(gs_CC,verbose = F)
    # pathcl_KEGG <- aPEAR::findPathClusters(gs_KEGG,verbose = F)
    write.table(pathcl_BP$clusters,file = paste0(name,"_GO_BP_",names(geneSets)[i],"_",Sys.Date(),"_top100Network_PathClusters.xls"),quote = F,sep = "\t",col.names = T,row.names = F)
    write.table(pathcl_MF$clusters,file = paste0(name,"_GO_MF_",names(geneSets)[i],"_",Sys.Date(),"_top100Network_PathClusters.xls"),quote = F,sep = "\t",col.names = T,row.names = F)
    write.table(pathcl_CC$clusters,file = paste0(name,"_GO_CC_",names(geneSets)[i],"_",Sys.Date(),"_top100Network_PathClusters.xls"),quote = F,sep = "\t",col.names = T,row.names = F)
    # write.table(pathcl_KEGG$clusters,file = paste0(name,"_KEGG_",names(geneSets)[i],"_",Sys.Date(),"_top100Network_PathClusters.xls"),quote = F,sep = "\t",col.names = T,row.names = F)    
  }
  print(paste0("[",date(),"] ","plot done!"))
}

### main 
genelist <- read.table(args[1])
name <- paste("6.DEGcompare/",basename(args[1]),sep="")
geneSets <- getGeneSetList(genelist)
list <- getIntersect(geneSets,name)
compare(list,name)
print(paste0("[",date(),"] ","All Finished!"))
