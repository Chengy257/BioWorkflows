#!/usr/bin/Rscript
## Usage:
## /usr/bin/Rscript /home/chengyu/workflows/enrich_GO_KEGG_clusterProfiler_gProfilerGO.R [gene list] [output dir prefix]     
.libPaths("/home/chengyu/R/Rlib_4.2.3")
.libPaths()
args <- commandArgs(T)
pkgs <- c('clusterProfiler','ggplot2','aPEAR','svglite','magrittr')
lapply(pkgs, function(x){
   suppressMessages(library(x, character.only = T))})

if (!require("org.Osativa.eg.db", quietly = TRUE))
    install.packages("/home/chengyu/references/osa/GO/gProfiler/org.Osativa.eg.db_0.1.tar.gz",repos=NULL)
suppressMessages(library(org.Osativa.eg.db))
genelist <- read.table(args[1])[,1]
head(genelist)
output_dir <- args[2]
out <- paste(output_dir,"/",sep="")
# prefix <- gsub(".DEGs.txt","",basename(args[1]))
prefix <- basename(args[1])
###############################################
## KEGG
R.utils::setOption("clusterProfiler.download.method",'auto')
#T2G <- read.table("/share/data/reference/osa/GO/KEGG/KEGG_TranscriptID2GeneIDs",header=F,sep="\t")
## translist <- T2G[(which( genelist %in% T2G[,2])),1] 
#tmp <- merge(x=T2G,y=genelist,by.x=2,by.y=1)
#translist <- unique(tmp[,2])
# head(translist)
###############################################
## GO
for( ont in c("MF","BP","CC","KEGG")){
    if(ont == "KEGG"){
		#ego <- enrichKEGG(translist, organism = "dosa",keyType = "kegg",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1)
        ego <- enrichKEGG(genelist, organism = "dosa",keyType = "kegg",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1)
        outname_prefix <- "KEGG"
        newnames <- c()
        for(i in 1:length(ego@result$Description)){
            name <- unlist(strsplit(ego@result$Description[i],split = " - Oryza sativa japonica")[[1]][1])
            newnames <- c(newnames,name)
        }
        ego@result$Description <- newnames
    }
    else{
        ego <- enrichGO(genelist,OrgDb = org.Osativa.eg.db,ont=ont,keyType = "GID",pvalueCutoff = 1,qvalueCutoff = 1,pAdjustMethod = "BH")
        outname_prefix <- paste0("GO:",ont)
    }
    write.table(as.data.frame(ego@result),file = paste0(out,prefix,"_",outname_prefix,"_AllResults_",Sys.Date(),".xls"),quote = F,sep = "\t",col.names = T,row.names = F)
	# write.table(as.data.frame(ego@result %>% filter(  pvalue <= 0.05)),file = paste0(out,prefix,"_",outname_prefix,"_PvalueSignificant_",Sys.Date(),".xls"),quote = F,sep = "\t",col.names = T,row.names = F)
    write.table(as.data.frame(ego@result %>% filter( p.adjust <= 0.05)),file = paste0(out,prefix,"_",outname_prefix,"_PadjSignificant_",Sys.Date(),".xls"),quote = F,sep = "\t",col.names = T,row.names = F)
    
    dotplot(ego,showCategory=10)+labs(title = paste0(outname_prefix," top10 Terms"))
    ggsave(filename = paste0(out,prefix, "_",outname_prefix, "_",Sys.Date(),"_dotplot.pdf"),device = "pdf",width = 5,height = 5)
    
    tryCatch({
#aPEAR::enrichmentNetwork(ego@result %>% head(100),colorBy = 'p.adjust', colorType = 'pval', drawEllipses = FALSE,repelLabels = TRUE,verbose = F,nodeSize = "Count")+labs(title = paste0(outname_prefix,"top 100 Terms Network"))
    aPEAR::enrichmentNetwork(ego@result %>% filter(  pvalue <= 0.05),colorBy = 'p.adjust', colorType = 'pval', drawEllipses = FALSE,repelLabels = TRUE,verbose = F,nodeSize = "Count")+labs(title = paste0(outname_prefix,"top 100 Terms Network"))
    ggsave(filename = paste0(out,prefix,"_",outname_prefix,"_",Sys.Date(),"_sigNetwork.pdf"),device = "pdf",width = 6,height = 6)
    
	#pathcl <- aPEAR::findPathClusters(ego@result %>% head(100),verbose = F)
    pathcl <- aPEAR::findPathClusters(ego@result %>% filter(  pvalue <= 0.05),verbose = F)
    write.table(pathcl$clusters,file = paste0(out,prefix,"_",outname_prefix,"_",Sys.Date(),"_sigNetwork_PathClusters.xls"),quote = F,sep = "\t",col.names = T,row.names = F)
    
    },error = function(e) {
        print(paste0("Error while running aPEAR!", outname_prefix, conditionMessage(e)))
    })
}
# print(paste0("[",date(),"] ","plot done!"))
print(paste0("[",date(),"] ","All Finished!"))
