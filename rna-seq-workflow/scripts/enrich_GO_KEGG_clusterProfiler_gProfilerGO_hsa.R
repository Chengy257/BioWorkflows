#!/usr/bin/Rscript
## Usage:
## /usr/bin/Rscript /home/chengyu/workflows/enrich_GO_KEGG_clusterProfiler_gProfilerGO.R [gene list] [output dir prefix]     
.libPaths("/share/R/library/4.2.3/")
args <- commandArgs(T)
pkgs <- c('clusterProfiler','ggplot2','aPEAR','svglite','GOSemSim','enrichplot')
lapply(pkgs, function(x){
   suppressMessages(library(x, character.only = T))})

genelist <- read.table(args[1])
output_dir <- args[2]
out <- paste(output_dir,"/",sep="")
prefix <- gsub(".DEGs.txt","",basename(args[1]))
spe <- args[3]  ## species: osa or hsa 

if( spe == "osa"){
    genelist <- genelist[,1]
    if (!require("org.Osativa.eg.db", quietly = TRUE))
        install.packages("/share/data/reference/osa/GO/gProfiler/org.Osativa.eg.db_0.1.tar.gz",repos=NULL)
    suppressMessages(library("org.Osativa.eg.db"))
    db <- org.Osativa.eg.db
    keyType <- "GID"
} else if( spe == "hsa" ){
    suppressMessages(library("org.Hs.eg.db"))
	genelist$ENSEMBL <- unlist(strsplit(genelist[,1],split = "\\."))[seq(from=1,to=2*length(genelist[,1]),by=2)]
    ENTREZID <- bitr(genelist$ENSEMBL,fromType = 'ENSEMBL',toType = 'ENTREZID',OrgDb='org.Hs.eg.db')
    ENTREZID <- ENTREZID[!duplicated(ENTREZID$ENSEMBL), ]
	ENTREZID <- na.omit(ENTREZID)
	genelist <- ENTREZID$ENTREZID
	# suppressMessages(library("org.Hs.eg.db"))
    db <- org.Hs.eg.db  
    keyType <- "ENTREZID"
} else{
    print("Please chose species from hsa/osa for human or rice!")
    break
}

for( ont in c("MF","BP","CC","KEGG")){
    if(ont == "KEGG"){
        outname_prefix <- "KEGG"
        ego <- enrichKEGG(genelist, organism = spe ,keyType = "kegg",pvalueCutoff=1,pAdjustMethod="BH",qvalueCutoff=1)
        ego <- setReadable(ego, OrgDb = db, keyType = keyType)
        newnames <- c()
        for(i in 1:length(ego@result$Description)){
            name <- unlist(strsplit(ego@result$Description[i],split = " - Oryza sativa japonica")[[1]][1])
            newnames <- c(newnames,name)
        }
        ego@result$Description <- newnames
    }
    else{
        outname_prefix <- paste0("GO:",ont)
        ego <- enrichGO(genelist,OrgDb = db, ont= ont ,keyType = keyType, pvalueCutoff = 0.05,qvalueCutoff = 1,pAdjustMethod = "BH",readable = TRUE)
        
    }
    write.table(as.data.frame(ego@result),file = paste0(out,prefix,"_",outname_prefix,"_AllResults_",Sys.Date(),".xls"),quote = F,sep = "\t",col.names = T,row.names = F)
    write.table(as.data.frame(ego@result %>% filter(  pvalue <= 0.05)),file = paste0(out,prefix,"_",outname_prefix,"_Pvalue0.05_",Sys.Date(),".xls"),quote = F,sep = "\t",col.names = T,row.names = F)
    write.table(as.data.frame(ego@result %>% filter( p.adjust <= 0.05)),file = paste0(out,prefix,"_",outname_prefix,"_Padj0.05_",Sys.Date(),".xls"),quote = F,sep = "\t",col.names = T,row.names = F)
    tryCatch({ 
		dotplot(ego,showCategory=10)+labs(title = paste0(outname_prefix," top10 Terms"))
		ggsave(filename = paste0(out,prefix, "_",outname_prefix, "_",Sys.Date(),"_Top10_dotplot.pdf"),device = "pdf",width = 5,height = 5)
		barplot(ego,showCategory=10)+labs(title = paste0(outname_prefix," top10 Terms"))
		ggsave(filename = paste0(out,prefix, "_",outname_prefix, "_",Sys.Date(),"_Top10_barplot.pdf"),device = "pdf",width = 5,height = 5)
	},error = function(e) {
		print(paste0("Error while running plot!", outname_prefix, conditionMessage(e)))
	})

    
    tryCatch({
        edox <- pairwise_termsim(ego)
        treeplot(edox, cluster.params = list(method="ward.D", n = 5, label_words_n = 5))+labs(title=outname_prefix)
        ggsave(filename = paste0(out,prefix, "_",outname_prefix, "_",Sys.Date(),"_treeplot.pdf"),device = "pdf",width = 14,height = 7)
        emapplot(edox, showCategory = 10, layout="kk", cex_category=1.5,min_edge = 0.8)+labs(title=outname_prefix)
        ggsave(filename = paste0(out,prefix, "_",outname_prefix, "_",Sys.Date(),"_emapplot.pdf"),device = "pdf",width = 12,height = 12)
        cnetplot(ego,showCategory = 5, colorEdge = T, node_label="gene",circular = F)+labs(title = paste0(outname_prefix," top5 Terms"))
        ggsave(filename = paste0(out,prefix, "_",outname_prefix, "_",Sys.Date(),"_Top5_cnetplot.pdf"),device = "pdf",width = 12,height = 12)

    },error = function(e) {
        print(paste0("Error while running plot!", outname_prefix, conditionMessage(e)))
    })
    
    tryCatch({
    aPEAR::enrichmentNetwork(ego@result %>% head(100),colorBy = 'p.adjust', colorType = 'pval', drawEllipses = FALSE,repelLabels = TRUE,verbose = F,nodeSize = "Count")+labs(title = paste0(outname_prefix," top 100 Terms Network"))
    ggsave(filename = paste0(out,prefix,"_",outname_prefix,"_",Sys.Date(),"_top100Network.pdf"),device = "pdf",width = 6,height = 6)
    pathcl <- aPEAR::findPathClusters(ego@result %>% head(100),verbose = F)
    write.table(pathcl$clusters,file = paste0(out,prefix,"_",outname_prefix,"_",Sys.Date(),"_top100Network_PathClusters.xls"),quote = F,sep = "\t",col.names = T,row.names = F)
    },error = function(e) {
        print(paste0("Error while running aPEAR!", outname_prefix, conditionMessage(e)))
    })
}
# print(paste0("[",date(),"] ","plot done!"))
print(paste0("[",date(),"] ","All Finished!"))
