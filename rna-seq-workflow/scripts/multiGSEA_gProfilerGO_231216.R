#!/usr/bin/Rscript
args <- commandArgs(T)
.libPaths("/home/chengyu/R/Rlib_4.2.3")
.libPaths()
pkgs <- c('clusterProfiler','ggplot2','circlize','GseaVis','enrichplot','dplyr')
lapply(pkgs, function(x){
   suppressMessages(library(x, character.only = T))})

if (!require("org.Osativa.eg.db", quietly = TRUE))
    install.packages("/home/chengyu/references/osa/GO/gProfiler/org.Osativa.eg.db_0.1.tar.gz",repos=NULL)
suppressMessages(library(org.Osativa.eg.db))

filenamelist <- read.table(args[1])
print(filenamelist)
# mygmt <- read.gmt(args[2])
out_prefix <- as.character(args[2])

spe <- args[3]  ## species: osa or hsa 


if(dir.exists(out_prefix) ==  FALSE){
  system(paste0("mkdir -p ",out_prefix))
}






# mylist = list()
# listname = c()
ONT=c("BP","MF","CC")
for(i in 1:length(filenamelist[,1])){
  filename <- filenamelist[i,1]
  print(paste0("parsing ",filename," ..."))
  df <- read.table(filename,header=T) 
  groupname <- strsplit(basename(filename),split = "_FoldChange")[[1]][1]
  # genelist <- df
  # if( spe == "osa"){
  #   genelist <- genelist[,1]
  #   if (!require("org.Osativa.eg.db", quietly = TRUE))
  #       install.packages("/share/data/reference/osa/GO/gProfiler/org.Osativa.eg.db_0.1.tar.gz",repos=NULL)
  #   suppressMessages(library("org.Osativa.eg.db"))
  #   db <- org.Osativa.eg.db
  #   keyType <- "GID"
  # } else if( spe == "hsa" ){
  #   genelist$ENSEMBL <- unlist(strsplit(genelist[,1],split = "\\."))[seq(from=1,to=2*length(genelist[,1]),by=2)]
  #   ENTREZID <- bitr(genelist$ENSEMBL,fromType = 'ENSEMBL',toType = 'ENTREZID',OrgDb='org.Hs.eg.db')
  #   genelist <- ENTREZID$ENTREZID
  #   suppressMessages(library("org.Hs.eg.db"))
  #   db <- org.Hs.eg.db  
  #   keyType <- "ENTREZID"
  # } else{
  #   print("Please chose species from hsa/osa for human or rice!")
  #   break
  # }

  
  fc <- as.vector(df$log2FoldChange[order(df$log2FoldChange,decreasing = T)])  
  names(fc) <- df$gene_id[order(df$log2FoldChange,decreasing = T)]
  for(j in 1:3){
    tryCatch({

      # if(){
        res <- clusterProfiler::gseGO(geneList = fc ,OrgDb = org.Osativa.eg.db, ont = ONT[j],keyType = "GID",pvalueCutoff = 1)  
      # }else {
      #   res <- clusterProfiler::gseGO(geneList = fc ,OrgDb = org.Osativa.eg.db, ont = ONT[j],keyType = "GID",pvalueCutoff = 1)
      # }
      
      outname <- paste0(out_prefix,"/",groupname,"_",ONT[j],"_multiGSEA_",Sys.Date()) 
      write.table(as.data.frame(res@result),file=paste0(outname,"_ResultTable_all.xls"),col.name=T,row.name=F,sep="\t",quote=F)
      write.table(as.data.frame(res@result %>% filter( abs(NES)>=1 & pvalue <= 0.05 & qvalue<=0.25)),file=paste0(outname,"_ResultTable_significant.xls"),col.name=T,row.name=F,sep="\t",quote=F)
      # gseaNb(res,geneSetID = res@result$ID[1:10],curveCol = circlize::rand_color(10))
      # ggsave(file=paste0(outname,"_multiGSEA_lineplot_",Sys.Date(),".pdf"),device="pdf")
      dotplot(res,showCategory=20)+facet_grid(~.sign)+labs(title=paste0(groupname," GO:",ONT[j]))
      ggsave(file=paste0(outname,"_multiGSEA_dotplot_",Sys.Date(),".pdf"),device="pdf")
      print(paste0("Finished ",groupname," GO:",ONT[j],"..."))
    },error = function(e) {
      print(paste0("Error processing ", groupname, " GO:", ONT[j], ": ", conditionMessage(e)))
    })
  }
}

print("All done!")
