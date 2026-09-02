#!/usr/bin/R
# .libPaths("/home/chengyu/R/Rlib_4.2.3")
.libPaths(c("/home/chengyu/R/Rlib_4.2.3","/opt/R/4.2.3/lib/R/library"))
pkgs <- c('ChIPseeker','GenomicFeatures')
lapply(pkgs, function(x){
   suppressMessages(library(x, character.only = T))})
args <- commandArgs(T)
gtf_file <- args[1]
# gtf_file <- "/home/chengyu/references/osa/Oryza_sativa.IRGSP-1.0.56.Chr.gtf"
peak_files <- args[2]
my_txdb <- makeTxDbFromGFF(gtf_file)
peaks <- readPeakFile(peak_files)
peakAnno <- annotatePeak(peaks,TxDb=my_txdb, level="gene", verbose=FALSE, addFlankGeneInfo=TRUE, flankDistance=3000)
write.table(as.data.frame(peakAnno), paste0(peak_files,"_peakAnnoTable.xls"), sep = "\t", row.names = FALSE, col.names=TRUE, quote = FALSE)
TSSRegion <- getPromoters(TxDb=my_txdb, upstream=3000, downstream=3000)
tagMatrix <-  getTagMatrix(peaks,windows=TSSRegion)  # lapply(peaks, getTagMatrix, windows=TSSRegion)
pdf(paste0(peak_files,"_peakAnnoPlot.pdf"),height=4,width=6)
  plotAnnoBar(peakAnno)
  plotAnnoPie(peakAnno)
  vennpie(peakAnno)
  plotDistToTSS(peakAnno, title="Distribution of peak relative to TSS")
  plotAvgProf(tagMatrix, xlim=c(-3000, 3000),conf=0.95,resample=500, facet="row")
  # tagHeatmap(tagMatrix,xlim=c(-3000, 3000))
dev.off()