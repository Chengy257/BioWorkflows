#!/usr/bin/R
# .libPaths("/home/chengyu/R/Rlib_4.2.3")
.libPaths(c("/home/chengyu/R/Rlib_4.2.3","/opt/R/4.2.3/lib/R/library"))
# .libPaths()
# suppressMessages(library(ChIPseeker))
# suppressMessages(library(GenomicFeatures))
# getwd()

pkgs <- c('ChIPseeker','GenomicFeatures')
lapply(pkgs, function(x){
   suppressMessages(library(x, character.only = T))})

# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("pacman")
# suppressMessages(library(pacman))
# pacman::p_load('ChIPseeker','GenomicFeatures')

load_file <- function(peak_files) {  
  peaks <- list()  
  #peaks <- lapply(peak_files, function(x) readPeakFile(x))
  for ( i in 1:length(peak_files) ) { 
    filepath <- peak_files[i]
    print(filepath)
    sample <- gsub("_peaks.narrowPeak","", basename(peak_files[i]))
    if(file.info(filepath)$size > 0){
        peaks[sample] <- readPeakFile(filepath) 
    }
  }  	  
  return(peaks)
}

annoPeaks <- function(my_txdb, peaks) {
  system("mkdir -p 4.peak/anno_result/")
  peakAnnoList <- lapply(peaks, annotatePeak, TxDb=my_txdb, level="gene", verbose=FALSE, addFlankGeneInfo=TRUE, flankDistance=3000)
  outdir = "4.peak/anno_result"
  if (file.exists(outdir)) {
    cat("Output Dir existed!\n")
  } else {
    dir.create(outdir)
  }
  lapply(names(peakAnnoList), function(x) {
    file_name <- paste("4.peak/anno_result/", x, ".Anno.xls", sep = "")
    write.table(as.data.frame(peakAnnoList[[x]]), file_name, sep = "\t", row.names = FALSE, col.names=TRUE, quote = FALSE)
  })
  print("plot done!")
  return(peakAnnoList)
}

# main 
args <- commandArgs(T)
gtf_file <- args[1]
# gtf_file <- "/home/chengyu/references/osa/Oryza_sativa.IRGSP-1.0.56.Chr.gtf"
peak_files <- read.table(args[2])$V1

my_txdb <- makeTxDbFromGFF(gtf_file)
peaks <- load_file(peak_files)
peakAnnoList <- annoPeaks(my_txdb, peaks)

TSSRegion <- getPromoters(TxDb=my_txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(peaks, getTagMatrix, windows=TSSRegion)

pdf("4.peak/anno_result/Peakanno_PeakDistributions.pdf",height=8,width=8)
  plotAnnoBar(peakAnnoList)
  plotDistToTSS(peakAnnoList, title="Distribution of binding loci relative to TSS")
  plotAvgProf(tagMatrixList, xlim=c(-3000, 3000),conf=0.95,resample=500, facet="row")
  tagHeatmap(tagMatrixList)
dev.off()
