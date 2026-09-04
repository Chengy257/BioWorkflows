#!/usr/bin/env Rscript
## featureCounts 定量（Rsubread）+ FPKM/TPM 计算，输出 5 列表格
## 依赖环境见 envs/quant.yaml；不再依赖个人 R 库路径
suppressWarnings({
pkgs <- c('argparser','Rsubread','limma','edgeR','getopt')
lapply(pkgs, function(x){
   suppressMessages(library(x, character.only = T))})
})

spec <- matrix(c("bam","b",2,"character","Input bam file.",
				"gtf","g",2,"character","Input gtf file.",
				"strand","s",1,"numeric","Strandedness info: 0 (unstranded), 1 (stranded) and 2 (reversely stranded).",
				"isPairedEnd","i",1,"logical","Is library layout PairedEnd, True/False, default is True.",
				"attribute","a",1,"character","Specify attribute type in GTF annotation, gene_id/transcript_id, default is gene_id.",
				"featureType","f",1,"character","Specify feature type in GTF annotation, exon/CDS, default is exon.",
				"output","o",2,"character","Output filname prefix.",
				"threads","t",1,"numeric","Using CPU numbers, optional, default is 1.",
				"help", "h", 0, "logical", "Show this help information."),
				byrow=T,ncol=5)
opt <- getopt(spec = spec)
# print(opt)

if( !is.null(opt$help) || is.null(opt$bam) || is.null(opt$gtf)|| is.null(opt$output)){
	cat(paste(getopt(spec=spec, usage = T), "\n"))
	quit()
}
if ( is.null(opt$threads) ) {
	opt$threads <- as.numeric("1") }
if ( is.null(opt$strand) ) {
	opt$strand <- as.numeric("0") }
if ( is.null(opt$isPairedEnd) ) {
	opt$isPairedEnd <- as.logical("True") }
if ( is.null(opt$attribute) ) {
	opt$attribute <- as.character("gene_id") }
if ( is.null(opt$featureType) ) {
	opt$featureType <- as.character("exon") }
bamFile<- opt$bam
gtfFile<- opt$gtf
nthreads<- opt$threads
outFilePref<- opt$output
strand<- opt$strand  # 0 (unstranded), 1 (stranded) and 2 (reversely stranded)
layout<- opt$isPairedEnd
attribute <- opt$attribute
featureType <- opt$featureType

#######################
outStatsFilePath<- paste(outFilePref, '.log',sep='');
outCountsFilePath<- paste(outFilePref,'.count',sep='');
fCountsList=featureCounts(bamFile,annot.ext=gtfFile,isGTFAnnotationFile=TRUE,GTF.attrType=attribute,GTF.featureType=featureType,nthreads=nthreads,isPairedEnd=layout,strandSpecific=strand)
dgeList=DGEList(counts=fCountsList$counts,genes=fCountsList$annotation)
fpkm=rpkm(dgeList,dgeList$genes$Length)
tpm=exp(log(fpkm)-log(sum(fpkm))+log(1e6))
write.table(fCountsList$stat,outStatsFilePath,sep="\t",col.names = FALSE,row.names = FALSE,quote=FALSE)
featureCounts=cbind(fCountsList$annotation[,1],dgeList$genes$Length,fCountsList$counts,fpkm,tpm)
colnames(featureCounts)=c('id','effLength','counts','fpkm','tpm')
write.table(featureCounts,outCountsFilePath,sep="\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
