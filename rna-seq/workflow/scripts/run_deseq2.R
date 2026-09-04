#!/usr/bin/env Rscript
# Usage:
#       Rscript run_deseq2.R -c count.matrix.tsv -s sample_info.csv -o 5.DEG ...
# R runtime and library paths are provided by config/software.yaml
###############
# functions
preprocess_data <- function(count_data,sample_info,is_batch){
    if(is_batch){   #
        colnames(sample_info) <-  c("id","group","batch")
        sample_info$group <- as.factor(sample_info$group)
        sample_info$batch <- as.factor(sample_info$batch)
        dds <- DESeqDataSetFromMatrix(countData = count_data,
                                colData = sample_info,
                                design = ~ batch + group)
    }
    else{
        sample_info <- sample_info[,1:2]
        colnames(sample_info) <-  c("id","group")
        sample_info$group <- as.factor(sample_info$group)
        dds <- DESeqDataSetFromMatrix(countData = count_data,
                                colData = sample_info,
                                design = ~ group)
    }
    dds
}
runDEseq2 <- function(dds,sample_info,output_name,FDR,FoldChange){
    dds<- DESeq(dds)
    all_groups <- levels(dds$group)
    ## get all pairwise combinations
    pairwise_combs <- expand.grid(all_groups, all_groups, stringsAsFactors = FALSE)
    pairwise_combs <- pairwise_combs[pairwise_combs$Var1 != pairwise_combs$Var2,]
    print(pairwise_combs)
    norm_count <- round(DESeq2::counts(dds,normalized=T),3)
    for( i in 1:length(pairwise_combs[,1])){
        ref <- pairwise_combs[i,2]    # pair[2]
        treat <- pairwise_combs[i,1]  # pair[1]
        ids <- c()
        result_group <- paste0(treat, "_vs_", ref)
        ## Only output treatment_vs_control comparisons (exact match; fixes the old grepl substring mismatch)
        if(ref == control_name) {
            print(paste("[",date(),"] DESeq2: preprocess ",result_group,"...",sep=""))
            res_contrast <- results(dds, contrast = c("group",treat, ref))
            ids <- sample_info$id[which(sample_info$group %in% c(ref,treat))]
            print(ids)
            out_norm_count <- norm_count[,ids]
            out <- as.data.frame(res_contrast);
            out$FoldChange <- 2^(out$log2FoldChange);
            out <- cbind(out,out_norm_count);
            out$type <- NA ;
            out$type[which(out$padj <= FDR & out$FoldChange >= FoldChange)] <- "Up";
            out$type[which(out$padj <= FDR & out$FoldChange <= (1/FoldChange))] <- "Down";
            out$type[is.na(out$type)] <- "Unsig";
            head(out)
            write.table(out,file=paste(output_name,result_group,"_DESeq2.output.tsv",sep=""),quote = F,row.names = T,col.names = T,sep = "\t");
            system(paste("sed -i '1 s/^/gene_id\t/'",paste(output_name,result_group,"_DESeq2.output.tsv",sep="")))
            group_contrast_plot(res_contrast,out,paste(output_name,result_group,sep=""))
        }
    }
}
sample_plots <- function(dds,output_name,ntop){
    vst_dds <- vst(dds, blind=FALSE)
    vstMat <- assay(vst_dds)
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
    # pearson correlation
    pearson_cor <- as.matrix(cor(vstMat, method="pearson"))
    # cluster
    hc <- hcluster(t(vstMat), method="pearson")
    # heatmap
    pdf(paste(output_name,"DESeq2.normalized.vst.Pearson_heatmap.pdf",sep=""),height=14,width=12)
        heatmap.2(pearson_cor, Rowv=as.dendrogram(hc), symm=T, trace="none",col=hmcol, margins=c(12,12), main="Samples' pearson correlation");
    dev.off()

    DESeq2::plotPCA(vst_dds, intgroup=c("group"), returnData=FALSE, ntop=ntop);
    ggsave(paste(output_name,"DESeq2.normalized.vst.PCA_plot.pdf",sep=""),device="pdf",height=5,width=5);

}

group_contrast_plot <- function(res,data,name){
    ## Volcano plot (fix P1-8): dynamic coordinate ranges (no data points clipped), fontface,
    ## explicit color mapping, tolerance for empty categories, label positions scaled with the data
    anno <- as.data.frame(table(data$type))
    freq_of <- function(lv){
        v <- anno$Freq[as.character(anno$Var1) == lv]
        if (length(v) == 0) 0L else v
    }
    ys <- suppressWarnings(-log10(data$padj))
    ys <- ys[is.finite(ys)]
    ymax <- max(c(2, ys), na.rm = TRUE) * 1.05
    xs <- abs(data$log2FoldChange)
    xs <- xs[is.finite(xs)]
    xmax <- max(c(2, xs), na.rm = TRUE) * 1.05

    v_p <- ggplot(data,aes(x=log2FoldChange,y =-log10(padj),color=type))+
        geom_point(alpha=0.75, size=1.2, na.rm=TRUE)+
        labs(title=name, x=expression(log[2](FoldChange)),y=expression(-log[10](adjusted_Pvalue)))+
        theme(plot.title = element_text(hjust = 0.4))+
        geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
        geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+
        theme_bw()+
        scale_color_manual(values = c(Up="red", Down="blue", Unsig="grey"),
                           breaks = c("Down","Unsig","Up"), na.translate = FALSE)+
        coord_cartesian(xlim=c(-xmax,xmax), ylim=c(0,ymax))+
        annotate("text", x=-xmax*0.8, y=ymax*0.95,
                 label=paste0("Down: ", freq_of("Down")), color="black", family="mono", fontface="plain")+
        annotate("text", x=xmax*0.8, y=ymax*0.95,
                 label=paste0("Up: ", freq_of("Up")), color="black", family="mono", fontface="plain")
    ggsave(v_p,filename=paste(name,"_VolcanoPlot.pdf",sep=""),device="pdf",width=5,height=4)
    pdf(paste(name,"_MAPlot.pdf",sep=""))
        plotMA(res,main=name)
    dev.off()
}
# main function
main <- function(count,sample,is_batch,Nthreads,output_name,FDR,FoldChange,Ntop){
    #
    register(MulticoreParam(as.numeric(Nthreads)));
    print(paste("[",date(),"] Using ",Nthreads," threads.",sep=""))
    # dir
    output_name <- paste(output_name,"/Diff_Expr_Analysis_Results/",sep="")
    if(dir.exists(output_name)){
        print("dir existed!")
    } else {
        system(paste0("mkdir -p ",output_name))
    }
    print(paste("[",date(),"] Preprocess data...",sep=""))
    count_data <- read.csv(count,header = T,row.names = 1,sep="\t") # read in raw count data
    count_data <- count_data[rowSums(count_data)>0,]   ## remove un-expressed genes
    sample_info <- read.csv(sample,header = T,sep=",") ## read in sample info
    ## Backward compatibility with historical data: '-' in ids replaced with '_' (consistent with count matrix column names;
    ## the validate_samples.py validator already forbids '-' in new projects, so this is redundant protection)
    sample_info[,1] <- gsub("-","_",sample_info[,1])
    dds <- preprocess_data(count_data,sample_info,is_batch)
    print(dds)
    #
    print(paste("[",date(),"] Running sample cluster & plot...",sep=""))
    sample_plots(dds,output_name,Ntop)
    #
    print(paste("[",date(),"] Running DESeq2...",sep=""))
    runDEseq2(dds,sample_info,output_name,FDR,FoldChange)
    ## Record the R sessionInfo (P2-9)
    writeLines(capture.output(sessionInfo()), paste0(output_name,"sessionInfo.txt"))
    print(paste("[",date(),"] All done!",sep=""))
}

## call main function
pkgs <- c('DESeq2','ggplot2','BiocParallel','gplots','RColorBrewer','amap','getopt')
lapply(pkgs, function(x){
   suppressMessages(library(x, character.only = T))})
spec <- matrix(c("count","c",2,"character","Input raw count matrix, [filename, tab-separated file].",
                 "sample","s",2,"character","Input sample infomation matrix, [filename, comma separated csv file].",
                 "output","o",2,"character","Output filename prefix.",
                 "batch","b",1,"logical","If considering the batch effect, True or False, [optional, default is False].",
                 "control_name","r",1,"character","The control group name in sample csv, [optional, default is control].",
                 "fdr","p",1,"numeric","DEGs adjusted pvalue threshold, [optional, default is 0.05].",
                 "foldchange","f",1,"numeric","DEGs foldchange threshold, [optional, default is 2].",
                 "ntop","n",1,"numeric","Number of top variable genes used in PCA, [optional, default is 20000].",
                 "threads","t",1,"numeric","Using CPU numbers, [optional, default is 1].",
                 "help", "h", 0, "logical", "Show this help information."),
                 byrow=T,ncol=5)
opt <- getopt(spec = spec)
print(opt)
# check
if( !is.null(opt$help) || is.null(opt$count) || is.null(opt$sample)|| is.null(opt$output)){
    cat(paste(getopt(spec=spec, usage = T), "\n"))
    quit()
}
if ( is.null(opt$threads) ) {
    opt$threads <- as.numeric("1") }
if ( is.null(opt$batch) ) {
    opt$batch <- as.logical("False") }
if ( is.null(opt$fdr) ) {
    opt$fdr <- as.numeric("0.05") }
if ( is.null(opt$foldchange) ) {
    opt$foldchange <- as.numeric("2") }
if ( is.null(opt$ntop) ) {
    opt$ntop <- as.numeric("20000") }
if ( is.null(opt$control_name) ) {
    opt$control_name <- as.character("control")}
#
count <- as.character(opt$count)
sample <- as.character(opt$sample)
is_batch <- as.logical(opt$batch)
Nthreads <- as.numeric(opt$threads)
FDR <- as.numeric(opt$fdr)
FoldChange <- as.numeric(opt$foldchange)
Ntop <- as.numeric(opt$ntop)
control_name <- as.character(opt$control_name)
output_name <- trimws(as.character(opt$output), which = c("both", "left", "right"), whitespace = "[ \t\r\n]")

# running main function
main(count,sample,is_batch,Nthreads,output_name,FDR,FoldChange,Ntop)
