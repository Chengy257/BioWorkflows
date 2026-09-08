#!/usr/bin/env Rscript
# Differential expression for srna-seq (DESeq2), ported from the rna-seq
# workflow's run_deseq2.R. Reads one per-class count matrix
# (4.expression/{class}/{class}_counts.tsv, row.names=1 consumable) plus the
# srna-seq sample table (header sample_id[, group[, batch]]), fits
# design ~ batch + group (batch_correction "T") or ~ group, and writes, for
# every <treat> group against the control group, into the output directory:
#   {treat}_vs_{control}_DESeq2.output.tsv  (first header column: "feature")
#   {treat}_vs_{control}_DESeq2.output.tsv_VolcanoPlot.pdf / _MAPlot.pdf
#   {class}_DESeq2.normalized.vst.PCA_plot.pdf
#   {class}_DESeq2.normalized.vst.Pearson_heatmap.pdf
#   sessionInfo.txt
#
# Differences from the rna-seq original:
#   - the sample table is read by the srna-seq column names (sample_id/group/
#     batch) instead of id/group;
#   - the feature column is written natively as the first header column
#     instead of being prepended afterwards with sed;
#   - a -k/--class flag names the analyzed cascade class (plot titles/labels);
#   - count-matrix columns and sample ids are aligned explicitly (srna-seq
#     sample ids may contain "-"; reads use check.names=FALSE to keep them).
# R runtime and library paths are provided by config/software.yaml.
###############
# functions
preprocess_data <- function(count_data, sample_info, is_batch){
    ## colData rows must match the count-matrix columns; srna-seq sample ids
    ## may contain "-", so both sides are read/kept verbatim (check.names=F).
    if (!setequal(colnames(count_data), sample_info$sample_id)) {
        stop("Sample table sample_id values and count-matrix columns differ: '",
             paste(setdiff(sample_info$sample_id, colnames(count_data)), collapse = ", "),
             "' missing from the count matrix; '",
             paste(setdiff(colnames(count_data), sample_info$sample_id), collapse = ", "),
             "' missing from the sample table", call. = FALSE)
    }
    sample_info <- sample_info[match(colnames(count_data), sample_info$sample_id), , drop = FALSE]
    rownames(sample_info) <- sample_info$sample_id
    sample_info$group <- as.factor(sample_info$group)
    if (is_batch) {
        sample_info$batch <- as.factor(sample_info$batch)
        dds <- DESeqDataSetFromMatrix(countData = count_data,
                                      colData = sample_info,
                                      design = ~ batch + group)
    } else {
        sample_info$batch <- NULL
        dds <- DESeqDataSetFromMatrix(countData = count_data,
                                      colData = sample_info,
                                      design = ~ group)
    }
    dds
}
runDEseq2 <- function(dds, sample_info, output_name, FDR, FoldChange){
    dds <- DESeq(dds)
    all_groups <- levels(dds$group)
    ## enumerate all pairwise combinations; emit only <treat>_vs_<control>
    pairwise_combs <- expand.grid(all_groups, all_groups, stringsAsFactors = FALSE)
    pairwise_combs <- pairwise_combs[pairwise_combs$Var1 != pairwise_combs$Var2, ]
    print(pairwise_combs)
    norm_count <- round(DESeq2::counts(dds, normalized = TRUE), 3)
    n_contrasts <- 0
    for (i in seq_len(nrow(pairwise_combs))) {
        ref <- pairwise_combs[i, 2]    # pair[2]
        treat <- pairwise_combs[i, 1]  # pair[1]
        result_group <- paste0(treat, "_vs_", ref)
        ## Only output treatment_vs_control comparisons (exact match against
        ## the configured control group)
        if (ref == control_name) {
            print(paste("[", date(), "] DESeq2: preprocess ", result_group, "...", sep = ""))
            res_contrast <- results(dds, contrast = c("group", treat, ref))
            ids <- sample_info$sample_id[which(sample_info$group %in% c(ref, treat))]
            print(ids)
            out_norm_count <- norm_count[, ids, drop = FALSE]
            out <- as.data.frame(res_contrast)
            out$FoldChange <- 2^(out$log2FoldChange)
            ## first header column is the feature name, written natively
            out <- cbind(feature = rownames(out), out, out_norm_count)
            out$type <- "Unsig";
            out$type[which(out$padj <= FDR & out$FoldChange >= FoldChange)] <- "Up";
            out$type[which(out$padj <= FDR & out$FoldChange <= (1 / FoldChange))] <- "Down";
            out$type[is.na(out$type)] <- "Unsig";
            write.table(out, file = paste0(output_name, result_group, "_DESeq2.output.tsv"),
                        quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t");
            group_contrast_plot(res_contrast, out, paste0(output_name, result_group))
            n_contrasts <- n_contrasts + 1
        }
    }
    n_contrasts
}
sample_plots <- function(dds, output_name, ntop, klass){
    ## vst()/varianceStabilizingTransformation() stop on matrices where the
    ## default parametric dispersion fit fails ("all gene-wise dispersion
    ## estimates are within 2 orders of magnitude ...") — common for the
    ## small, narrow-range matrices of sRNA classes (a miRNA class can also
    ## hold fewer features than vst()'s nsub=1000 default, where vst()
    ## refuses directly). Try the default transform, then a local-fit
    ## dispersion fallback; if even that fails, skip the exploratory
    ## PCA/heatmap plots with a message — the contrast tables produced by
    ## DESeq() do not depend on them.
    vst_dds <- tryCatch(
        tryCatch({
            if (nrow(dds) >= 1000) {
                vst(dds, blind = FALSE)
            } else {
                varianceStabilizingTransformation(dds, blind = FALSE)
            }
        }, error = function(e) {
            message("vst/VST with the default dispersion fit failed (",
                    conditionMessage(e), "); retrying with a local fit")
            dds_loc <- estimateDispersionsGeneEst(dds, quiet = TRUE)
            dds_loc <- estimateDispersionsFit(dds_loc, fitType = "local", quiet = TRUE)
            varianceStabilizingTransformation(dds_loc, blind = FALSE)
        }),
        error = function(e) {
            message("Skipping the PCA/heatmap plots: the variance-stabilizing ",
                    "transform failed on this matrix (", conditionMessage(e), ")")
            NULL
        })
    if (is.null(vst_dds)) {
        return(invisible(NULL))
    }
    vstMat <- assay(vst_dds)
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
    # pearson correlation
    pearson_cor <- as.matrix(cor(vstMat, method = "pearson"))
    # cluster
    hc <- hcluster(t(vstMat), method = "pearson")
    # heatmap
    pdf(paste0(output_name, klass, "_DESeq2.normalized.vst.Pearson_heatmap.pdf"),
        height = 14, width = 12)
        heatmap.2(pearson_cor, Rowv = as.dendrogram(hc), symm = TRUE, trace = "none",
                  col = hmcol, margins = c(12, 12), main = "Samples' pearson correlation");
    dev.off()

    DESeq2::plotPCA(vst_dds, intgroup = c("group"), returnData = FALSE, ntop = ntop);
    ggsave(paste0(output_name, klass, "_DESeq2.normalized.vst.PCA_plot.pdf"),
           device = "pdf", height = 5, width = 5);
}
group_contrast_plot <- function(res, data, name){
    ## Volcano plot: dynamic coordinate ranges (no data points clipped),
    ## explicit color mapping, tolerance for empty categories, label
    ## positions scaled with the data (same shape as the rna-seq original).
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

    v_p <- ggplot(data, aes(x = log2FoldChange, y = -log10(padj), color = type)) +
        geom_point(alpha = 0.75, size = 1.2, na.rm = TRUE) +
        labs(title = name, x = expression(log[2](FoldChange)), y = expression(-log[10](adjusted_Pvalue))) +
        theme(plot.title = element_text(hjust = 0.4)) +
        geom_hline(yintercept = -log10(0.05), lty = 4, lwd = 0.6, alpha = 0.8) +
        geom_vline(xintercept = c(1, -1), lty = 4, lwd = 0.6, alpha = 0.8) +
        theme_bw() +
        scale_color_manual(values = c(Up = "red", Down = "blue", Unsig = "grey"),
                           breaks = c("Down", "Unsig", "Up"), na.translate = FALSE) +
        coord_cartesian(xlim = c(-xmax, xmax), ylim = c(0, ymax)) +
        annotate("text", x = -xmax * 0.8, y = ymax * 0.95,
                 label = paste0("Down: ", freq_of("Down")), color = "black", family = "mono", fontface = "plain") +
        annotate("text", x = xmax * 0.8, y = ymax * 0.95,
                 label = paste0("Up: ", freq_of("Up")), color = "black", family = "mono", fontface = "plain")
    ggsave(v_p, filename = paste0(name, "_VolcanoPlot.pdf"), device = "pdf", width = 5, height = 4)
    pdf(paste0(name, "_MAPlot.pdf"))
        plotMA(res, main = name)
    dev.off()
}
warn_small_groups <- function(sample_info){
    ## A single replicate per group is legal for DESeq2 but statistically
    ## weak: warn, do not stop (the parse-time validator warns the same way).
    tab <- table(sample_info$group)
    for (lv in names(tab)) {
        if (tab[[lv]] < 2) {
            warning(sprintf("Group '%s' has only %d replicate(s); DESeq2 will run but the statistics are weak",
                            lv, tab[[lv]]), call. = FALSE)
        }
    }
}
# main function
main <- function(count, sample, is_batch, Nthreads, output_name, FDR, FoldChange, Ntop, klass){
    register(MulticoreParam(as.numeric(Nthreads)));
    print(paste("[", date(), "] Using ", Nthreads, " threads.", sep = ""))
    # output directory (created when missing; no extra subdirectory)
    if (!dir.exists(output_name)) {
        dir.create(output_name, recursive = TRUE)
    }
    print(paste("[", date(), "] Preprocess data (", klass, ")...", sep = ""))
    count_data <- read.delim(count, header = TRUE, row.names = 1, sep = "\t",
                             check.names = FALSE)  # read the class count matrix
    count_data <- count_data[rowSums(count_data) > 0, , drop = FALSE]  ## drop un-expressed features
    sample_info <- read.csv(sample, header = TRUE, sep = ",",
                            colClasses = "character", check.names = FALSE)  ## srna-seq sample table
    required <- c("sample_id", "group")
    if (!all(required %in% colnames(sample_info))) {
        stop("Sample table must carry the columns sample_id and group (header: ",
             "sample_id[, group[, batch]]); got: ",
             paste(colnames(sample_info), collapse = ", "), call. = FALSE)
    }
    if (is_batch && !("batch" %in% colnames(sample_info))) {
        stop("batch_correction is T but the sample table has no 'batch' column", call. = FALSE)
    }
    dds <- preprocess_data(count_data, sample_info, is_batch)
    print(dds)
    warn_small_groups(sample_info)
    #
    print(paste("[", date(), "] Running sample cluster & plot...", sep = ""))
    sample_plots(dds, output_name, Ntop, klass)
    #
    print(paste("[", date(), "] Running DESeq2...", sep = ""))
    runDEseq2(dds, sample_info, output_name, FDR, FoldChange)
    ## Record the R sessionInfo
    writeLines(capture.output(sessionInfo()), paste0(output_name, "sessionInfo.txt"))
    print(paste("[", date(), "] All done!", sep = ""))
}

## call main function
pkgs <- c('DESeq2', 'ggplot2', 'BiocParallel', 'gplots', 'RColorBrewer', 'amap', 'getopt')
lapply(pkgs, function(x){
   suppressMessages(library(x, character.only = TRUE))})
spec <- matrix(c("count", "c", 2, "character", "Input per-class count matrix, [filename, tab-separated file].",
                 "sample", "s", 2, "character", "Input sample table (header sample_id[,group[,batch]]), [filename, comma separated csv file].",
                 "output", "o", 2, "character", "Output directory, [filename].",
                 "batch", "b", 1, "logical", "If considering the batch effect (design ~ batch + group), True or False, [optional, default is False].",
                 "control_name", "r", 1, "character", "The control group name in the sample table, [optional, default is control].",
                 "fdr", "p", 1, "numeric", "DEGs adjusted pvalue threshold, [optional, default is 0.05].",
                 "foldchange", "f", 1, "numeric", "DEGs foldchange threshold, [optional, default is 2].",
                 "ntop", "n", 1, "numeric", "Number of top variable features used in PCA, [optional, default is 2000].",
                 "class", "k", 1, "character", "Cascade class analyzed, used in plot names, [optional, default is miRNA].",
                 "threads", "t", 1, "numeric", "Using CPU numbers, [optional, default is 1].",
                 "help", "h", 0, "logical", "Show this help information."),
                 byrow = TRUE, ncol = 5)
opt <- getopt(spec = spec)
print(opt)
# check
if (!is.null(opt$help) || is.null(opt$count) || is.null(opt$sample) || is.null(opt$output)){
    cat(paste(getopt(spec = spec, usage = TRUE), "\n"))
    quit()
}
if (is.null(opt$threads)) {
    opt$threads <- as.numeric("1") }
if (is.null(opt$batch)) {
    opt$batch <- as.logical("False") }
if (is.null(opt$fdr)) {
    opt$fdr <- as.numeric("0.05") }
if (is.null(opt$foldchange)) {
    opt$foldchange <- as.numeric("2") }
if (is.null(opt$ntop)) {
    opt$ntop <- as.numeric("2000") }
if (is.null(opt$class)) {
    opt$class <- as.character("miRNA")}
if (is.null(opt$control_name)) {
    opt$control_name <- as.character("control")}
#
count <- as.character(opt$count)
sample <- as.character(opt$sample)
is_batch <- as.logical(opt$batch)
Nthreads <- as.numeric(opt$threads)
FDR <- as.numeric(opt$fdr)
FoldChange <- as.numeric(opt$foldchange)
Ntop <- as.numeric(opt$ntop)
klass <- as.character(opt$class)
control_name <- as.character(opt$control_name)
output_name <- trimws(as.character(opt$output), which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
if (!endsWith(output_name, "/")) {
    output_name <- paste0(output_name, "/")
}

# running main function
main(count, sample, is_batch, Nthreads, output_name, FDR, FoldChange, Ntop, klass)
