#!/usr/bin/env Rscript
# Usage:
#       Rscript run_dmr.R -f s1.CpG_merged.CpG_report.merged_CpG_evidence.cov.gz,... \
#           -s samples.csv -o 6.DMR -r control -q 0.01 -d 25 ...
#
# Differential methylation (methylKit) over the merged CpG tables produced by
# coverage2cytosine --merge_CpG. The tables carry the bismark 6-column
# coverage layout (chromosome, start, end, methylation percentage, count
# methylated, count unmethylated), read natively via pipeline="bismarkCoverage"
# (methylKit recomputes the percentage from the counts; gzipped input is fine).
#
# Per non-control sample-table group one pairwise contrast (treatment vector
# control=0 / treat=1) produces, under the output directory:
#   {treat}_vs_{control}_DMC_all.tsv / _DMC_hyper.tsv / _DMC_hypo.tsv
#       per-CpG differential sites: chr, start, end, strand, meth.diff
#       (percent, treatment minus control), pvalue, qvalue (SLIM)
#   {treat}_vs_{control}_DMR_tiles.tsv
#       the same columns for tiled regions (win.size/step.size/cov.bases)
#   DMR_summary.tsv   one row per contrast with the table counts
#   sessionInfo.txt   R session record
# Nothing is written outside the output directory.
# R runtime and library paths are provided by config/software.yaml.
###############
# functions
read_design <- function(sampleinfo) {
    ## Sample table (validated at workflow parse time): sample_id[,group[,batch]],
    ## read as plain strings so group names keep their exact spelling.
    design <- read.csv(sampleinfo, header = TRUE, stringsAsFactors = FALSE,
                       colClasses = "character")
    if (!"sample_id" %in% names(design)) {
        stop("sample table ", sampleinfo, " has no sample_id column")
    }
    design$sample_id <- trimws(design$sample_id)
    if (anyDuplicated(design$sample_id) > 0) {
        stop("sample table ", sampleinfo, " contains duplicate sample_id values")
    }
    design
}

calc_diff <- function(meth_base, batch_values, is_batch, Nthreads) {
    ## Differential methylation on one united methylBase; the optional batch
    ## covariate is a one-column data.frame whose rows follow the sample order
    ## of the base object (methylKit requires covariates+2 < sample count,
    ## guaranteed by the workflow's >= 2 replicates per group design check).
    if (is_batch) {
        calculateDiffMeth(meth_base,
                          covariates = data.frame(batch = batch_values),
                          mc.cores = Nthreads)
    } else {
        calculateDiffMeth(meth_base, mc.cores = Nthreads)
    }
}

write_dm <- function(obj, type, file, difference, qvalue) {
    ## Write one DMC/DMR table (see the column contract in the header) and
    ## return the row count for DMR_summary.tsv. getData() unwraps the
    ## methylDiff object into a plain data.frame: direct column subsetting on
    ## the S4 object is rejected by methylKit.
    dm <- getData(getMethylDiff(obj, difference = difference,
                                qvalue = qvalue, type = type))
    out <- dm[, c("chr", "start", "end", "strand", "meth.diff", "pvalue", "qvalue")]
    write.table(out, file = file, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = TRUE)
    nrow(out)
}

run_contrast <- function(files, ids, design, control, treat, is_batch, Nthreads,
                         output_name, qvalue, mindiff, tilelen, tilestep,
                         mincpg, mincov, maxcov) {
    ## One treat-vs-control contrast: read only the two groups' tables,
    ## treatment vector control=0 / treat=1 (methylKit is pairwise).
    result_group <- paste0(treat, "_vs_", control)
    sub_design <- design[design$group %in% c(control, treat), , drop = FALSE]
    idx <- match(sub_design$sample_id, ids)
    treatment <- ifelse(sub_design$group == control, 0L, 1L)
    sample_ids <- make.names(sub_design$sample_id, unique = TRUE)  # legal, unique R names

    print(paste0("[", date(), "] ", result_group, ": reading ", length(idx),
                 " merged CpG tables (control=", control, ", treat=", treat, ")..."))
    obj <- methRead(location = as.list(files[idx]),
                    sample.id = as.list(sample_ids),
                    assembly = "bsseq",
                    treatment = treatment,
                    context = "CpG",
                    pipeline = "bismarkCoverage",
                    mincov = as.integer(mincov))
    ## filterByCoverage: lo.count repeats the read-time minimum (harmless),
    ## hi.count is the absolute coverage ceiling that removes repeat-cluster
    ## artifacts (hi.perc would be a percentile, not a count).
    filtered <- filterByCoverage(obj, lo.count = as.integer(mincov),
                                 hi.count = as.integer(maxcov))
    meth <- unite(filtered, destrand = FALSE)
    batch_values <- if (is_batch) as.factor(sub_design$batch) else NULL

    print(paste0("[", date(), "] ", result_group, ": per-CpG differential test..."))
    myDiff <- calc_diff(meth, batch_values, is_batch, Nthreads)
    prefix <- file.path(output_name, result_group)
    n_all   <- write_dm(myDiff, "all",   paste0(prefix, "_DMC_all.tsv"),   mindiff, qvalue)
    n_hyper <- write_dm(myDiff, "hyper", paste0(prefix, "_DMC_hyper.tsv"), mindiff, qvalue)
    n_hypo  <- write_dm(myDiff, "hypo",  paste0(prefix, "_DMC_hypo.tsv"),  mindiff, qvalue)

    print(paste0("[", date(), "] ", result_group, ": tiled DMR test (win=",
                 tilelen, ", step=", tilestep, ", cov.bases=", mincpg, ")..."))
    tiles <- tileMethylCounts(filtered, win.size = tilelen, step.size = tilestep,
                              cov.bases = as.integer(mincpg))
    tile_meth <- unite(tiles, destrand = FALSE)
    tile_diff <- calc_diff(tile_meth, batch_values, is_batch, Nthreads)
    n_tiles <- write_dm(tile_diff, "all", paste0(prefix, "_DMR_tiles.tsv"),
                        mindiff, qvalue)

    print(paste0("[", date(), "] ", result_group, ": DMC all=", n_all,
                 " hyper=", n_hyper, " hypo=", n_hypo, " DMR tiles=", n_tiles))
    list(contrast = result_group, dmc_all = n_all, dmc_hyper = n_hyper,
         dmc_hypo = n_hypo, dmr_tiles = n_tiles)
}

main <- function(files, sampleinfo, control, qvalue, mindiff, tilelen,
                 tilestep, mincpg, mincov, maxcov, is_batch, Nthreads,
                 output_name) {
    print(paste("[", date(), "] methylKit", as.character(packageVersion("methylKit")),
                "with", Nthreads, "thread(s).", sep = " "))
    design <- read_design(sampleinfo)
    has_group <- "group" %in% names(design)
    if (is_batch && !"batch" %in% names(design)) {
        stop("batch correction requested (-b T) but the sample table has no batch column")
    }

    ## Match the input tables to the sample-table rows by sample_id (derived
    ## from the coverage2cytosine file names); the contrast loop then uses the
    ## sample-table order.
    ids <- sub("\\.CpG_merged.*$", "", basename(files))
    if (anyDuplicated(ids) > 0) {
        stop("input file list contains duplicate sample ids: ",
             paste(ids[duplicated(ids)], collapse = ", "))
    }
    unknown <- !(ids %in% design$sample_id)
    if (any(unknown)) {
        stop("sample table has no row for input file(s): ",
             paste(files[unknown], collapse = ", "))
    }
    unmatched <- setdiff(design$sample_id, ids)
    if (length(unmatched) > 0) {
        stop("no merged CpG table for sample(s): ",
             paste(unmatched, collapse = ", "),
             " (run with methylation_extractor.merge_cpg: true)")
    }

    dir.create(output_name, recursive = TRUE, showWarnings = FALSE)
    treats <- if (has_group) unique(design$group[design$group != control]) else character(0)
    if (length(treats) == 0) {
        stop("no non-control group in the sample table: contrasts need at least ",
             "one group other than control_group=", control)
    }

    summary_rows <- do.call(rbind, lapply(treats, function(treat) {
        run_contrast(files, ids, design, control, treat, is_batch, Nthreads,
                     output_name, qvalue, mindiff, tilelen, tilestep,
                     mincpg, mincov, maxcov)
    }))
    summary <- data.frame(matrix(unlist(summary_rows), ncol = 5, byrow = TRUE),
                          stringsAsFactors = FALSE)
    names(summary) <- c("contrast", "dmc_all", "dmc_hyper", "dmc_hypo", "dmr_tiles")
    write.table(summary, file = file.path(output_name, "DMR_summary.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    writeLines(capture.output(sessionInfo()),
               file.path(output_name, "sessionInfo.txt"))
    print(paste("[", date(), "] All done!", sep = ""))
}

## call main function
pkgs <- c('methylKit', 'getopt')
lapply(pkgs, function(x){
   suppressMessages(library(x, character.only = T))})
spec <- matrix(c("files",     "f", 2, "character", "Input merged CpG tables, comma-separated .cov.gz files (bismark 6-column coverage layout).",
                 "sampleinfo","s", 2, "character", "Input sample table, [filename, comma separated csv file].",
                 "outdir",    "o", 2, "character", "Output directory for the DMR results.",
                 "control",   "r", 1, "character", "The control group name in the sample table, [optional, default is control].",
                 "qvalue",    "q", 1, "numeric",   "getMethylDiff qvalue cutoff (SLIM adjusted), [optional, default is 0.01].",
                 "mindiff",   "d", 1, "numeric",   "Minimum percent methylation difference (treatment - control), [optional, default is 25].",
                 "tilelen",   "l", 1, "numeric",   "DMR tiling window size in bp, [optional, default is 1000].",
                 "tilestep",  "k", 1, "numeric",   "DMR tiling step size in bp, [optional, default is 100].",
                 "mincpg",    "g", 1, "numeric",   "Minimum covered CpGs per tile (cov.bases), [optional, default is 3].",
                 "mincov",    "c", 1, "numeric",   "Minimum per-site coverage (read-time mincov + filterByCoverage lo.count), [optional, default is 10].",
                 "maxcov",    "m", 1, "numeric",   "Absolute per-site coverage ceiling (filterByCoverage hi.count), [optional, default is 500].",
                 "batch",     "b", 1, "logical",   "If considering the batch column as a covariate, True or False, [optional, default is False].",
                 "threads",   "t", 1, "numeric",   "Using CPU numbers (methylKit mc.cores), [optional, default is 1].",
                 "help",      "h", 0, "logical",   "Show this help information."),
                byrow = T, ncol = 5)
opt <- getopt(spec = spec)
print(opt)
# check
if (!is.null(opt$help) || is.null(opt$files) || is.null(opt$sampleinfo) || is.null(opt$outdir)) {
    cat(paste(getopt(spec = spec, usage = T), "\n"))
    quit()
}
if (is.null(opt$control))  { opt$control  <- as.character("control") }
if (is.null(opt$qvalue))   { opt$qvalue   <- as.numeric("0.01") }
if (is.null(opt$mindiff))  { opt$mindiff  <- as.numeric("25") }
if (is.null(opt$tilelen))  { opt$tilelen  <- as.numeric("1000") }
if (is.null(opt$tilestep)) { opt$tilestep <- as.numeric("100") }
if (is.null(opt$mincpg))   { opt$mincpg   <- as.numeric("3") }
if (is.null(opt$mincov))   { opt$mincov   <- as.numeric("10") }
if (is.null(opt$maxcov))   { opt$maxcov   <- as.numeric("500") }
if (is.null(opt$batch))    { opt$batch    <- as.logical("False") }
if (is.null(opt$threads))  { opt$threads  <- as.numeric("1") }
#
files        <- strsplit(as.character(opt$files), ",", fixed = TRUE)[[1]]
sample       <- as.character(opt$sampleinfo)
control      <- as.character(opt$control)
qvalue       <- as.numeric(opt$qvalue)
mindiff      <- as.numeric(opt$mindiff)
tilelen      <- as.numeric(opt$tilelen)
tilestep     <- as.numeric(opt$tilestep)
mincpg       <- as.numeric(opt$mincpg)
mincov       <- as.numeric(opt$mincov)
maxcov       <- as.numeric(opt$maxcov)
is_batch     <- as.logical(opt$batch)
Nthreads     <- as.numeric(opt$threads)
output_name  <- trimws(as.character(opt$outdir), which = c("both", "left", "right"),
                       whitespace = "[ \t\r\n]")

# running main function
main(files, sample, control, qvalue, mindiff, tilelen, tilestep, mincpg,
     mincov, maxcov, is_batch, Nthreads, output_name)
