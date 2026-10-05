process RUN_DIFFBIND {
    tag "${meta.contrast}"
    label 'process_high'
    container 'quay.io/biocontainers/bioconductor-diffbind:3.20.0--r45ha27e39d_0'

    input:
    tuple val(meta), path(samplesheet)
    path bams, stageAs: 'bams/*'
    path peaks, stageAs: 'peaks/*'
    val analysis
    val summit_flank
    val fdr
    val foldchange
    val batch_correction

    output:
    tuple val(meta), path("DB_results.tsv"), emit: results
    tuple val(meta), path("DB_significant.tsv"), emit: significant
    tuple val(meta), path("MA_plot.png"), emit: ma_plot
    tuple val(meta), path("Volcano_plot.png"), emit: volcano_plot
    tuple val(meta), path("PCA_plot.png"), emit: pca_plot
    tuple val(meta), path("sessionInfo.txt"), emit: session_info
    tuple val("${task.process}"), val('diffbind'), eval("Rscript -e 'cat(as.character(packageVersion(\"DiffBind\")))'"), topic: versions, emit: versions_diffbind

    script:
    def block_flag = batch_correction ? "T" : "F"
    """
    run_diffbind.R \\
        ${samplesheet} \\
        ${analysis} \\
        ${summit_flank} \\
        ${fdr} \\
        ${foldchange} \\
        ${block_flag}
    """

    stub:
    """
    touch DB_results.tsv DB_significant.tsv MA_plot.png Volcano_plot.png PCA_plot.png sessionInfo.txt
    """
}
