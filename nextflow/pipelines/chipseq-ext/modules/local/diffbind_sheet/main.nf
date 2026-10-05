process DIFFBIND_SHEET {
    tag "${meta.contrast}"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.12.12'

    input:
    val   meta
    path  input_sheet
    path  bams, stageAs: 'bams/*'
    path  peaks, stageAs: 'peaks/*'
    val use_controls
    val bam_suffix
    val peak_suffix
    val peak_type

    output:
    tuple val(meta), path("${meta.contrast}_samplesheet.tsv"), emit: sheet

    script:
    def ctrl_flag = use_controls ? "T" : "F"
    """
    diffbind_sheet.py \\
        --contrast ${meta.contrast} \\
        --sheet ${input_sheet} \\
        --bam-suffix '${bam_suffix}' \\
        --peak-suffix '${peak_suffix}' \\
        --peak-type ${peak_type} \\
        --use-controls ${ctrl_flag} \\
        --out ${meta.contrast}_samplesheet.tsv \\
        bams/* peaks/*
    """

    stub:
    """
    touch ${meta.contrast}_samplesheet.tsv
    """
}
