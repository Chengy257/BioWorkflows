process MERGE_LNCRNA_COUNTS {
    tag "merge-counts"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.12.12'

    input:
    tuple val(meta), path(counts, stageAs: 'counts/*')
    val sample_ids

    output:
    tuple val(meta), path("lncrna.count.matrix.tsv")  , emit: counts
    tuple val(meta), path("lncrna.tpm.matrix.tsv")    , emit: tpm
    tuple val(meta), path("lncrna.fpkm.matrix.tsv")   , emit: fpkm
    tuple val(meta), path("lncrna.count_stats.tsv")   , emit: stats

    script:
    """
    echo '${sample_ids}' | tr ',' '\\n' > ids.txt
    merge_lncrna_counts.py counts ids.txt 'lncrna.'
    """

    stub:
    """
    touch lncrna.count.matrix.tsv lncrna.tpm.matrix.tsv lncrna.fpkm.matrix.tsv lncrna.count_stats.tsv
    """
}
