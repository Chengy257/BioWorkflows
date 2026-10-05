process IDR_UNION {
    tag "${group}"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.12.12'

    input:
    tuple val(group), path(idr_files, stageAs: 'pairs/*')
    val min_replicates

    output:
    path "${group}_idr_consensus.narrowPeak", emit: narrowpeak

    script:
    """
    idr_union.py \\
        --min-replicates ${min_replicates} \\
        --out ${group}_idr_consensus.narrowPeak \\
        pairs/*
    """

    stub:
    """
    touch ${group}_idr_consensus.narrowPeak
    """
}
