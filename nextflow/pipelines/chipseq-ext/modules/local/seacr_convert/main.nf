process SEACR_CONVERT {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.12.12'

    input:
    tuple val(meta), path(bedgraph), val(prefix)

    output:
    tuple val(meta), path("${prefix}_peaks.narrowPeak"), emit: narrowpeak

    script:
    """
    seacr_to_narrowpeak.py ${bedgraph} ${prefix}_peaks.narrowPeak ${prefix}
    """

    stub:
    """
    touch ${prefix}_peaks.narrowPeak
    """
}
