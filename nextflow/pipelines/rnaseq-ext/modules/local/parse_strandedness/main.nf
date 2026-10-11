process PARSE_STRANDEDNESS {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.12.12'

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path("${meta.id}.strand.txt"), emit: strand
    path "${meta.id}.strand.txt"                  , emit: strand_file

    script:
    """
    parse_strandedness.py \\
        ${report} \\
        ${meta.id}.strand.txt
    """

    stub:
    """
    echo 0 > ${meta.id}.strand.txt
    """
}
