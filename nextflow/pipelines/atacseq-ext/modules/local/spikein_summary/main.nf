process SPIKEIN_SUMMARY {
    tag "spike_in"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.12.12'

    input:
    val  patterns
    val  name
    path idxstats

    output:
    path "spikein_summary.tsv"    , emit: tsv
    path "spikein_summary_mqc.tsv", emit: mqc

    script:
    """
    spikein_summary.py \\
        --patterns '${patterns}' \\
        --name '${name}' \\
        --out spikein_summary.tsv \\
        --mqc spikein_summary_mqc.tsv \\
        ${idxstats}
    """

    stub:
    """
    touch spikein_summary.tsv spikein_summary_mqc.tsv
    """
}
