process ORGANELLE_SUMMARY {
    tag "organelle_qc"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.12.12'

    input:
    val   patterns
    path  idxstats

    output:
    path "organelle_summary.tsv"     , emit: tsv
    path "organelle_summary_mqc.tsv" , emit: mqc

    script:
    """
    organelle_summary.py \\
        --patterns '${patterns}' \\
        --out organelle_summary.tsv \\
        --mqc organelle_summary_mqc.tsv \\
        ${idxstats}
    """

    stub:
    """
    touch organelle_summary.tsv organelle_summary_mqc.tsv
    """
}
