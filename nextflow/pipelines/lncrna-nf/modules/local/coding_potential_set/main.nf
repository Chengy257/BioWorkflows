process CODING_POTENTIAL_SET {
    tag "coding-potential"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.12.12'

    input:
    tuple val(meta), path(cpc2_out)
    tuple val(meta2), path(cpat_out)
    tuple val(meta3), path(novel_gtf)
    val min_length

    output:
    tuple val(meta), path("noncoding.transcriptIDs.txt"), emit: noncoding_ids
    tuple val(meta), path("noncoding.transcripts.gtf")  , emit: noncoding_gtf

    script:
    """
    coding_potential_intersect.py \\
        --cpc2 ${cpc2_out} \\
        --cpat ${cpat_out} \\
        --gtf ${novel_gtf} \\
        --min-length ${min_length} \\
        --out-ids noncoding.transcriptIDs.txt \\
        --out-gtf noncoding.transcripts.gtf
    """

    stub:
    """
    touch noncoding.transcriptIDs.txt noncoding.transcripts.gtf
    """
}
