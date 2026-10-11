process NOVEL_TRANSCRIPTS {
    tag "novel-transcripts"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.12.12'

    input:
    tuple val(meta), path(annotated_gtf)

    output:
    tuple val(meta), path("novel.transcripts.gtf"), emit: novel_gtf
    tuple val(meta), path("novel.transcriptIDs.txt"), emit: novel_ids

    script:
    """
    parse_classcode_u.py \\
        --gtf ${annotated_gtf} \\
        --out-gtf novel.transcripts.gtf \\
        --out-ids novel.transcriptIDs.txt
    """

    stub:
    """
    touch novel.transcripts.gtf novel.transcriptIDs.txt
    """
}
