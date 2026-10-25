process PREPARE_ENRICHMENT_INPUTS {
    tag "${meta.contrast}"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.12.12'

    input:
    tuple val(meta), path(de_table)

    output:
    tuple val(meta), path("${meta.contrast}_UP.DEGs.txt"), path("${meta.contrast}_DOWN.DEGs.txt"), \
          path("${meta.contrast}_ALL.DEGs.txt"), path("${meta.contrast}_ranked.tsv"), emit: gene_lists

    script:
    """
    prepare_enrichment_inputs.py \\
        ${de_table} . ${meta.contrast} \\
        --id-col '${meta.id_col}' \\
        --log2fc-col '${meta.log2fc_col}' \\
        --padj-col '${meta.padj_col}' \\
        --stat-col '${meta.stat_col}' \\
        --padj-max ${meta.padj_max instanceof CharSequence ? meta.padj_max : meta.padj_max} \\
        --log2fc-min ${meta.log2fc_min}
    """

    stub:
    """
    echo "" > ${meta.contrast}_UP.DEGs.txt
    echo "" > ${meta.contrast}_DOWN.DEGs.txt
    echo "" > ${meta.contrast}_ALL.DEGs.txt
    echo -e "gene_id\\tstat" > ${meta.contrast}_ranked.tsv
    """
}
