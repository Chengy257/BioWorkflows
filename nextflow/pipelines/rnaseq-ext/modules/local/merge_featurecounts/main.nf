process MERGE_FEATURECOUNTS {
    tag "count_matrix"
    label 'process_low'
    container 'quay.io/biocontainers/python:3.12.12'

    input:
    val sample_ids
    path count_files,   stageAs: 'fc/*'
    path summary_files, stageAs: 'fc/*'

    output:
    path "count.matrix.tsv", emit: count_matrix
    path "count_stats.tsv"  , emit: count_stats
    path "sample_ids.txt"   , emit: sample_ids

    script:
    """
    cat > sample_ids.txt <<'IDS'
${sample_ids instanceof List ? sample_ids.join('\n') : sample_ids}
IDS
    merge_featurecounts.py fc sample_ids.txt
    mv fc/count.matrix.tsv fc/count_stats.tsv .
    """

    stub:
    """
    echo -e "gene_id" > count.matrix.tsv
    echo -e "sample" > count_stats.tsv
    echo "" > sample_ids.txt
    """
}
