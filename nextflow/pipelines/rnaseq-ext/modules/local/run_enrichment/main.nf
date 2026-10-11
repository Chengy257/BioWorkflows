process RUN_ENRICHMENT {
    tag "${meta.contrast}:${meta.direction}"
    container 'https://depot.galaxyproject.org/singularity/bioconductor-clusterprofiler:4.14.0--r44hdfd78af_0'

    input:
    tuple val(meta), path(up_list), path(down_list), path(all_list), path(ranked)
    path orgdb_files // possibly-empty list: the optional OrgDb tarball

    output:
    path "*.{tsv,pdf}"  , emit: enrichment
    path "versions.yml" , emit: versions

    script:
    def orgdb_file = orgdb_files ? "${orgdb_files[0]}" : ''
    def kegg_arg = meta.kegg_org ? "${meta.kegg_org}" : ''
    """
    run_enrichment.R \\
        ${meta.direction == 'UP' ? up_list : meta.direction == 'DOWN' ? down_list : all_list} \\
        . ${meta.direction} \\
        '${meta.orgdb}' \\
        '${meta.keytype}' \\
        '${kegg_arg}' \\
        '${orgdb_file}'

    cat > versions.yml <<EOF
"${task.process}":
  r-enrichment: "\$(Rscript -e 'cat(as.character(packageVersion("clusterProfiler")))' 2>/dev/null || echo NA)"
  orgdb: "${meta.orgdb}"
EOF
    """

    stub:
    """
    echo "" > ${meta.direction}_GO-BP_AllResults.tsv
    echo "" > versions.yml
    """
}
