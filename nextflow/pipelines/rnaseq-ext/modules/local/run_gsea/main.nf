process RUN_GSEA {
    tag "${meta.contrast}"
    container 'https://depot.galaxyproject.org/singularity/bioconductor-clusterprofiler:4.14.0--r44hdfd78af_0'

    input:
    tuple val(meta), path(ranked)
    path orgdb_files // possibly-empty list: the optional OrgDb tarball

    output:
    path "*.{tsv,pdf}"  , emit: gsea
    path "versions.yml" , emit: versions

    script:
    def orgdb_file = orgdb_files ? "${orgdb_files[0]}" : ''
    """
    run_gsea.R \\
        ${ranked} \\
        . ${meta.contrast} \\
        '${meta.orgdb}' \\
        '${meta.keytype}' \\
        '${orgdb_file}'

    cat > versions.yml <<EOF
"${task.process}":
  r-gsea: "\$(Rscript -e 'cat(as.character(packageVersion("clusterProfiler")))' 2>/dev/null || echo NA)"
  orgdb: "${meta.orgdb}"
EOF
    """

    stub:
    """
    echo "" > ${meta.contrast}_BP_multiGSEA_ResultTable_all.tsv
    echo "" > versions.yml
    """
}
