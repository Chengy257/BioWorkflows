process HOMER_FINDMOTIFS {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/homer:5.1--pl5321hc52dbad_1'

    input:
    tuple val(meta), path(bed)
    val  genome
    val  size
    val  extra_args
    path bg

    output:
    path "${meta.id}_motifs", emit: motifs
    path "${meta.id}_motifs.log", emit: log

    script:
    def bg_arg = bg ? "-bg ${bg}" : ''
    """
    findMotifsGenome.pl \\
        ${bed} \\
        ${genome} \\
        ${meta.id}_motifs \\
        -size ${size} \\
        ${bg_arg} \\
        ${extra_args} \\
        > ${meta.id}_motifs.log 2>&1
    """

    stub:
    """
    mkdir -p ${meta.id}_motifs
    touch ${meta.id}_motifs.log
    """
}

// Custom-FASTA variant: the genome is staged as a path input so it is
// visible inside the container (a val string pointing at a host path is
// not). Installed HOMER genome tags stay on HOMER_FINDMOTIFS.
process HOMER_FINDMOTIFS_CUSTOM {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/homer:5.1--pl5321hc52dbad_1'

    input:
    tuple val(meta), path(bed)
    path genome_fasta
    val  size
    val  extra_args
    path bg

    output:
    path "${meta.id}_motifs", emit: motifs
    path "${meta.id}_motifs.log", emit: log

    script:
    def bg_arg = bg ? "-bg ${bg}" : ''
    """
    findMotifsGenome.pl \\
        ${bed} \\
        ${genome_fasta} \\
        ${meta.id}_motifs \\
        -size ${size} \\
        ${bg_arg} \\
        ${extra_args} \\
        > ${meta.id}_motifs.log 2>&1
    """

    stub:
    """
    mkdir -p ${meta.id}_motifs
    touch ${meta.id}_motifs.log
    """
}
