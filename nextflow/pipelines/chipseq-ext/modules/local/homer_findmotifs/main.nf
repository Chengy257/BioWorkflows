process HOMER_FINDMOTIFS {
    tag "${meta.id}"
    label 'process_medium'

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
