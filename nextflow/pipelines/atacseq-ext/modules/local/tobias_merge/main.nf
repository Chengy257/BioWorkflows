process TOBIAS_MERGE_GROUPS {
    tag "${meta.id}"
    label 'process_medium'
    // Local module (not central samtools/merge): the modern central interface
    // requires fasta/fai/gzi inputs and CRAM handling this stage does not
    // need - a plain merge + index of the group's treat BAMs is the whole job
    // (legacy chip v0.6 footprint stage semantics). Same samtools image the
    // vendored samtools/idxstats module uses.
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    tuple val(meta), path(bams, stageAs: 'group_bams/*')

    output:
    tuple val(meta), path("${meta.id}.pooled.bam"), path("${meta.id}.pooled.bam.bai"), emit: merged

    script:
    """
    samtools merge \\
        -f \\
        -@ ${task.cpus} \\
        ${meta.id}.pooled.bam \\
        group_bams/*
    samtools index ${meta.id}.pooled.bam
    """

    stub:
    """
    touch ${meta.id}.pooled.bam ${meta.id}.pooled.bam.bai
    """
}
