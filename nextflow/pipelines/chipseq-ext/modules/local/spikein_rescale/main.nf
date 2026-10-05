process SPIKEIN_RESCALE {
    tag "${meta.id}"
    label 'process_medium'
    // Same mulled samtools+deeptools image the ecosystem bamCoverage module
    // uses; a local process is required here because the scale factor is a
    // per-sample runtime value that the module's static ext.args cannot
    // express.
    container 'quay.io/biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0'

    input:
    tuple val(meta), path(bam), path(bai), path(idxstats)
    val patterns

    output:
    tuple val(meta), path("${meta.id}.spikein_scaled.bigWig"), emit: bigwig

    script:
    """
    factor=\$(spikein_factor.py --patterns '${patterns}' ${idxstats} | cut -f2)
    bamCoverage \\
        -b ${bam} \\
        -o ${meta.id}.spikein_scaled.bigWig \\
        --normalizeUsing None \\
        --scaleFactor \${factor} \\
        --binSize 25
    """

    stub:
    """
    touch ${meta.id}.spikein_scaled.bigWig
    """
}
