include { HOMER_FINDMOTIFS } from '../../../modules/local/homer_findmotifs/main'

workflow HOMER_MOTIF_EXT {
    take:
    ch_consensus    // channel: [ val(meta: id=antibody), path(consensus_peaks.bed) ]

    main:
    ch_bg = params.homer_motif.background
        ? Channel.value([file(params.homer_motif.background, checkIfExists: true)])
        : Channel.value([])
    HOMER_FINDMOTIFS(
        ch_consensus,
        params.homer_motif.genome,
        params.homer_motif.size,
        params.homer_motif.extra,
        ch_bg,
    )

    emit:
    motifs = HOMER_FINDMOTIFS.out.motifs
}
