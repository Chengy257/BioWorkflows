include { HOMER_FINDMOTIFS } from '../../../modules/local/homer_findmotifs/main'
include { HOMER_FINDMOTIFS_CUSTOM } from '../../../modules/local/homer_findmotifs/main'

workflow HOMER_MOTIF_EXT {
    take:
    ch_consensus    // channel: [ val(meta: id=antibody), path(consensus_peaks.bed) ]

    main:
    ch_bg = params.homer_motif.background
        ? Channel.value([file(params.homer_motif.background, checkIfExists: true)])
        : Channel.value([])
    // A genome value pointing at an existing file (custom FASTA) is staged as
    // a path input so it is visible inside the container; installed HOMER
    // genome tags stay string-valued.
    def genome_spec = params.homer_motif.genome
    def genome_file = genome_spec ? file(genome_spec) : null
    def ch_motifs
    if (genome_file && genome_file.isFile()) {
        HOMER_FINDMOTIFS_CUSTOM(ch_consensus, genome_file, params.homer_motif.size, params.homer_motif.extra, ch_bg)
        ch_motifs = HOMER_FINDMOTIFS_CUSTOM.out.motifs
    } else {
        HOMER_FINDMOTIFS(ch_consensus, genome_spec, params.homer_motif.size, params.homer_motif.extra, ch_bg)
        ch_motifs = HOMER_FINDMOTIFS.out.motifs
    }

    emit:
    motifs = ch_motifs
}
