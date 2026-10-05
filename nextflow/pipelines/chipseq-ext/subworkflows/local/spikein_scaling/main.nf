include { SAMTOOLS_IDXSTATS } from '../../../modules/nf-core/samtools/idxstats/main'
include { SPIKEIN_SUMMARY } from '../../../modules/local/spikein_summary/main'
include { SPIKEIN_RESCALE } from '../../../modules/local/spikein_rescale/main'

workflow SPIKEIN_SCALING {
    take:
    ch_bams         // channel: [ val(meta), path(bam) ] (combined-reference run)

    main:
    ch_with_idx = ch_bams.map { meta, bam ->
        tuple(meta, bam, file("${bam}.bai", checkIfExists: true))
    }
    SAMTOOLS_IDXSTATS(ch_with_idx)
    SPIKEIN_SUMMARY(
        params.spike_in.patterns,
        params.spike_in.name,
        SAMTOOLS_IDXSTATS.out.idxstats.collect(),
    )
    if (params.spike_in.rescale_bigwigs) {
        ch_idx_by_meta = SAMTOOLS_IDXSTATS.out.idxstats.map { meta, idx -> tuple(meta.id, idx) }
        ch_rescale_in = ch_with_idx
            .map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }
            .join(ch_idx_by_meta)
            .map { id, meta, bam, bai, idx -> tuple(meta, bam, bai, idx) }
        SPIKEIN_RESCALE(ch_rescale_in, params.spike_in.patterns)
    }

    emit:
    summary = SPIKEIN_SUMMARY.out.tsv
    multiqc = SPIKEIN_SUMMARY.out.mqc
    bigwigs = (params.spike_in.rescale_bigwigs ? SPIKEIN_RESCALE.out.bigwig : channel.empty())
}
