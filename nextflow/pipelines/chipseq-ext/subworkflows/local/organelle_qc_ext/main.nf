include { SAMTOOLS_IDXSTATS } from '../../../modules/nf-core/samtools/idxstats/main'
include { ORGANELLE_SUMMARY } from '../../../modules/local/organelle_summary/main'

workflow ORGANELLE_QC_EXT {
    take:
    ch_bams   // channel: [ val(meta), path(bam) ]

    main:
    ch_with_idx = ch_bams.map { meta, bam ->
        tuple(meta, bam, file("${bam}.bai", checkIfExists: true))
    }
    SAMTOOLS_IDXSTATS(ch_with_idx)
    ORGANELLE_SUMMARY(params.organelle_qc.patterns, SAMTOOLS_IDXSTATS.out.idxstats.collect())

    emit:
    summary = ORGANELLE_SUMMARY.out.tsv
    multiqc = ORGANELLE_SUMMARY.out.mqc
}
