include { IDR } from '../../../modules/nf-core/idr/main'
include { IDR_UNION } from '../../../modules/local/idr_union/main'

workflow IDR_REPRODUCIBILITY {
    take:
    ch_pairs        // channel: [ val(meta: id=<group>__<a>_vs_<b>, group=group), path(peak_a), path(peak_b) ]

    main:
    IDR(ch_pairs.map { meta, peak_a, peak_b -> tuple(meta, [peak_a, peak_b], 'narrowPeak') })
    IDR_UNION(
        IDR.out.idr
            .map { meta, idr_file -> tuple(meta.group, idr_file) }
            .groupTuple()
            .map { group, files -> tuple(group, files) },
        params.idr.min_replicates,
    )

    emit:
    narrowpeak = IDR_UNION.out.narrowpeak
    pairwise   = IDR.out.idr
}
