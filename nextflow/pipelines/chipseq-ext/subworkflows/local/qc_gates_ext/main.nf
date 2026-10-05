include { QC_GATES_SUMMARY } from '../../../modules/local/qc_gates/main'

workflow QC_GATES_EXT {
    take:
    ch_flagstats   // channel: [ val(sample_id), path(flagstat) ]
    ch_dup         // channel: [ val(sample_id), path(metrics) ] (may be empty)
    ch_org         // channel: [ path(organelle_summary.tsv) ] (single task, or empty)

    main:
    QC_GATES_SUMMARY(
        params.qc_gates,
        params.ext_flagstat_suffix,
        params.ext_dup_suffix,
        ch_flagstats.map { id, flagstat -> flagstat }.collect(),
        ch_dup.map { id, metrics -> metrics }.collect(),
        ch_org,
    )

    emit:
    tsv = QC_GATES_SUMMARY.out.tsv
    multiqc = QC_GATES_SUMMARY.out.mqc
}
