include { DIFFBIND_SHEET } from '../../../modules/local/diffbind_sheet/main'
include { RUN_DIFFBIND } from '../../../modules/local/run_diffbind/main'

workflow DIFFBIND_EXT {
    take:
    ch_bams         // channel: [ val(meta), path(bam) ]
    ch_peaks        // channel: [ val(meta), path(peak) ]

    main:
    def peak_suffix = "_peaks.${params.ext_peak_type == 'narrow_peak' ? 'narrowPeak' : 'broadPeak'}"
    def peak_type = params.ext_peak_type == 'narrow_peak' ? 'narrowpeak' : 'broadpeak'
    ch_contrasts = Channel.fromList(params.diffbind.contrasts)
        .map { c -> [id: c, contrast: c] }
    ch_sheet_file = Channel.value(file(params.ext_sample_sheet, checkIfExists: true))
    ch_bam_files = ch_bams.map { meta, bam -> bam }.collect()
    ch_peak_files = ch_peaks.map { meta, peak -> peak }.collect()

    DIFFBIND_SHEET(
        ch_contrasts,
        ch_sheet_file,
        ch_bam_files,
        ch_peak_files,
        params.diffbind.use_controls,
        params.ext_bam_suffix,
        peak_suffix,
        peak_type,
    )
    RUN_DIFFBIND(
        DIFFBIND_SHEET.out.sheet,
        ch_bam_files,
        ch_peak_files,
        params.diffbind.analysis,
        params.diffbind.summit_flank,
        params.diffbind.fdr,
        params.diffbind.foldchange,
        params.diffbind.batch_correction,
    )

    emit:
    results     = RUN_DIFFBIND.out.results
    significant = RUN_DIFFBIND.out.significant
}
