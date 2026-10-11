/*
 * featureCounts quantification on the adopted run's published STAR BAMs,
 * with optional strand inference (decision record D1: the legacy
 * featureCounts counting semantics, retained as a local stage).
 *
 * Strand handling:
 * - strand='auto' (default): RSeQC infer_experiment on each BAM, parsed by
 *   PARSE_STRANDEDNESS into the legacy featureCounts codes (secondstrand ->
 *   1, firststrand -> 2, unstranded -> 0; the exact awk mapping of the
 *   legacy check_strandedness rule, ported verbatim).
 * - '0'/'1'/'2': used directly (never hardcoded - D7).
 * The strand code rides the meta as strand_code and reaches featureCounts
 * through the SUBREAD_FEATURECOUNTS ext.args closure (nextflow.config);
 * -p/-s come from meta.single_end (module convention) - never hardcoded.
 *
 * Outputs: per-sample featureCounts files plus the merged gene-count matrix
 * (the authoritative matrix when quantification_route=featurecounts).
 */

include { SUBREAD_FEATURECOUNTS } from '../../../modules/nf-core/subread/featurecounts'
include { RSEQC_INFEREXPERIMENT } from '../../../modules/nf-core/rseqc/inferexperiment'
include { PARSE_STRANDEDNESS } from '../../../modules/local/parse_strandedness'
include { MERGE_FEATURECOUNTS } from '../../../modules/local/merge_featurecounts'

workflow FEATURECOUNTS_QUANT {
    take:
    ch_meta_bams // tuple(meta, bam, bai)
    ch_gtf       // value channel: the reference GTF
    ch_bed       // value channel: gene BED, or null when strand != 'auto'

    main:
    if (params.featurecounts_quant.strand == 'auto') {
        RSEQC_INFEREXPERIMENT(ch_meta_bams, ch_bed)
        PARSE_STRANDEDNESS(RSEQC_INFEREXPERIMENT.out.txt)
        ch_meta_bams_stranded = ch_meta_bams
            .map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }
            .join(PARSE_STRANDEDNESS.out.strand.map { meta, f -> tuple(meta.id, f) }, failOnMismatch: true)
            .map { id, meta, bam, bai, strand_file ->
                tuple(meta + [strand_code: strand_file.text.trim(), single_end: !params.featurecounts_quant.paired], bam)
            }
    } else {
        ch_meta_bams_stranded = ch_meta_bams.map { meta, bam, bai ->
            tuple(meta + [strand_code: params.featurecounts_quant.strand,
                          single_end  : !params.featurecounts_quant.paired], bam)
        }
    }

    // The module declares ONE tuple input: tuple(meta, bams, annotation);
    // the value-channel combine appends the GTF as the third element.
    SUBREAD_FEATURECOUNTS(ch_meta_bams_stranded.combine(ch_gtf))

    MERGE_FEATURECOUNTS(
        SUBREAD_FEATURECOUNTS.out.counts.map { meta, f -> meta.id }.collect(),
        SUBREAD_FEATURECOUNTS.out.counts.map { meta, f -> f }.toList().map { files -> files.sort { it.name } },
        SUBREAD_FEATURECOUNTS.out.summary.map { meta, f -> f }.toList().map { files -> files.sort { it.name } }
    )

    emit:
    counts        = SUBREAD_FEATURECOUNTS.out.counts
    strandedness  = params.featurecounts_quant.strand == 'auto' ? PARSE_STRANDEDNESS.out.strand : channel.empty()
    count_matrix  = MERGE_FEATURECOUNTS.out.count_matrix
    count_stats   = MERGE_FEATURECOUNTS.out.count_stats
}
