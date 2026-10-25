/*
 * featureCounts quantification on the adopted run's published STAR BAMs,
 * with optional strand inference (decision record D1: the legacy
 * featureCounts counting semantics, retained as a local stage).
 *
 * Strand handling (S2 convention, shared with lncrna-nf):
 * - strand='auto' (default): RSeQC infer_experiment on each BAM, parsed by
 *   PARSE_STRANDEDNESS into the legacy featureCounts codes (secondstrand ->
 *   1, firststrand -> 2, unstranded -> 0; the exact awk mapping of the
 *   legacy check_strandedness rule, ported verbatim).
 * - '0'/'1'/'2': used directly (never hardcoded - D7).
 * The code is mapped to the module-convention name (forward/reverse/
 * unstranded) and rides the meta as strandedness; the vendored
 * subread/featurecounts module derives -s from it as the SINGLE source of
 * truth (injecting a second -s via ext.args would strand-disable every
 * stranded project - getopt last-wins). -p comes from meta.single_end.
 *
 * Outputs: per-sample featureCounts files plus the merged gene-count matrix
 * (the authoritative matrix when quantification_route=featurecounts).
 */

include { SUBREAD_FEATURECOUNTS } from '../../../modules/nf-core/subread/featurecounts'
include { RSEQC_INFEREXPERIMENT } from '../../../modules/nf-core/rseqc/inferexperiment'
include { PARSE_STRANDEDNESS } from '../../../modules/local/parse_strandedness'
include { MERGE_FEATURECOUNTS } from '../../../modules/local/merge_featurecounts'

def strandName(code) {
    // legacy featureCounts code -> the module-convention strandedness name
    code == '1' ? 'forward' : code == '2' ? 'reverse' : 'unstranded'
}

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
                tuple(meta + [strandedness: strandName(strand_file.text.trim()),
                              single_end  : !params.featurecounts_quant.paired], bam)
            }
    } else {
        ch_meta_bams_stranded = ch_meta_bams.map { meta, bam, bai ->
            tuple(meta + [strandedness: strandName(params.featurecounts_quant.strand),
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
