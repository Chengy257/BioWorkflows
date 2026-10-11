/*
 * FOOTPRINT_EXT - per-group TOBIAS footprinting chain (legacy chip v0.6
 * stage semantics, ported onto the central nf-core tobias modules).
 *
 * Per group: pool the treat BAMs -> TOBIAS ATACorrect (Tn5 bias correction
 * over the consensus peaks) -> TOBIAS ScoreBigwig (footprint scores).
 * Then ONE TOBIAS BINDetect over all groups' footprint bigWigs (needs >= 2
 * groups; enforced loudly in workflows/ext.nf).
 *
 * Decision-record D2 override: nf-core/atacseq ships no TOBIAS of its own
 * (verified by source inspection 2026-10-12); this port vendors the central
 * nf-core tobias module code verbatim as local modules (containers
 * tobias:0.17.5) - stored under modules/local/ because this development box
 * intermittently serves stale file contents to new processes on the
 * nf-core/ paths it first created with broken fetches (verified by md5 on
 * both sides; re-verify on the HPC at the deferred execution round). Known
 * trap carried from the legacy stage: BINDetect can exit on a sklearn NaN
 * with degenerate (toy/degenerate) signal - re-validate on real research
 * data.
 */

include { TOBIAS_MERGE_GROUPS } from '../../../modules/local/tobias_merge/main'
include { TOBIAS_ATACORRECT  } from '../../../modules/local/tobias_atacorrect/main'
include { TOBIAS_SCOREBIGWIG  } from '../../../modules/local/tobias_scorebigwig/main'
include { TOBIAS_BINDETECT    } from '../../../modules/local/tobias_bindetect/main'

workflow FOOTPRINT_EXT {
    take:
    ch_groups_bams   // channel: [ val(meta), path(treat bams) ] - one entry per group
    consensus_bed    // path: consensus peak BED (single file, shared by all groups)
    genome_fasta     // path: genome FASTA (sibling .fai must exist)
    motifs_pfm       // path: JASPAR PFM motif file

    main:
    TOBIAS_MERGE_GROUPS(ch_groups_bams)

    TOBIAS_ATACORRECT(
        TOBIAS_MERGE_GROUPS.out.merged.map { meta, bam, bai ->
            tuple(meta, bam, bai, consensus_bed, genome_fasta)
        }
    )

    TOBIAS_SCOREBIGWIG(
        TOBIAS_ATACORRECT.out.corrected.map { meta, bw ->
            tuple(meta, bw, consensus_bed)
        }
    )

    TOBIAS_BINDETECT(
        TOBIAS_SCOREBIGWIG.out.footprints
            .map { meta, bw -> bw }
            .toList()
            .map { bws -> tuple([id: 'footprint_bindetect'], bws, consensus_bed, motifs_pfm, genome_fasta) }
    )

    emit:
    corrected  = TOBIAS_ATACORRECT.out.corrected
    footprints = TOBIAS_SCOREBIGWIG.out.footprints
    bindetect  = TOBIAS_BINDETECT.out.outdir
}
