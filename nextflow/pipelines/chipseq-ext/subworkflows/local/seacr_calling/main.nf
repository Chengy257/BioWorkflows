include { SAMTOOLS_MERGE } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_CTL } from '../../../modules/nf-core/samtools/merge/main'
include { DEEPTOOLS_BAMCOVERAGE as DEEPTOOLS_BDGC_RAW } from '../../../modules/nf-core/deeptools/bamcoverage/main'
include { DEEPTOOLS_BAMCOVERAGE as DEEPTOOLS_BDGC_RAW_CTL } from '../../../modules/nf-core/deeptools/bamcoverage/main'
include { SEACR_CALLPEAK } from '../../../modules/nf-core/seacr/callpeak/main'
include { SEACR_CONVERT } from '../../../modules/local/seacr_convert/main'

workflow SEACR_CALLING {
    take:
    ch_groups       // channel: [ val(meta: id=group, group, control), val(treat_metas) ]
    ch_bams         // channel: [ val(meta), path(bam) ]

    main:
    // All samples of a group, BAM + BAI side by side.
    ch_group_bams = ch_bams
        .map { meta, bam -> tuple(meta.group, meta.id, bam, file("${bam}.bai", checkIfExists: true)) }
        .groupTuple()

    ch_treat_ids = ch_groups.map { meta_g, treat_metas ->
        tuple(meta_g.id, treat_metas.collect { it.id }.toSet())
    }

    // One tuple per group: [meta, treat_bams, treat_bais, ctl_bams, ctl_bais]
    ch_prepared = ch_group_bams
        .join(ch_treat_ids)
        .flatMap { group, ids, bams, bais, treat_ids ->
            def treats = []
            def treats_bai = []
            def ctls = []
            def ctls_bai = []
            [ids, bams, bais].transpose().each { id, bam, bai ->
                if (treat_ids.contains(id)) {
                    treats << bam
                    treats_bai << bai
                } else {
                    ctls << bam
                    ctls_bai << bai
                }
            }
            def meta = [id: group, group: group]
            [tuple(meta, treats, treats_bai, ctls, ctls_bai)]
        }

    ch_fasta_stub = channel.value(tuple([:], [], [], []))
    ch_blacklist_stub = channel.value(tuple([:], []))

    // Pooled treat BAM -> raw-depth bedGraph (zero-omitted, bin 1), the
    // coverage SEACR integrates (legacy recipe).
    SAMTOOLS_MERGE(
        ch_prepared.map { meta, treats, treats_bai, ctls, ctls_bai ->
            tuple(meta, treats, treats_bai)
        },
        ch_fasta_stub,
        'bai',
    )
    SAMTOOLS_MERGE.out.bam.set { ch_treat_merged }
    DEEPTOOLS_BDGC_RAW(
        ch_treat_merged.map { meta, bam ->
            tuple(meta, bam, file("${bam}.bai", checkIfExists: true))
        },
        channel.value([]),
        channel.value([]),
        ch_blacklist_stub,
    )
    DEEPTOOLS_BDGC_RAW.out.bedgraph.set { ch_treat_bdg_raw }

    // Pooled control BAM, only for groups that declare a control.
    SAMTOOLS_MERGE_CTL(
        ch_prepared
            .flatMap { meta, treats, treats_bai, ctls, ctls_bai ->
                ctls ? [tuple([id: "${meta.id}_control", group: meta.group], ctls, ctls_bai)] : []
            },
        ch_fasta_stub,
        'bai',
    )
    SAMTOOLS_MERGE_CTL.out.bam.set { ch_ctl_merged }
    DEEPTOOLS_BDGC_RAW_CTL(
        ch_ctl_merged.map { meta, bam ->
            tuple(meta, bam, file("${bam}.bai", checkIfExists: true))
        },
        channel.value([]),
        channel.value([]),
        ch_blacklist_stub,
    )
    DEEPTOOLS_BDGC_RAW_CTL.out.bedgraph.set { ch_ctl_bdg_raw }

    // Join treat and control bedGraphs; groups without a control fall back
    // to the numeric FDR threshold inside SEACR_CALLPEAK.
    ch_treat_bdg = ch_treat_bdg_raw.map { meta, bdg -> tuple(meta.id, meta, bdg) }
    ch_ctl_bdg = ch_ctl_bdg_raw.map { meta, bdg -> tuple(meta.group, bdg) }
    ch_seacr_in = ch_treat_bdg
        .join(ch_ctl_bdg, remainder: true)
        .map { row ->
            row.size() == 4 ? tuple(row[1], row[2], row[3]) : tuple(row[1], row[2], [])
        }

    SEACR_CALLPEAK(ch_seacr_in, params.seacr.fdr_threshold)
    SEACR_CONVERT(
        SEACR_CALLPEAK.out.bed.map { meta, bed -> tuple(meta, bed, meta.id) },
    )

    emit:
    narrowpeak = SEACR_CONVERT.out.narrowpeak
    bedgraph   = ch_treat_bdg.map { id, meta, bdg -> tuple(meta, bdg) }
}
