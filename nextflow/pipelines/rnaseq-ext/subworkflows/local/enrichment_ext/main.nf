/*
 * Enrichment stage (decision record D3): clusterProfiler GO/KEGG
 * over-representation per contrast and direction (UP/DOWN/ALL DEG lists),
 * plus GSEA on the full ranked gene list. Species-neutral by construction:
 * the OrgDb package name, keyType and KEGG code are all params.
 */

include { PREPARE_ENRICHMENT_INPUTS } from '../../../modules/local/prepare_enrichment_inputs'
include { RUN_ENRICHMENT } from '../../../modules/local/run_enrichment'
include { RUN_GSEA } from '../../../modules/local/run_gsea'

workflow ENRICHMENT_EXT {
    take:
    ch_de_tables // one tuple(meta, de_table) per contrast (see workflows/ext.nf)

    main:
    PREPARE_ENRICHMENT_INPUTS(ch_de_tables)

    // One enrichment task per (contrast, direction) gene list. The optional
    // OrgDb tarball rides as a possibly-empty path list on a separate,
    // positionally-aligned input (Nextflow 26's replacement for optional
    // inputs; one entry per emitted list).
    ch_lists = PREPARE_ENRICHMENT_INPUTS.out.gene_lists
    ch_gene_lists = ch_lists.flatMap { meta, up, down, all, ranked ->
        [
            tuple(meta + [direction: 'UP'], up, down, all, ranked),
            tuple(meta + [direction: 'DOWN'], up, down, all, ranked),
            tuple(meta + [direction: 'ALL'], up, down, all, ranked),
        ]
    }
    ch_orgdb_lists = ch_lists.flatMap { meta, up, down, all, ranked ->
        def f = meta.orgdb_file ? [meta.orgdb_file] : []
        [f, f, f]
    }
    RUN_ENRICHMENT(ch_gene_lists, ch_orgdb_lists)

    RUN_GSEA(ch_lists.map { meta, up, down, all, ranked -> tuple(meta + [direction: 'gsea'], ranked) },
             ch_lists.map { meta, up, down, all, ranked -> meta.orgdb_file ? [meta.orgdb_file] : [] })

    emit:
    gene_lists = ch_gene_lists
    enrichment = RUN_ENRICHMENT.out.enrichment
    gsea       = RUN_GSEA.out.gsea
}
