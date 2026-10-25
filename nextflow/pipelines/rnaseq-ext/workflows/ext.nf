/*
 * EXT - resolve the adopted nf-core/rnaseq run's published outputs and wire
 * the default-off extension stages.
 *
 * Two files are parsed synchronously in the workflow body (project files,
 * not pipeline config):
 * - ext_sample_sheet: the adopted run's nf-core/rnaseq sample sheet (header
 *   sample,fastq_1,fastq_2,...). Only the sample column is consumed (BAM
 *   resolution); existence is validated, content is not re-parsed.
 * - ext_metadata: analysis metadata CSV with header sample,group[,batch]
 *   (one row per sample; feeds the DA adapter's contrast derivation and the
 *   enrichment stage grouping). 'sample' must match the adopted sample ids.
 *
 * The quantification route switch (decision record D1) selects the
 * authoritative gene matrix: 'salmon' resolves the adopted run's native
 * merged matrix; 'featurecounts' forces the featureCounts stage on and
 * emits its matrix. EXT().out.quant_matrix carries the authoritative
 * matrix, so downstream consumers never pick the wrong file.
 */

include { FEATURECOUNTS_QUANT } from '../subworkflows/local/featurecounts_quant'
include { ENRICHMENT_EXT } from '../subworkflows/local/enrichment_ext'

def require_params() {
    if (!params.ext_results_dir) {
        error "params.ext_results_dir is required: point it at the adopted nf-core/rnaseq run's --outdir"
    }
    if (!params.ext_sample_sheet) {
        error "params.ext_sample_sheet is required: the adopted run's sample sheet (sample,fastq_1,...)"
    }
    if (!params.ext_metadata) {
        error "params.ext_metadata is required: CSV with sample,group[,batch]"
    }
    if (!['salmon', 'featurecounts'].contains(params.quantification_route)) {
        error "params.quantification_route must be 'salmon' or 'featurecounts', got '${params.quantification_route}'"
    }
    if (params.quantification_route == 'featurecounts' && !params.ext_gtf) {
        error "quantification_route=featurecounts requires params.ext_gtf (the reference GTF used for counting)"
    }
    if (params.featurecounts_quant.enabled && params.featurecounts_quant.strand == 'auto' && !params.ext_gene_bed) {
        error "featurecounts_quant.strand='auto' requires params.ext_gene_bed (RSeQC infer_experiment reference BED); or set strand explicitly to '0', '1' or '2'"
    }
    if (params.enrichment.enabled) {
        if (!params.enrichment.de_dir) {
            error "enrichment is enabled but params.enrichment.de_dir is empty (the differentialabundance outdir's tables/differential directory)"
        }
        if (!params.enrichment.orgdb) {
            error "enrichment is enabled but params.enrichment.orgdb is empty (OrgDb R package name - the pipeline never binds a species)"
        }
    }
}

workflow EXT {

    main:
    require_params()

    def results = params.ext_results_dir.replaceAll('/+$', '')
    def aligner = params.ext_aligner

    // Defaults mirror the pinned nf-core/rnaseq 3.27.0 publish layout
    // (conf/modules/*.config).
    def bam_dir = params.ext_bam_dir ?: "${results}/${aligner}"
    def salmon_dir = params.ext_salmon_dir ?: "${results}/salmon"

    def bam_suffix = params.ext_bam_suffix

    // ---- Adopted sample sheet (synchronous parse; sample ids only) ----------
    def sheet_lines = file(params.ext_sample_sheet, checkIfExists: true)
        .readLines()
        .findAll { it.trim() }
    if (sheet_lines.size() < 2) {
        error "sample sheet ${params.ext_sample_sheet} has no data rows"
    }
    def sheet_header = sheet_lines[0].split(',').collect { it.trim() }
    if (!sheet_header.contains('sample')) {
        error "sample sheet is missing the required column 'sample' (header: ${sheet_header.join(',')})"
    }
    def sample_idx = sheet_header.indexOf('sample')
    def fastq2_idx = sheet_header.indexOf('fastq_2')
    def adopted_ids = []
    def paired_seen = false
    sheet_lines.tail().each { line ->
        def id = line.split(',')[sample_idx]?.trim()
        if (!id) {
            error "sample sheet row missing sample id: ${line}"
        }
        if (id in adopted_ids) {
            error "duplicate sample in adopted sample sheet: ${id}"
        }
        adopted_ids << id
        if (fastq2_idx >= 0) {
            def cells = line.split(',')
            if (fastq2_idx < cells.size() && cells[fastq2_idx].trim()) {
                paired_seen = true
            }
        }
    }

    // ---- Analysis metadata sheet (synchronous plain-Groovy parse) -----------
    def meta_lines = file(params.ext_metadata, checkIfExists: true)
        .readLines()
        .findAll { it.trim() }
    if (meta_lines.size() < 2) {
        error "metadata sheet ${params.ext_metadata} has no data rows"
    }
    def meta_header = meta_lines[0].split(',').collect { it.trim() }
    if (!meta_header.contains('sample') || !meta_header.contains('group')) {
        error "metadata sheet must have 'sample' and 'group' columns (header: ${meta_header.join(',')})"
    }
    def metaById = [:]
    meta_lines.tail().each { line ->
        def cells = line.split(',')
        def row = [:]
        meta_header.eachWithIndex { col, i -> row[col] = i < cells.size() ? cells[i].trim() : '' }
        def id = row.sample
        def group = row.group
        if (!id || !group) {
            error "metadata sheet row missing sample or group: ${line}"
        }
        if (metaById.containsKey(id)) {
            error "duplicate sample in metadata sheet: ${id}"
        }
        metaById[id] = [id: id, group: group, batch: (row.batch ?: '')]
    }
    def missing_meta = adopted_ids.findAll { !metaById.containsKey(it) }
    if (missing_meta) {
        error "metadata sheet is missing adopted sample(s): ${missing_meta.sort().join(', ')}"
    }

    // ---- Per-sample BAMs (strict: one BAM per adopted sample, no strays) ----
    def bam_dir_path = file(bam_dir)
    if (!bam_dir_path.exists()) {
        error "BAM directory not found: ${bam_dir}"
    }
    def bam_ids = []
    bam_dir_path.list().findAll { it.endsWith(bam_suffix) }.sort().each { name ->
        def id = name.substring(0, name.length() - bam_suffix.length())
        if (!(id in adopted_ids)) {
            error "found BAM ${name} with no adopted sample sheet row"
        }
        bam_ids << id
    }
    def missing_bams = adopted_ids.findAll { !(it in bam_ids) }
    if (missing_bams) {
        error "no BAM found for adopted sample(s) in ${bam_dir}: ${missing_bams.sort().join(', ')}"
    }
    ch_meta_bams = Channel.fromList(bam_ids.collect { id ->
        def bam = file("${bam_dir}/${id}${bam_suffix}", checkIfExists: true)
        def bai = file("${bam_dir}/${id}${bam_suffix}.bai", checkIfExists: true)
        tuple(metaById[id], bam, bai)
    })

    // ---- Stage: featureCounts quantification + route switch (D1) ------------
    def fc_route = params.quantification_route == 'featurecounts' || params.featurecounts_quant.enabled
    ch_quant_matrix = channel.empty()
    if (fc_route) {
        if (!params.ext_gtf) {
            error "the featureCounts stage is enabled but params.ext_gtf is empty"
        }
        if (paired_seen && !params.featurecounts_quant.paired) {
            log.warn "the sample sheet declares paired reads (fastq_2) but featurecounts_quant.paired=false - fragments will be counted as two single reads; set featurecounts_quant.paired=true unless that is intended"
        }
        FEATURECOUNTS_QUANT(
            ch_meta_bams,
            Channel.value(file(params.ext_gtf, checkIfExists: true)),
            params.featurecounts_quant.strand == 'auto'
                ? Channel.value(file(params.ext_gene_bed, checkIfExists: true))
                : null,
            Channel.value(adopted_ids),
        )
        ch_quant_matrix = FEATURECOUNTS_QUANT.out.count_matrix
    }

    // The authoritative matrix per the D1 switch: when route=salmon the
    // adopted run's native merged matrix is emitted (verified to exist so
    // downstream consumers never start from a missing authoritative file).
    if (params.quantification_route == 'salmon') {
        def salmon_matrix = file("${salmon_dir}/${params.ext_salmon_matrix}", checkIfExists: true)
        ch_quant_matrix = Channel.value(salmon_matrix)
    }

    // ---- Stage: enrichment on differentialabundance results (D3) ------------
    if (params.enrichment.enabled) {
        def de_path = file(params.enrichment.de_dir)
        if (!de_path.exists()) {
            error "enrichment is enabled but the DE directory does not exist: ${params.enrichment.de_dir}"
        }
        def de_tables = de_path.listFiles()
            .findAll { it.name.endsWith(params.enrichment.de_suffix) }
            .sort { it.name }
        if (!de_tables) {
            error "enrichment is enabled but no *${params.enrichment.de_suffix} tables were found under ${params.enrichment.de_dir}"
        }
        def orgdb_file = params.enrichment.orgdb_file
            ? file(params.enrichment.orgdb_file, checkIfExists: true)
            : null
        // ONE (meta, table) tuple per contrast: staging every table in every
        // task (a raw list would) breaks argparse for >= 2 contrasts.
        ENRICHMENT_EXT(
            Channel.fromList(de_tables).map { t ->
                tuple(
                    [
                        contrast  : t.name.substring(0, t.name.length() - params.enrichment.de_suffix.length()),
                        group     : '',
                        batch     : '',
                        orgdb     : params.enrichment.orgdb,
                        keytype   : params.enrichment.keytype,
                        kegg_org  : params.enrichment.kegg_organism,
                        padj_max  : params.enrichment.padj_max,
                        log2fc_min: params.enrichment.log2fc_min,
                        id_col    : params.enrichment.id_col,
                        log2fc_col: params.enrichment.log2fc_col,
                        padj_col  : params.enrichment.padj_col,
                        stat_col  : params.enrichment.stat_col,
                        orgdb_file: orgdb_file,
                    ],
                    t)
            }
        )
    }

    emit:
    quant_matrix = ch_quant_matrix
}
