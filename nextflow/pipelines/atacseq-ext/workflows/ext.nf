/*
 * EXT - resolve the adopted nf-core/atacseq run's published outputs and wire
 * the default-off extension stages.
 *
 * The sample sheet is parsed synchronously in the workflow body (it is a
 * project file, not pipeline config): CSV with header
 * sample_id,group,condition,control,batch. The control column names the
 * group's control sample (empty when the group has none); one control may
 * serve several treat groups (shared IgG, legacy semantics). Controls may be
 * their own rows - they are processed but excluded from the footprint group
 * pooling.
 *
 * The consumed published layout (pinned atacseq dev HEAD 6a930770,
 * conf/modules.config): filtered BAMs + BAI under {aligner}/merged_library
 * (unconditional), samtools stats/flagstat/idxstats under
 * {aligner}/merged_library/samtools_stats (clN prefix, unconditional),
 * picard metrics under {aligner}/merged_library/picard_metrics, MACS3 peaks
 * under {aligner}/merged_library/macs3/{peak_type} and the consensus BED
 * under .../consensus/.
 */

include { ORGANELLE_QC_EXT } from '../subworkflows/local/organelle_qc_ext'
include { QC_GATES_EXT } from '../subworkflows/local/qc_gates_ext'
include { SPIKEIN_SCALING } from '../subworkflows/local/spikein_scaling'
include { FOOTPRINT_EXT } from '../subworkflows/local/footprint_ext'

def require_params() {
    if (!params.ext_results_dir) {
        error "params.ext_results_dir is required: point it at the adopted nf-core/atacseq run's --outdir"
    }
    if (!params.ext_sample_sheet) {
        error "params.ext_sample_sheet is required: CSV with sample_id,group,condition,control,batch"
    }
    if (params.footprint.enabled) {
        if (!params.footprint.genome) {
            error "footprint is enabled but params.footprint.genome is empty (reference FASTA; a sibling .fai must exist)"
        }
        if (!params.footprint.motifs) {
            error "footprint is enabled but params.footprint.motifs is empty (JASPAR PFM motif file)"
        }
    }
}

workflow EXT {

    main:
    require_params()

    def results = params.ext_results_dir.replaceAll('/+$', '')
    def aligner = params.ext_aligner

    // Defaults mirror the pinned nf-core/atacseq dev HEAD publish layout
    // (conf/modules.config).
    def bam_dir = params.ext_bam_dir ?: "${results}/${aligner}/merged_library"
    def flagstat_dir = params.ext_flagstat_dir ?: "${results}/${aligner}/merged_library/samtools_stats"
    def consensus_bed = params.ext_consensus_bed
        ?: "${results}/${aligner}/merged_library/macs3/${params.ext_peak_type}/consensus/consensus_peaks.mLb.clN.bed"

    def bam_suffix = params.ext_bam_suffix
    def flagstat_suffix = params.ext_flagstat_suffix

    // ---- Sample sheet (synchronous plain-Groovy parse) ----------------------
    def csv_lines = file(params.ext_sample_sheet, checkIfExists: true)
        .readLines()
        .findAll { it.trim() }
    if (csv_lines.size() < 2) {
        error "sample sheet ${params.ext_sample_sheet} has no data rows"
    }
    def csv_header = csv_lines[0].split(',').collect { it.trim() }
    ['sample_id', 'group'].each { col ->
        if (!csv_header.contains(col)) {
            error "sample sheet is missing the required column '${col}' (header: ${csv_header.join(',')})"
        }
    }
    def metaById = [:]
    def groups = [:]
    csv_lines.tail().each { line ->
        def cells = line.split(',')
        def row = [:]
        csv_header.eachWithIndex { col, i -> row[col] = i < cells.size() ? cells[i].trim() : '' }
        def id = row.sample_id
        def group = row.group
        if (!id || !group) {
            error "sample sheet row missing sample_id or group: ${line}"
        }
        if (metaById.containsKey(id)) {
            error "duplicate sample_id in sample sheet: ${id}"
        }
        metaById[id] = [
            id       : id,
            group    : group,
            condition: (row.condition ?: group),
            control  : (row.control ?: ''),
            batch    : (row.batch ?: ''),
        ]
        if (!groups.containsKey(group)) groups[group] = []
        groups[group] << id
    }
    def groupControl = [:]
    groups.each { group, ids ->
        def declared = ids.collect { metaById[it].control }.findAll { it }.unique()
        if (declared.size() > 1) {
            error "group ${group} declares more than one control: ${declared}"
        }
        groupControl[group] = declared ? declared[0] : ''
        if (groupControl[group] && !metaById.containsKey(groupControl[group])) {
            error "group ${group} declares control '${groupControl[group]}' which is not a sample_id in the sheet"
        }
    }

    // ---- Per-sample BAMs (strict: one BAM per sheet sample, no strays).
    // Resolved synchronously: validation happens before any task is wired.
    def bam_dir_path = file(bam_dir)
    if (!bam_dir_path.exists()) {
        error "BAM directory not found: ${bam_dir}"
    }
    def bam_ids = []
    bam_dir_path.list().findAll { it.endsWith(bam_suffix) }.sort().each { name ->
        def id = name.substring(0, name.length() - bam_suffix.length())
        if (!metaById.containsKey(id)) {
            error "found BAM ${name} with no sample sheet row"
        }
        bam_ids << id
    }
    def missing_bams = metaById.keySet().findAll { !(it in bam_ids) }
    if (missing_bams) {
        error "no BAM found for sheet sample(s) in ${bam_dir}: ${missing_bams.sort().join(', ')}"
    }
    ch_meta_bams = Channel.fromList(bam_ids.collect { id ->
        tuple(metaById[id], file("${bam_dir}/${id}${bam_suffix}", checkIfExists: true))
    })

    // ---- Optional QC inputs (empty when the glob matches nothing) -----------
    ch_flagstats = Channel
        .fromPath("${flagstat_dir}/*${flagstat_suffix}")
        .map { p -> tuple(p.name.substring(0, p.name.length() - flagstat_suffix.length()), p) }
    def dup_glob = params.ext_dup_glob
        ?: "${results}/${aligner}/merged_library/picard_metrics/*.metrics.txt"
    ch_dup_metrics = Channel.fromPath(dup_glob)
        .map { p -> tuple(p.name.endsWith(params.ext_dup_suffix)
            ? p.name.substring(0, p.name.length() - params.ext_dup_suffix.length()) : p.simpleName, p) }

    // ---- Stage: organelle QC ------------------------------------------------
    ch_organelle_summary = channel.empty()
    if (params.organelle_qc.enabled) {
        ORGANELLE_QC_EXT(ch_meta_bams)
        ch_organelle_summary = ORGANELLE_QC_EXT.out.summary
    }

    // ---- Stage: QC gate table ------------------------------------------------
    if (params.qc_gates.enabled) {
        def flagstat_dir_path = file(flagstat_dir)
        if (!flagstat_dir_path.exists()) {
            error "qc_gates is enabled but the flagstat directory does not exist: ${flagstat_dir}"
        }
        if (!flagstat_dir_path.list().findAll { it.endsWith(flagstat_suffix) }) {
            error "qc_gates is enabled but no *${flagstat_suffix} files were found under ${flagstat_dir}"
        }
        QC_GATES_EXT(
            ch_flagstats,
            ch_dup_metrics,
            ch_organelle_summary.toList(),
            metaById.keySet().sort().join(','),
        )
    }

    // ---- Stage: spike-in scaling (combined-reference runs) -------------------
    ch_spikein_summary = channel.empty()
    if (params.spike_in.enabled) {
        SPIKEIN_SCALING(ch_meta_bams)
        ch_spikein_summary = SPIKEIN_SCALING.out.summary
    }

    // ---- Stage: per-group TOBIAS footprinting --------------------------------
    if (params.footprint.enabled) {
        def footprint_group_list = []
        groups.each { group, ids ->
            def treats = ids.findAll { it != groupControl[group] }.sort()
            if (!treats) return // continue
            footprint_group_list << tuple(
                [id: group, group: group],
                treats.collect { id -> file("${bam_dir}/${id}${bam_suffix}", checkIfExists: true) },
            )
        }
        if (footprint_group_list.size() < 2) {
            error "footprint is enabled but TOBIAS BINDetect needs >= 2 groups with treat samples (found ${footprint_group_list.size()})"
        }
        def consensus_path = file(consensus_bed, checkIfExists: true)
        def genome = file(params.footprint.genome, checkIfExists: true)
        file("${genome}.fai", checkIfExists: true) // TOBIAS needs the indexed genome
        def motifs = file(params.footprint.motifs, checkIfExists: true)
        FOOTPRINT_EXT(Channel.fromList(footprint_group_list), consensus_path, genome, motifs)
    }

    emit:
    organelle_summary = ch_organelle_summary
    spikein_summary   = ch_spikein_summary
}
