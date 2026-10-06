/*
 * EXT - resolve the adopted nf-core/chipseq run's published outputs and wire
 * the default-off extension stages.
 *
 * The sample sheet is parsed synchronously in the workflow body (it is a
 * project file, not pipeline config): CSV with header
 * sample_id,group,condition,control,batch. The control column names the
 * group's control sample (empty when the group has none); one control may
 * serve several treat groups (shared IgG, legacy semantics). Controls may be
 * their own rows - they are processed but excluded from IDR treat pairs and
 * from DiffBind condition rows.
 */

include { ORGANELLE_QC_EXT } from '../subworkflows/local/organelle_qc_ext'
include { QC_GATES_EXT } from '../subworkflows/local/qc_gates_ext'
include { SPIKEIN_SCALING } from '../subworkflows/local/spikein_scaling'
include { IDR_REPRODUCIBILITY } from '../subworkflows/local/idr_reproducibility'
include { SEACR_CALLING } from '../subworkflows/local/seacr_calling'
include { HOMER_MOTIF_EXT } from '../subworkflows/local/homer_motif_ext'
include { DIFFBIND_EXT } from '../subworkflows/local/diffbind_ext'

def require_params() {
    if (!params.ext_results_dir) {
        error "params.ext_results_dir is required: point it at the adopted nf-core/chipseq run's --outdir"
    }
    if (!params.ext_sample_sheet) {
        error "params.ext_sample_sheet is required: CSV with sample_id,group,condition,control,batch"
    }
    if (params.diffbind.enabled) {
        if (!params.diffbind.contrasts) {
            error "diffbind is enabled but params.diffbind.contrasts is empty"
        }
        params.diffbind.contrasts.each { c ->
            if (!(c =~ /.+_vs_.+/)) {
                error "invalid diffbind contrast '${c}': expected '<conditionA>_vs_<conditionB>'"
            }
        }
    }
    if (params.homer_motif.enabled && !params.homer_motif.genome) {
        error "homer_motif is enabled but params.homer_motif.genome is empty (HOMER genome tag or custom FASTA path)"
    }
    if (params.seacr.enabled && !['stringent', 'relaxed'].contains(params.seacr.mode)) {
        error "params.seacr.mode must be 'stringent' or 'relaxed', got '${params.seacr.mode}'"
    }
    if (params.seacr.enabled && !['norm', 'non'].contains(params.seacr.normalize)) {
        error "params.seacr.normalize must be 'norm' or 'non', got '${params.seacr.normalize}'"
    }
}

workflow EXT {

    main:
    require_params()

    def results = params.ext_results_dir.replaceAll('/+$', '')
    def aligner = params.ext_aligner
    def peak_ext = params.ext_peak_type == 'narrow_peak' ? 'narrowPeak' : 'broadPeak'

    def bam_dir = params.ext_bam_dir ?: "${results}/${aligner}/filtered_bam"
    def peaks_dir = params.ext_peaks_dir ?: "${results}/${aligner}/merged_library/macs3/${params.ext_peak_type}"
    def consensus_dir = params.ext_consensus_dir ?: "${peaks_dir}/consensus"
    def bigwig_dir = params.ext_bigwig_dir ?: "${results}/${aligner}/bigwig"
    def flagstat_dir = params.ext_flagstat_dir ?: "${results}/${aligner}/flagstat"

    def bam_suffix = params.ext_bam_suffix
    def bigwig_suffix = params.ext_bigwig_suffix
    def flagstat_suffix = params.ext_flagstat_suffix
    def peak_suffix = "_peaks.${peak_ext}"

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
        if (id.contains('_vs_')) {
            error "sample id '${id}' contains '_vs_' (reserved for the IDR pair slug); rename it"
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
    ch_dup_metrics = params.ext_dup_glob
        ? Channel.fromPath(params.ext_dup_glob)
            .map { p -> tuple(p.name.endsWith(params.ext_dup_suffix)
                ? p.name.substring(0, p.name.length() - params.ext_dup_suffix.length()) : p.simpleName, p) }
        : channel.empty()
    ch_bigwigs = Channel
        .fromPath("${bigwig_dir}/*${bigwig_suffix}")
        .map { p -> tuple(p.name.substring(0, p.name.length() - bigwig_suffix.length()), p) }

    // ---- Stage: organelle QC ------------------------------------------------
    ch_organelle_summary = channel.empty()
    if (params.organelle_qc.enabled) {
        ORGANELLE_QC_EXT(ch_meta_bams)
        ch_organelle_summary = ORGANELLE_QC_EXT.out.summary
    }

    // ---- Stage: QC gate table ------------------------------------------------
    if (params.qc_gates.enabled) {
        QC_GATES_EXT(
            ch_flagstats,
            ch_dup_metrics,
            ch_organelle_summary.toList(),
        )
    }

    // ---- Stage: spike-in scaling (combined-reference runs) -------------------
    if (params.spike_in.enabled) {
        SPIKEIN_SCALING(ch_meta_bams)
    }

    // ---- Stage: narrow-peak IDR on native per-sample calls -------------------
    if (params.idr.enabled) {
        def idr_pair_list = []
        groups.each { group, ids ->
            def treats = ids.findAll { it != groupControl[group] }.sort()
            if (treats.size() < params.idr.min_replicates) return // continue
            (0..<treats.size()).each { i ->
                ((i + 1)..<treats.size()).each { j ->
                    def a = treats[i]
                    def b = treats[j]
                    idr_pair_list << tuple(
                        [id: "${group}__${a}_vs_${b}", group: group],
                        file("${peaks_dir}/${a}${peak_suffix}", checkIfExists: true),
                        file("${peaks_dir}/${b}${peak_suffix}", checkIfExists: true),
                    )
                }
            }
        }
        if (!idr_pair_list) {
            error "idr is enabled but no group has >= ${params.idr.min_replicates} treat replicates"
        }
        IDR_REPRODUCIBILITY(Channel.fromList(idr_pair_list))
    }

    // ---- Stage: pooled per-group SEACR calling -------------------------------
    if (params.seacr.enabled) {
        def seacr_group_list = []
        groups.each { group, ids ->
            def treats = ids.findAll { it != groupControl[group] }.sort()
            if (!treats) return // continue
            // one control sample may serve several treat groups (shared IgG);
            // resolve its BAM synchronously so every group gets its own copy
            def ctl_id = groupControl[group]
            def ctl_bam = ctl_id
                ? file("${bam_dir}/${ctl_id}${bam_suffix}", checkIfExists: true)
                : null
            seacr_group_list << tuple(
                [id: group, group: group, control: ctl_id],
                treats.collect { id -> metaById[id] },
                ctl_bam,
            )
        }
        if (!seacr_group_list) {
            error "seacr is enabled but the sample sheet has no groups"
        }
        SEACR_CALLING(Channel.fromList(seacr_group_list), ch_meta_bams)
    }

    // ---- Stage: HOMER motif discovery on consensus peaks ---------------------
    if (params.homer_motif.enabled) {
        ch_consensus = Channel
            .fromPath("${consensus_dir}/*/*.consensus_peaks.bed")
            .map { p -> tuple([id: p.parent.name], p) }
        HOMER_MOTIF_EXT(ch_consensus)
    }

    // ---- Stage: DiffBind differential binding --------------------------------
    if (params.diffbind.enabled) {
        ch_meta_peaks = Channel
            .fromPath("${peaks_dir}/*${peak_suffix}", checkIfExists: true)
            .map { p ->
                def id = p.name.substring(0, p.name.length() - peak_suffix.length())
                if (!metaById.containsKey(id)) {
                    error "found peak file ${p.name} with no sample sheet row"
                }
                tuple(metaById[id], p)
            }
        DIFFBIND_EXT(ch_meta_bams, ch_meta_peaks)
    }

    emit:
    organelle_summary = ch_organelle_summary
}
