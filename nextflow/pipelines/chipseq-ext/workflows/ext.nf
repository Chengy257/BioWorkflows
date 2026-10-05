/*
 * EXT - resolve the adopted nf-core/chipseq run's published outputs and wire
 * the default-off extension stages.
 *
 * The sample sheet is parsed at runtime (it is a project file, not pipeline
 * config): CSV with header sample_id,group,condition,control,batch. The
 * control column names the group's control sample (empty when the group has
 * none); controls may be their own rows (they are processed but excluded
 * from IDR treat pairs and from condition contrasts).
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

    // ---- Sample sheet (runtime) -------------------------------------------
    // Emits a single value: [ metaById: Map, groups: Map<group, List<sample>> ]
    ch_sheet = Channel
        .fromPath(params.ext_sample_sheet, checkIfExists: true)
        .splitCsv(header: true)
        .collect()
        .map { rows ->
            if (!rows) {
                error "sample sheet ${params.ext_sample_sheet} has no data rows"
            }
            def metaById = [:]
            def groups = [:]
            rows.each { row ->
                def id = (row.sample_id ?: '').trim()
                def group = (row.group ?: '').trim()
                if (!id || !group) {
                    error "sample sheet row missing sample_id or group: ${row}"
                }
                if (metaById.containsKey(id)) {
                    error "duplicate sample_id in sample sheet: ${id}"
                }
                def control = (row.control ?: '').trim()
                if (control && !groups.containsKey(group)) {
                    // remembered below after all rows are seen; nothing here
                }
                metaById[id] = [
                    id        : id,
                    group     : group,
                    condition : (row.condition ?: group).trim(),
                    control   : control,
                    batch     : (row.batch ?: '').trim(),
                    is_control: false,
                ]
                if (!groups.containsKey(group)) groups[group] = []
                groups[group] << id
            }
            // A sample named as a control must not simultaneously be a treat
            // inside the same group; mark controls and validate consistency.
            def groupControl = [:]
            groups.each { group, ids ->
                def controls = ids.collect { metaById[it].control }.findAll { it }
                def distinct = controls.unique()
                if (distinct.size() > 1) {
                    error "group ${group} declares more than one control: ${distinct}"
                }
                groupControl[group] = distinct ? distinct[0] : ''
            }
            // Control samples declared per group: the control row may live in
            // its own group (typical) or inside the treat group. Consistency
            // check: a control id must exist in the sheet.
            groupControl.each { group, ctl ->
                if (ctl && !metaById.containsKey(ctl)) {
                    error "group ${group} declares control '${ctl}' which is not a sample_id in the sheet"
                }
                if (ctl && groupControl.values().count(ctl) > 1) {
                    error "control '${ctl}' is declared by more than one group"
                }
            }
            [metaById, groups, groupControl]
        }

    // ---- Per-sample inputs -------------------------------------------------
    // BAMs: strict resolution - every sheet sample must have exactly one BAM
    // and every BAM must belong to a sheet sample.
    ch_meta_bams = ch_sheet
        .combine(Channel.fromPath("${bam_dir}/*${bam_suffix}", checkIfExists: true).collect())
        .flatMap { sheet, bamFiles ->
            def (metaById, groups, groupControl) = sheet
            def byId = [:]
            bamFiles.each { p ->
                def id = p.name.substring(0, p.name.length() - bam_suffix.length())
                if (!metaById.containsKey(id)) {
                    error "found BAM ${p.name} with no sample sheet row"
                }
                if (byId.containsKey(id)) {
                    error "duplicate BAM for sample ${id}"
                }
                byId[id] = p
            }
            metaById.keySet().each { id ->
                if (!byId.containsKey(id)) {
                    error "no BAM found for sheet sample ${id} in ${bam_dir}"
                }
            }
            metaById.collect { id, meta -> tuple(meta, byId[id]) }
        }

    ch_meta_peaks = ch_sheet
        .flatMap { sheet ->
            def (metaById, groups, groupControl) = sheet
            metaById.collect { id, meta ->
                tuple(meta, file("${peaks_dir}/${id}${peak_suffix}", checkIfExists: true))
            }
        }

    // Consensus peaks per antibody (optional: single-sample-no-control runs
    // publish no consensus; stages consuming them simply do not run).
    ch_consensus = Channel
        .fromPath("${consensus_dir}/*/*.consensus_peaks.bed")
        .map { p -> tuple([id: p.parent.name], p) }

    // ---- Optional QC inputs (empty when the glob matches nothing) ----------
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
        ch_idr_pairs = ch_sheet
            .flatMap { sheet ->
                def (metaById, groups, groupControl) = sheet
                def out = []
                groups.each { group, ids ->
                    def treats = ids.findAll { it != groupControl[group] }.sort()
                    if (treats.size() < 2) return // continue
                    (0..<treats.size()).each { i ->
                        ((i + 1)..<treats.size()).each { j ->
                            def a = treats[i]
                            def b = treats[j]
                            out << tuple(
                                [id: "${group}__${a}_vs_${b}", group: group, pair: [a, b]],
                                file("${peaks_dir}/${a}${peak_suffix}", checkIfExists: true),
                                file("${peaks_dir}/${b}${peak_suffix}", checkIfExists: true),
                            )
                        }
                    }
                }
                if (!out) {
                    error "idr is enabled but no group has >= ${params.idr.min_replicates} treat replicates"
                }
                out
            }
        IDR_REPRODUCIBILITY(ch_idr_pairs)
    }

    // ---- Stage: pooled per-group SEACR calling -------------------------------
    if (params.seacr.enabled) {
        ch_seacr_groups = ch_sheet
            .flatMap { sheet ->
                def (metaById, groups, groupControl) = sheet
                def out = []
                groups.each { group, ids ->
                    def treats = ids.findAll { it != groupControl[group] }.sort()
                    if (!treats) return // continue
                    def ctl = groupControl[group]
                    out << tuple(
                        [id: group, group: group, control: ctl],
                        treats.collect { id -> metaById[id] },
                    )
                }
                if (!out) {
                    error "seacr is enabled but the sample sheet has no groups"
                }
                out
            }
        SEACR_CALLING(ch_seacr_groups, ch_meta_bams)
    }

    // ---- Stage: HOMER motif discovery on consensus peaks ---------------------
    if (params.homer_motif.enabled) {
        HOMER_MOTIF_EXT(ch_consensus)
    }

    // ---- Stage: DiffBind differential binding --------------------------------
    if (params.diffbind.enabled) {
        DIFFBIND_EXT(ch_sheet, ch_meta_bams, ch_meta_peaks)
    }

    emit:
    organelle_summary = ch_organelle_summary
}
