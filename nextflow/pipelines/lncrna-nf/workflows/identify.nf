/*
 * IDENTIFY - the lncRNA identification chain.
 *
 * Input resolution mirrors the Stage-2 handoff contract: published
 * {aligner}-route BAMs under params.bam_dir (markdup BAMs by default;
 * bam_suffix selects the pre-dedup STAR BAMs when the adopted run passed
 * --save_align_intermeds). No alignment happens here.
 *
 * Stage map (decision record Section 5 dispositions):
 *   STRINGTIE_STRINGTIE     per-sample --conservative assembly   (kept)
 *   STRINGTIE_MERGE         merge vs gtf_pcg                     (kept)
 *   GFFCOMPARE              classcode assignment vs gtf_pcg      (kept)
 *   NOVEL_TRANSCRIPTS       classcode-u extraction               (kept, local)
 *   GFFREAD                 novel FASTA                          (kept)
 *   CPC2 / CPAT_PREDICT     coding potential                     (kept / replaced CNCI)
 *   CODING_POTENTIAL_SET    intersection + length >= 200         (kept, local)
 *   GFFREAD (2nd)           noncoding FASTA for the filters
 *   PFAM_SCAN               protein-domain filter                (kept, local module)
 *   DIAMOND_BLASTX          NR homology filter                   (kept)
 *   FINAL_LNCRNA_SET        union-filter + ref_withLncRNA.gtf    (kept, local)
 *   SUBREAD_FEATURECOUNTS   expression vs ref_withLncRNA         (kept)
 *   MERGE_LNCRNA_COUNTS     count/TPM matrices                   (kept, local)
 */

include { STRINGTIE_STRINGTIE } from '../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_MERGE } from '../modules/nf-core/stringtie/merge/main'
include { GFFCOMPARE } from '../modules/nf-core/gffcompare/main'
include { GFFREAD as GFFREAD_NOVEL } from '../modules/nf-core/gffread/main'
include { GFFREAD as GFFREAD_NONCODING } from '../modules/nf-core/gffread/main'
include { CPC2 } from '../modules/nf-core/cpc2/main'
include { DIAMOND_BLASTX } from '../modules/nf-core/diamond/blastx/main'
include { SUBREAD_FEATURECOUNTS } from '../modules/nf-core/subread/featurecounts/main'
include { NOVEL_TRANSCRIPTS } from '../modules/local/novel_transcripts/main'
include { CPAT_PREDICT } from '../modules/local/cpat_predict/main'
include { CODING_POTENTIAL_SET } from '../modules/local/coding_potential_set/main'
include { PFAM_SCAN } from '../modules/local/pfam_scan/main'
include { FINAL_LNCRNA_SET } from '../modules/local/final_lncrna_set/main'
include { MERGE_LNCRNA_COUNTS } from '../modules/local/merge_lncrna_counts/main'

def require_params() {
    if (!params.bam_dir) {
        error "params.bam_dir is required: the adopted nf-core/rnaseq run's {outdir}/{aligner} directory"
    }
    if (!params.sample_sheet) {
        error "params.sample_sheet is required: CSV with a sample_id column (the Stage-2 samplesheet works as-is)"
    }
    if (!params.genome) {
        error "params.genome is required: reference genome FASTA (a sibling .fai must exist)"
    }
    if (!params.gtf) {
        error "params.gtf is required: reference annotation GTF"
    }
    if (!params.gtf_pcg) {
        error "params.gtf_pcg is required: protein-coding CDS GTF (rice: the MSU+RAP combined CDS GTF - see README for the build recipe)"
    }
    if (params.strand && !(params.strand in ['rf', 'fr'])) {
        error "params.strand must be 'rf', 'fr' or '' (got '${params.strand}')"
    }
    if (params.pfam.enabled && !params.pfam.db_dir) {
        error "pfam is enabled but params.pfam.db_dir is empty (Pfam database directory for pfam_scan.pl)"
    }
    if (params.diamond.enabled && !params.diamond.db) {
        error "diamond is enabled but params.diamond.db is empty (path to nr.dmnd)"
    }
    if (params.cpat.enabled && (!params.cpat.model || !params.cpat.hexamer)) {
        error "cpat is enabled but params.cpat.model/params.cpat.hexamer are empty (species logit model + hexamer table; human files ship with CPAT, a rice model must be user-provided - see README)"
    }
}

workflow IDENTIFY {

    main:
    require_params()

    // ---- Sample sheet: sample_id column (legacy single-column rule; the
    // Stage-2 samplesheet's fastq_2 column, when present, derives PE/SE) ----
    def sheet_lines = file(params.sample_sheet, checkIfExists: true)
        .readLines()
        .findAll { it.trim() }
    if (sheet_lines.size() < 2) {
        error "sample sheet ${params.sample_sheet} has no data rows"
    }
    def header = sheet_lines[0].split(',').collect { it.trim() }
    if (!header.contains('sample_id')) {
        error "sample sheet is missing the required column 'sample_id' (header: ${header})"
    }
    def id_col = header.indexOf('sample_id')
    def fastq2_col = header.indexOf('fastq_2')
    def sample_ids = []
    def paired_seen = false
    sheet_lines.tail().each { line ->
        def cells = line.split(',')
        def id = (id_col < cells.size() ? cells[id_col] : '').trim()
        if (!id) {
            error "sample sheet row missing sample_id: ${line}"
        }
        if (sample_ids.contains(id)) {
            error "duplicate sample_id in sample sheet: ${id}"
        }
        sample_ids << id
        if (fastq2_col >= 0 && fastq2_col < cells.size() && cells[fastq2_col].trim()) {
            paired_seen = true
        }
    }

    // PE/SE layout is project-wide (the Stage-2 sheet is homogeneous): a
    // fastq_2 column decides it; sheets without one fall back to the param.
    def single_end = fastq2_col >= 0 ? !paired_seen : params.single_end

    // Strandedness: params.strand carries the StringTie spelling (rf/fr);
    // the featureCounts module consumes the Salmon/legacy spelling.
    // rf = fr-firststrand -> -s 2 (reverse); fr = fr-secondstrand -> -s 1.
    def strandedness = params.strand == 'rf' ? 'reverse'
        : (params.strand == 'fr' ? 'forward' : 'unstranded')

    // ---- BAM resolution (strict: one BAM per sheet sample, no strays) -------
    def bam_dir_path = file(params.bam_dir)
    if (!bam_dir_path.exists()) {
        error "BAM directory not found: ${params.bam_dir}"
    }
    def bam_suffix = params.bam_suffix
    def bam_ids = []
    bam_dir_path.list().findAll { it.endsWith(bam_suffix) }.sort().each { name ->
        def id = name.substring(0, name.length() - bam_suffix.length())
        if (!sample_ids.contains(id)) {
            error "found BAM ${name} with no sample sheet row"
        }
        bam_ids << id
    }
    def missing = sample_ids.findAll { !(it in bam_ids) }
    if (missing) {
        error "no BAM found for sheet sample(s) in ${params.bam_dir}: ${missing.sort().join(', ')}"
    }

    def genome = file(params.genome, checkIfExists: true)
    file("${genome}.fai", checkIfExists: true)
    def gtf = file(params.gtf, checkIfExists: true)
    def gtf_pcg = file(params.gtf_pcg, checkIfExists: true)

    ch_meta_bams = Channel.fromList(bam_ids.collect { id ->
        tuple([id: id, single_end: single_end, strandedness: strandedness],
              file("${params.bam_dir}/${id}${bam_suffix}", checkIfExists: true))
    })

    // ---- Assembly + merge + novel extraction (core chain) -------------------
    STRINGTIE_STRINGTIE(
        ch_meta_bams.map { meta, bam -> tuple(meta, bam, []) },
        '',                        // mode: plain assembly (legacy --conservative route)
        gtf,
    )
    STRINGTIE_MERGE(
        // .toList() emits ONE item (a list of (meta, path) tuples); reshape
        // into the (meta, [paths]) the module input declares - handing the
        // raw list of tuples to a path input fails task binding.
        STRINGTIE_STRINGTIE.out.transcript_gtf.toList()
            .map { rows -> tuple([id: 'lncrna'], rows.collect { row -> row[1] }) },
        tuple([id: 'lncrna'], gtf_pcg),
    )
    GFFCOMPARE(
        STRINGTIE_MERGE.out.merged_gtf.map { meta, g -> tuple(meta, [g]) },
        tuple([id: 'lncrna'], genome, file("${genome}.fai", checkIfExists: true)),
        tuple([id: 'lncrna'], gtf_pcg),
    )
    NOVEL_TRANSCRIPTS(
        GFFCOMPARE.out.annotated_gtf,
    )
    GFFREAD_NOVEL(
        NOVEL_TRANSCRIPTS.out.novel_gtf.map { meta, g -> tuple(meta, g) },
        genome,
    )

    // ---- Coding potential: CPC2 AND CPAT AND length >= min_length -----------
    CPC2(GFFREAD_NOVEL.out.gffread_fasta)
    CPAT_PREDICT(
        GFFREAD_NOVEL.out.gffread_fasta,
        file(params.cpat.model, checkIfExists: true),
        file(params.cpat.hexamer, checkIfExists: true),
    )
    CODING_POTENTIAL_SET(
        CPC2.out.txt,
        CPAT_PREDICT.out.prob_best,
        CPAT_PREDICT.out.no_orf,
        NOVEL_TRANSCRIPTS.out.novel_gtf,
        params.min_length,
    )
    GFFREAD_NONCODING(
        CODING_POTENTIAL_SET.out.noncoding_gtf.map { meta, g -> tuple(meta, g) },
        genome,
    )

    // ---- Homology filters (default ON, per-project OFF switches) ------------
    // An OFF switch must not starve FINAL_LNCRNA_SET: the placeholder is an
    // empty-list VALUE channel (the documented replacement for the removed
    // input optional:true), not Channel.empty().
    ch_pfam_hits = Channel.value([[id: 'lncrna'], []])
    if (params.pfam.enabled) {
        PFAM_SCAN(
            CODING_POTENTIAL_SET.out.noncoding_gtf.map { meta, g -> tuple(meta, g) },
            GFFREAD_NONCODING.out.gffread_fasta,
            file(params.pfam.db_dir, checkIfExists: true),
            params.pfam.evalue,
        )
        ch_pfam_hits = PFAM_SCAN.out.hit_ids
    }
    ch_diamond_hits = Channel.value([[id: 'lncrna'], []])
    if (params.diamond.enabled) {
        DIAMOND_BLASTX(
            // the query is the NONCODING FASTA (GFFREAD's second pass), not
            // the GTF - diamond blastx reads sequences
            GFFREAD_NONCODING.out.gffread_fasta,
            tuple([id: 'lncrna'], file(params.diamond.db, checkIfExists: true)),
            'txt',                   // outfmt 6 tabular: column 1 = query id
            'qseqid sseqid',
        )
        ch_diamond_hits = DIAMOND_BLASTX.out.txt
    }

    // ---- Final set + expression ---------------------------------------------
    FINAL_LNCRNA_SET(
        ch_pfam_hits,
        ch_diamond_hits,
        CODING_POTENTIAL_SET.out.noncoding_gtf,
        gtf,
        genome,
    )
    // The annotation is one value channel for every per-sample task (reading
    // a channel inside a map closure does not compile in Nextflow).
    ch_ref_gtf = FINAL_LNCRNA_SET.out.ref_with_lncrna_gtf.first().map { _meta, g -> g }
    SUBREAD_FEATURECOUNTS(
        ch_meta_bams.combine(ch_ref_gtf).map { meta, bam, ref -> tuple(meta, [bam], ref) },
    )
    MERGE_LNCRNA_COUNTS(
        // .toList() emits ONE item (a list of (meta, path) tuples); reshape
        // into the (meta, [paths]) the module input declares - and stage the
        // .summary files too (the helper reads Assigned reads from them).
        SUBREAD_FEATURECOUNTS.out.counts.toList()
            .map { rows -> tuple([id: 'lncrna'], rows.collect { row -> row[1] }) },
        SUBREAD_FEATURECOUNTS.out.summary.toList()
            .map { rows -> tuple([id: 'lncrna'], rows.collect { row -> row[1] }) },
        sample_ids.sort().join(','),
    )

    emit:
    novel_gtf    = NOVEL_TRANSCRIPTS.out.novel_gtf
    noncoding    = CODING_POTENTIAL_SET.out.noncoding_ids
    final_set    = FINAL_LNCRNA_SET.out.final_fasta
    final_gtf    = FINAL_LNCRNA_SET.out.final_gtf
    counts       = MERGE_LNCRNA_COUNTS.out.counts
    tpm          = MERGE_LNCRNA_COUNTS.out.tpm
}
