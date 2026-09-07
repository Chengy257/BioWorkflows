# Cross-sample reproducible peaks (v0.2, optional) and GTF-based peak
# annotation (v0.2, optional). Both stages are default-off and parse-time
# gated in the same style as the callpeak_clipper block:
#
#   reproducible_peaks.enabled  -> consensus_peaks groups the ip-role
#       PureCLIP beds of one condition into a bedtools multiinter consensus
#       (a site is kept when present in >= min_replicates beds; column 4 of
#       the consensus BED reports the support).
#   annotate_peaks.enabled -> gtf_gene_regions derives gene/exon BEDs plus a
#       gene attribute table from the already-required GTF (stdlib parser);
#       the annotate rules then classify every peak of a set via
#       bedtools intersect (exon / gene overlap) and bedtools closest
#       (nearest gene + signed distance) and merge everything into one TSV
#       per peak set with the stdlib script annotate_peaks.py.
#
# Depends on workflow/rules/common.smk for SAMPLES / SAMPLE_CONDITIONS /
# SAMPLE_ROLES, CONDITIONS, REPRODUCIBLE_PEAKS_ENABLED / ANNOTATE_PEAKS_ENABLED,
# MIN_REPLICATES, condition_ip_samples(), R(), and the resource helpers.
# Helper functions live in common.smk to keep this module rules-only.

_SCRIPTS_DIR = os.path.join(WORKFLOW_DIR, "scripts")

if ANNOTATE_PEAKS_ENABLED:
    rule gtf_gene_regions:
        input:
            gtf=config["gtf"],
        output:
            genes=R("6.annotation/_ref/genes.bed"),
            exons=R("6.annotation/_ref/exons.bed"),
            table=R("6.annotation/_ref/genes.tsv"),
        params:
            script=os.path.join(_SCRIPTS_DIR, "gtf_to_gene_regions.py"),
            python=os.environ.get("SECLIP_PYTHON", "python3"),
        log:
            R("logs/gtf_gene_regions.log"),
        threads: rthreads("gtf_gene_regions")
        resources:
            mem_mb=rmem("gtf_gene_regions"),
            runtime_min=rruntime("gtf_gene_regions"),
            runtime_sec=rruntime_sec("gtf_gene_regions"),
        shell:
            "{params.python} {params.script} --gtf {input.gtf} "
            "--genes-out {output.genes} --exons-out {output.exons} "
            "--table-out {output.table} >> {log} 2>&1"

    rule annotate_sample_peaks:
        input:
            peaks=R("5.callpeak/{sample}.pureclip.bed"),
            genes=R("6.annotation/_ref/genes.bed"),
            exons=R("6.annotation/_ref/exons.bed"),
            table=R("6.annotation/_ref/genes.tsv"),
        output:
            R("6.annotation/{sample}.annotation.tsv"),
        params:
            script=os.path.join(_SCRIPTS_DIR, "annotate_peaks.py"),
            python=os.environ.get("SECLIP_PYTHON", "python3"),
            # PureCLIP writes BED6 (chrom start end name score strand): the
            # peak occupies 6 columns and the crosslink-site score is column 5.
            peak_cols=6,
            score_col=5,
        log:
            R("logs/annotate_sample_peaks/{sample}.log"),
        threads: rthreads("annotate_peaks")
        resources:
            mem_mb=rmem("annotate_peaks"),
            runtime_min=rruntime("annotate_peaks"),
            runtime_sec=rruntime_sec("annotate_peaks"),
        shell:
            """
            ## Classify each peak against the GTF-derived features and find
            ## its nearest gene, then merge into the final TSV (stdlib
            ## script): intersect -u gives the exon/gene overlaps,
            ## closest -d -t first the nearest gene and the signed distance
            ## (0 = overlapping).
            tmpdir={resources.tmpdir}/annotate_{wildcards.sample}
            rm -rf "$tmpdir" && mkdir -p "$tmpdir"
            bedtools intersect -a {input.peaks} -b {input.exons} -u > "$tmpdir/exons.u.bed" 2>> {log}
            bedtools intersect -a {input.peaks} -b {input.genes} -u > "$tmpdir/genes.u.bed" 2>> {log}
            bedtools closest -a {input.peaks} -b {input.genes} -d -t first > "$tmpdir/closest.tsv" 2>> {log}
            {params.python} {params.script} --peaks {input.peaks} \\
                --exon-hits "$tmpdir/exons.u.bed" --gene-hits "$tmpdir/genes.u.bed" \\
                --closest "$tmpdir/closest.tsv" --gene-table {input.table} \\
                --peak-cols {params.peak_cols} --score-col {params.score_col} \\
                --out {output} >> {log} 2>&1
            rm -rf "$tmpdir"
            """

    if REPRODUCIBLE_PEAKS_ENABLED:
        rule annotate_consensus_peaks:
            input:
                peaks=R("6.reproducible_peaks/{condition}.consensus.bed"),
                genes=R("6.annotation/_ref/genes.bed"),
                exons=R("6.annotation/_ref/exons.bed"),
                table=R("6.annotation/_ref/genes.tsv"),
            output:
                R("6.annotation/{condition}.consensus.annotation.tsv"),
            params:
                script=os.path.join(_SCRIPTS_DIR, "annotate_peaks.py"),
                python=os.environ.get("SECLIP_PYTHON", "python3"),
                # Consensus BEDs are 4-column (chrom start end support): the
                # support doubles as the annotation score column.
                peak_cols=4,
                score_col=4,
            log:
                R("logs/annotate_consensus_peaks/{condition}.log"),
            threads: rthreads("annotate_peaks")
            resources:
                mem_mb=rmem("annotate_peaks"),
                runtime_min=rruntime("annotate_peaks"),
                runtime_sec=rruntime_sec("annotate_peaks"),
            shell:
                """
                ## Same pipeline as the per-sample annotation, applied to the
                ## condition consensus BED (4 columns; support as score).
                tmpdir={resources.tmpdir}/annotate_consensus_{wildcards.condition}
                rm -rf "$tmpdir" && mkdir -p "$tmpdir"
                bedtools intersect -a {input.peaks} -b {input.exons} -u > "$tmpdir/exons.u.bed" 2>> {log}
                bedtools intersect -a {input.peaks} -b {input.genes} -u > "$tmpdir/genes.u.bed" 2>> {log}
                bedtools closest -a {input.peaks} -b {input.genes} -d -t first > "$tmpdir/closest.tsv" 2>> {log}
                {params.python} {params.script} --peaks {input.peaks} \\
                    --exon-hits "$tmpdir/exons.u.bed" --gene-hits "$tmpdir/genes.u.bed" \\
                    --closest "$tmpdir/closest.tsv" --gene-table {input.table} \\
                    --peak-cols {params.peak_cols} --score-col {params.score_col} \\
                    --out {output} >> {log} 2>&1
                rm -rf "$tmpdir"
                """


if REPRODUCIBLE_PEAKS_ENABLED:
    rule consensus_peaks:
        input:
            beds=lambda wc: [R(f"5.callpeak/{s}.pureclip.bed")
                             for s in condition_ip_samples(wc.condition)],
        output:
            R("6.reproducible_peaks/{condition}.consensus.bed"),
        params:
            min_replicates=MIN_REPLICATES,
        log:
            R("logs/consensus_peaks/{condition}.log"),
        threads: rthreads("consensus_peaks")
        resources:
            mem_mb=rmem("consensus_peaks"),
            runtime_min=rruntime("consensus_peaks"),
            runtime_sec=rruntime_sec("consensus_peaks"),
        shell:
            """
            ## multiinter needs coordinate-sorted inputs; sort each per-sample
            ## bed into the per-job temp dir (see the STAR rules: use the
            ## scheduler-provided tmpdir, not the working directory).
            tmpdir={resources.tmpdir}/consensus_{wildcards.condition}
            rm -rf "$tmpdir" && mkdir -p "$tmpdir"
            for bed in {input.beds}; do
                sort -k1,1 -k2,2 "$bed" > "$tmpdir/$(basename "$bed").sorted"
            done
            ## multiinter columns: chrom, start, end, nclust (number of input
            ## files carrying the feature), then per-file membership lists.
            ## Keep sites present in >= min_replicates beds; the consensus BED
            ## reports chrom, start, end, support.
            bedtools multiinter -i "$tmpdir"/*.sorted 2>> {log} \\
                | awk '$4 >= {params.min_replicates}' | cut -f1-4 \\
                > {output} 2>> {log}
            rm -rf "$tmpdir"
            """
