# Cross-sample reproducible peaks (v0.2, optional) and GTF-based peak
# annotation (v0.2, optional). Both stages are default-off and parse-time
# gated in the same style as the callpeak_clipper block:
#
#   reproducible_peaks.enabled  -> consensus_peaks groups the ip-role
#       PureCLIP beds of one condition into a bedtools multiinter consensus
#       (a site is kept when present in >= min_replicates beds; column 4 of
#       the consensus BED reports the support).
#   reproducible_peaks.input_control (requires enabled) -> the role=input
#       samples of a condition also run the peak-calling chain (their target
#       scoping lives in common.smk); input_background unions their PureCLIP
#       beds into a per-condition background, flag_input_background appends
#       the binary in_input_background column (0/1) to the final consensus
#       BED (BED5), and filter_input_background (filter_by_input=true) drops
#       the flagged sites into {condition}.consensus.filtered.bed (BED4).
#       Conditions without input samples skip the whole branch and keep the
#       plain W7 consensus pipeline.
#   annotate_peaks.enabled -> gtf_gene_regions derives gene/exon BEDs plus a
#       gene attribute table from the already-required GTF (stdlib parser);
#       the annotate rules then classify every peak of a set via
#       bedtools intersect (exon / gene overlap) and bedtools closest
#       (nearest gene + signed distance) and merge everything into one TSV
#       per peak set with the stdlib script annotate_peaks.py.
#
# Depends on workflow/rules/common.smk for SAMPLES / SAMPLE_CONDITIONS /
# SAMPLE_ROLES, CONDITIONS, CONDITIONS_WITH_INPUT / CONDITIONS_WITHOUT_INPUT,
# REPRODUCIBLE_PEAKS_ENABLED / ANNOTATE_PEAKS_ENABLED, INPUT_CONTROL_ENABLED /
# FILTER_BY_INPUT, MIN_REPLICATES, condition_ip_samples() /
# condition_input_samples(), consensus_annotated_bed() / consensus_peak_cols(),
# R(), and the resource helpers. Helper functions live in common.smk to keep
# this module rules-only.

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
                peaks=lambda wc: consensus_annotated_bed(wc.condition),
                genes=R("6.annotation/_ref/genes.bed"),
                exons=R("6.annotation/_ref/exons.bed"),
                table=R("6.annotation/_ref/genes.tsv"),
            output:
                R("6.annotation/{condition}.consensus.annotation.tsv"),
            params:
                script=os.path.join(_SCRIPTS_DIR, "annotate_peaks.py"),
                python=os.environ.get("SECLIP_PYTHON", "python3"),
                # Consensus BEDs are 4-column (chrom start end support; the
                # support doubles as the annotation score column) or, with
                # input_control and no filtering, 5-column (binary
                # in_input_background flag appended as column 5, which the
                # annotation drops per its column contract).
                peak_cols=lambda wc: consensus_peak_cols(wc.condition),
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
                ## condition consensus BED (4 columns, or 5 with the input
                ## flag; support as score in both shapes).
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
    if INPUT_CONTROL_ENABLED:
        # With input_control the raw (pre-flagging) consensus of the ip
        # samples is an intermediate; the final {condition}.consensus.bed is
        # the flagged BED5 for conditions with input samples and the plain
        # BED4 for conditions without. The two producing rules of the final
        # pattern are disambiguated by per-rule wildcard constraints.
        _INPUT_CONDITIONS_RE = (
            "(?:" + "|".join(re.escape(c) for c in CONDITIONS_WITH_INPUT) + ")"
        ) if CONDITIONS_WITH_INPUT else "(?!x)x"
        _NO_INPUT_CONDITIONS_RE = (
            "(?:" + "|".join(re.escape(c) for c in CONDITIONS_WITHOUT_INPUT) + ")"
        ) if CONDITIONS_WITHOUT_INPUT else "(?!x)x"

        if CONDITIONS_WITH_INPUT:
            rule consensus_peaks_raw:
                # Identical multiinter pipeline to the plain consensus, kept
                # as an intermediate; flag_input_background turns it into the
                # final consensus BED5.
                input:
                    beds=lambda wc: [R(f"5.callpeak/{s}.pureclip.bed")
                                     for s in condition_ip_samples(wc.condition)],
                output:
                    R("6.reproducible_peaks/{condition}.consensus.raw.bed"),
                params:
                    min_replicates=MIN_REPLICATES,
                log:
                    R("logs/consensus_peaks_raw/{condition}.log"),
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
                    tmpdir={resources.tmpdir}/consensus_raw_{wildcards.condition}
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

            rule input_background:
                # Union consensus of the condition's input-control PureCLIP
                # beds: every multiinter row carries support >= 1 by
                # construction, so the whole union is kept (column 4 reports
                # how many input beds cover the feature).
                input:
                    beds=lambda wc: [R(f"5.callpeak/{s}.pureclip.bed")
                                     for s in condition_input_samples(wc.condition)],
                output:
                    R("6.reproducible_peaks/{condition}.input_background.bed"),
                log:
                    R("logs/input_background/{condition}.log"),
                threads: rthreads("input_background")
                resources:
                    mem_mb=rmem("input_background"),
                    runtime_min=rruntime("input_background"),
                    runtime_sec=rruntime_sec("input_background"),
                shell:
                    """
                    ## multiinter needs coordinate-sorted inputs; sort each
                    ## per-sample bed into the per-job temp dir (see the STAR
                    ## rules: use the scheduler-provided tmpdir, not the
                    ## working directory).
                    tmpdir={resources.tmpdir}/input_background_{wildcards.condition}
                    rm -rf "$tmpdir" && mkdir -p "$tmpdir"
                    for bed in {input.beds}; do
                        sort -k1,1 -k2,2 "$bed" > "$tmpdir/$(basename "$bed").sorted"
                    done
                    ## All multiinter rows have support >= 1, so the union
                    ## keeps every row (chrom, start, end, support).
                    bedtools multiinter -i "$tmpdir"/*.sorted 2>> {log} \\
                        | cut -f1-4 > {output} 2>> {log}
                    rm -rf "$tmpdir"
                    """

            rule flag_input_background:
                input:
                    raw=R("6.reproducible_peaks/{condition}.consensus.raw.bed"),
                    background=R("6.reproducible_peaks/{condition}.input_background.bed"),
                output:
                    R("6.reproducible_peaks/{condition}.consensus.bed"),
                wildcard_constraints:
                    condition=_INPUT_CONDITIONS_RE,
                log:
                    R("logs/flag_input_background/{condition}.log"),
                threads: rthreads("flag_input_background")
                resources:
                    mem_mb=rmem("flag_input_background"),
                    runtime_min=rruntime("flag_input_background"),
                    runtime_sec=rruntime_sec("flag_input_background"),
                shell:
                    """
                    ## Flag every consensus site overlapping the condition's
                    ## input background: intersect -c appends the number of
                    ## overlapping background intervals, awk collapses it to
                    ## the binary in_input_background column (1 = in
                    ## background). The final consensus BED keeps the W7
                    ## columns 1-4 and appends the flag as column 5.
                    bedtools intersect -a {input.raw} -b {input.background} -c 2>> {log} \\
                        | awk 'BEGIN {{OFS = "\\t"}} {{$5 = ($5 > 0) ? 1 : 0; print}}' \\
                        > {output} 2>> {log}
                    """

            if FILTER_BY_INPUT:
                rule filter_input_background:
                    # Drop consensus sites flagged as input background (keep
                    # in_input_background == 0); the filtered BED returns to
                    # the plain W7 BED4 shape (support in column 4) and
                    # replaces the flagged consensus in the annotation stage.
                    input:
                        flagged=R("6.reproducible_peaks/{condition}.consensus.bed"),
                    output:
                        R("6.reproducible_peaks/{condition}.consensus.filtered.bed"),
                    log:
                        R("logs/filter_input_background/{condition}.log"),
                    threads: rthreads("filter_input_background")
                    resources:
                        mem_mb=rmem("filter_input_background"),
                        runtime_min=rruntime("filter_input_background"),
                        runtime_sec=rruntime_sec("filter_input_background"),
                    shell:
                        """
                        ## Keep unflagged consensus sites only; cut returns
                        ## the BED4 shape (the flag column is dropped).
                        awk '$5 == 0' {input.flagged} 2>> {log} \\
                            | cut -f1-4 > {output} 2>> {log}
                        """

        if CONDITIONS_WITHOUT_INPUT:
            rule consensus_peaks:
                # Conditions without input samples keep the exact W7 pipeline
                # and byte-identical BED4 output; restricted to those
                # conditions so it never collides with the flagging rule over
                # the same output pattern.
                input:
                    beds=lambda wc: [R(f"5.callpeak/{s}.pureclip.bed")
                                     for s in condition_ip_samples(wc.condition)],
                output:
                    R("6.reproducible_peaks/{condition}.consensus.bed"),
                wildcard_constraints:
                    condition=_NO_INPUT_CONDITIONS_RE,
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
    else:
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
