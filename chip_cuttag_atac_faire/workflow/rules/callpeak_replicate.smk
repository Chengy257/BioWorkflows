# Replicate-aware peak stage (v0.5; included only when peak.replicate.enabled):
#   1. per-replicate peak calling for every treat sample (relaxed narrow
#      cutoff, the group's pooled control as the MACS2 -c input — ENCODE
#      style; no pileup/summit side outputs, they exist only for the pooled
#      calls),
#   2. pairwise IDR for narrow groups with >= 2 treats (classic idr tool,
#      rank/threshold from config),
#   3. final reproducible peak set per group: exactly 2 treats -> the single
#      pair result; > 2 treats -> union of pairwise IDR peaks kept where the
#      pairwise support >= consensus_min_replicates (a deliberate, documented
#      simplification of the ENCODE rescue/self-consistency scheme); broad
#      groups -> bedtools multiinter consensus over the per-replicate broad
#      peaks with the same support cutoff (multiinter column 4 = number of
#      input files carrying the feature),
#   4. a project-wide summary table (per-group replicate peak counts, final
#      peak count, retained fraction) injected into MultiQC.
# The pooled peak files (4.peak/{group}_peaks.*) keep their exact names and
# stay the FRiP basis unless peak.replicate.frip_on switches to consensus.

_NARROW_REP_TREATS = [s for g in NARROW_REP_GROUPS for s in REPLICATE_TREATS[g]]
_BROAD_REP_TREATS = [s for g in BROAD_REP_GROUPS for s in REPLICATE_TREATS[g]]
_ALL_PAIR_SLUGS = [idr_pair_slug(a, b)
                   for g in IDR_GROUPS for a, b in IDR_PAIRS[g]]


rule callpeak_narrow_replicate:
    input:
        treat=lambda wc: sample_bam(wc.sample),
        control=lambda wc: group_bams(wc.group, "control"),
    output:
        R("4.peak/replicates/{group}/{sample}_peaks.narrowPeak"),
    wildcard_constraints:
        group=_group_regex(NARROW_REP_GROUPS),
        sample=_group_regex(_NARROW_REP_TREATS),
    params:
        control=group_control_arg,
        # same format/nomodel routing as the pooled narrow/atac rules
        fmt=lambda wc: ("BAMPE" if GROUPS[wc.group]["layout"] == "PE" else "BAM")
        if (GROUPS[wc.group]["seqtype"] not in ("atac", "faire")
            or config["peak"]["atac"]["mode"] == "bampe") else "BAM",
        nomodel=lambda wc: "" if (GROUPS[wc.group]["seqtype"] not in ("atac", "faire")
                                  or config["peak"]["atac"]["mode"] == "bampe") else
            "--nomodel --shift {} --extsize {}".format(
                config["peak"]["atac"]["shift"], config["peak"]["atac"]["extsize"]),
        gsize=config["genome_size"],
        keepdup=config["peak"]["keepdup"],
        qvalue=REPLICATE["qvalue"],
        outdir=lambda wc, output: os.path.dirname(str(output)),
    log:
        R("logs/callpeak_replicate/{group}__{sample}_narrow.log"),
    threads: rthreads("callpeak_narrow_replicate")  # MACS2 is single-threaded
    resources:
        mem_mb=rmem("callpeak_narrow_replicate"),
        runtime_min=rruntime("callpeak_narrow_replicate"),
        runtime_sec=rruntime_sec("callpeak_narrow_replicate"),
    shell:
        """
        macs2 callpeak \
            -t {input.treat} {params.control} \
            -f {params.fmt} -g {params.gsize} \
            --keep-dup {params.keepdup} -q {params.qvalue} \
            {params.nomodel} \
            --outdir {params.outdir} -n {wildcards.sample} > {log} 2>&1
        # the unmanaged _peaks.xls side artifact MACS2 always writes is
        # removed to keep 4.peak tidy (only the narrowPeak output is part
        # of the contract)
        rm -f {params.outdir}/{wildcards.sample}_peaks.xls
        """


rule callpeak_broad_replicate:
    input:
        treat=lambda wc: sample_bam(wc.sample),
        control=lambda wc: group_bams(wc.group, "control"),
    output:
        R("4.peak/replicates/{group}/{sample}_peaks.broadPeak"),
    wildcard_constraints:
        group=_group_regex(BROAD_REP_GROUPS),
        sample=_group_regex(_BROAD_REP_TREATS),
    params:
        control=group_control_arg,
        fmt=lambda wc: "BAMPE" if GROUPS[wc.group]["layout"] == "PE" else "BAM",
        gsize=config["genome_size"],
        keepdup=config["peak"]["keepdup"],
        broad_cutoff=config["peak"]["broad_cutoff"],
        outdir=lambda wc, output: os.path.dirname(str(output)),
    log:
        R("logs/callpeak_replicate/{group}__{sample}_broad.log"),
    threads: rthreads("callpeak_broad_replicate")  # MACS2 is single-threaded
    resources:
        mem_mb=rmem("callpeak_broad_replicate"),
        runtime_min=rruntime("callpeak_broad_replicate"),
        runtime_sec=rruntime_sec("callpeak_broad_replicate"),
    shell:
        """
        macs2 callpeak \
            -t {input.treat} {params.control} \
            -f {params.fmt} -g {params.gsize} \
            --keep-dup {params.keepdup} --broad --broad-cutoff {params.broad_cutoff} \
            --outdir {params.outdir} -n {wildcards.sample} > {log} 2>&1
        rm -f {params.outdir}/{wildcards.sample}_peaks.xls
        """


rule idr_pair:
    input:
        a=lambda wc: replicate_peak_file(wc.group, parse_idr_pair_slug(wc.pair)[0]),
        b=lambda wc: replicate_peak_file(wc.group, parse_idr_pair_slug(wc.pair)[1]),
    output:
        R("4.peak/idr/{group}/{pair}.narrowPeak"),
    wildcard_constraints:
        group=_group_regex(IDR_GROUPS),
        pair=_group_regex(_ALL_PAIR_SLUGS),
    params:
        idr=IDR_BIN,
        rank=REPLICATE["idr_rank"],
        threshold=REPLICATE["idr_threshold"],
        full=lambda wc, output: str(output) + ".idr.full",
    log:
        R("logs/idr/{group}__{pair}.log"),
    threads: rthreads("idr_pair")
    resources:
        mem_mb=rmem("idr_pair"),
        runtime_min=rruntime("idr_pair"),
        runtime_sec=rruntime_sec("idr_pair"),
    shell:
        """
        {params.idr} --samples {input.a} {input.b} \
            --input-file-type narrowPeak --rank {params.rank} \
            --idr-threshold {params.threshold} \
            --output-file {params.full} --log-output-file {log} > {log} 2>&1
        # idr appends its per-peak IDR columns after the narrowPeak columns;
        # trim back to the 10-column contract every downstream consumer uses
        awk -v OFS="\\t" '{{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' \
            {params.full} > {output} 2>> {log}
        rm -f {params.full}
        """


rule idr_final:
    input:
        pairs=lambda wc: [idr_pair_file(wc.group, a, b)
                          for a, b in IDR_PAIRS[wc.group]],
    output:
        peaks=R("4.peak/{group}_IDR_peaks.narrowPeak"),
        support=R("4.peak/{group}_IDR_support.bed"),
    wildcard_constraints:
        group=_group_regex(IDR_GROUPS),
    params:
        npairs=lambda wc: len(IDR_PAIRS[wc.group]),
        minrep=REPLICATE["consensus_min_replicates"],
    log:
        R("logs/idr/{group}__final.log"),
    threads: rthreads("idr_final")
    resources:
        mem_mb=rmem("idr_final"),
        runtime_min=rruntime("idr_final"),
        runtime_sec=rruntime_sec("idr_final"),
    shell:
        """
        tmpdir={resources.tmpdir}/idr_final_{wildcards.group}
        rm -rf "$tmpdir" && mkdir -p "$tmpdir"
        if [ "{params.npairs}" -eq 1 ]; then
            awk -v OFS="\\t" '{{print $1,$2,$3,2}}' {input.pairs[0]} > {output.support}
            cp {input.pairs[0]} {output.peaks}
        else
            # multiinter needs coordinate-sorted inputs; sort each pairwise
            # result into the scheduler-provided tmpdir
            for peaks in {input.pairs}; do
                LC_COLLATE=C sort -k1,1 -k2,2n "$peaks" > "$tmpdir/$(basename "$peaks").sorted"
            done
            # multiinter columns: chrom, start, end, support (number of input
            # files carrying the feature), membership lists; column 4 is the
            # replicate support used for the cutoff
            bedtools multiinter -i "$tmpdir"/*.sorted > "$tmpdir/multiinter.tsv" 2>> {log}
            awk -v OFS="\\t" -v m={params.minrep} '$4 >= m {{print $1,$2,$3,$4}}' \
                "$tmpdir/multiinter.tsv" > {output.support} 2>> {log}
            # fill the kept intervals into a standard 10-column narrowPeak
            # (name/score/signal carry the support count; statistical columns
            # are placeholders — consensus intervals carry no MACS2 statistics)
            awk -v OFS="\\t" '{{n++; print $1,$2,$3,"idr_" n,$4,".",$4,0,0,0}}' \
                {output.support} > {output.peaks} 2>> {log}
        fi
        rm -rf "$tmpdir"
        """


rule broad_consensus:
    input:
        lambda wc: [replicate_peak_file(wc.group, s)
                    for s in REPLICATE_TREATS[wc.group]],
    output:
        peaks=R("4.peak/{group}_consensus_peaks.broadPeak"),
        support=R("4.peak/{group}_consensus_support.bed"),
    wildcard_constraints:
        group=_group_regex(BROAD_CONSENSUS_GROUPS),
    params:
        minrep=REPLICATE["consensus_min_replicates"],
    log:
        R("logs/broad_consensus/{group}.log"),
    threads: rthreads("broad_consensus")
    resources:
        mem_mb=rmem("broad_consensus"),
        runtime_min=rruntime("broad_consensus"),
        runtime_sec=rruntime_sec("broad_consensus"),
    shell:
        """
        tmpdir={resources.tmpdir}/broad_consensus_{wildcards.group}
        rm -rf "$tmpdir" && mkdir -p "$tmpdir"
        for peaks in {input}; do
            LC_COLLATE=C sort -k1,1 -k2,2n "$peaks" > "$tmpdir/$(basename "$peaks").sorted"
        done
        bedtools multiinter -i "$tmpdir"/*.sorted > "$tmpdir/multiinter.tsv" 2> {log}
        awk -v OFS="\\t" -v m={params.minrep} '$4 >= m {{print $1,$2,$3,$4}}' \
            "$tmpdir/multiinter.tsv" > {output.support} 2>> {log}
        # standard 9-column broadPeak with the support count as name/score/signal
        awk -v OFS="\\t" '{{n++; print $1,$2,$3,"cons_" n,$4,".",$4,0,0}}' \
            {output.support} > {output.peaks} 2>> {log}
        rm -rf "$tmpdir"
        """


rule replicate_summary:
    input:
        samples=lambda wc: _resolve_sample_table(config["grouplist"]),
    output:
        tsv=R("5.QC/replicate_peaks/Replicate_summary.tsv"),
        mqc=R("5.QC/replicate_peaks/Replicate_summary_mqc.tsv"),
    params:
        script=os.path.join(WORKFLOW_DIR, "scripts", "replicate_summary.py"),
        python=PYTHON_BIN,
        results_dir=RD.rstrip(os.sep),
    log:
        R("logs/replicate_summary.log"),
    threads: rthreads("replicate_summary")
    resources:
        mem_mb=rmem("replicate_summary"),
        runtime_min=rruntime("replicate_summary"),
        runtime_sec=rruntime_sec("replicate_summary"),
    shell:
        """
        {params.python} {params.script} \
            --samples {input.samples} --results-dir {params.results_dir} \
            --out {output.tsv} --mqc {output.mqc} > {log} 2>&1
        """
