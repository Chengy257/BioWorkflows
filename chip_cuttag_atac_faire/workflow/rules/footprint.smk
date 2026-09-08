# TOBIAS footprinting on the ATAC/FAIRE groups (v0.6; included only when
# footprint.enabled is true). Footprinting models the Tn5 cut-site bias of
# open-chromatin assays, so the rules enumerate the atac/faire groups only —
# chip/cuttag enrichment groups carry no usable footprint signal. Per group:
#   1. ATACorrect: the group's treat analysis BAMs are pooled (samtools merge,
#      the same recipe the SEACR bedGraph rules use; the merged BAM is a
#      dot-prefixed temp removed after the run) and bias-corrected against the
#      reference genome over the group's final peak set
#      (group_final_peak_file(): IDR/consensus when the replicate stage is on,
#      pooled otherwise) -> {group}_corrected.bw plus QC plots,
#   2. ScoreBigwig: footprint scores over the same peak regions from the
#      corrected cutsite track (TOBIAS 0.8 API: --signal/--regions/--output;
#      motif scanning is NOT part of this step) -> one score bigWig,
#   3. BINDetect (optional, footprint.bindetect): per-TF binding detection
#      over the same corrected track and peak set (--skip-excel keeps the
#      output dependency-free).
# TOBIAS is an external install with a heavyweight dependency set, resolved
# like idr/HOMER/SEACR via the software.yaml paths: mechanism (CHIP_TOBIAS).
# Deliverables land under results/7.footprint/{group}/: the ATACorrect rule
# declares the bias-corrected bigWig as a tracked FILE (snakemake cannot link
# inputs that live inside another rule's directory() output, and the corrected
# track feeds both downstream stages); the QC plots ATACorrect drops next to it
# stay untracked; ScoreBigwig and BINDetect declare one directory output each
# (motif.smk precedent) because nothing consumes files inside them.

rule footprint_ataccorrect:
    input:
        bams=lambda wc: group_bams(wc.group, "treat"),
        peaks=lambda wc: group_final_peak_file(wc.group),
        genome=config["genome_fa"],
    output:
        bw=R("7.footprint/{group}/ataccorrect/{group}_corrected.bw"),
    wildcard_constraints:
        group=_group_regex(FOOTPRINT_GROUPS),
    params:
        tobias=TOBIAS_BIN,
        outdir=lambda wc: R(f"7.footprint/{wc.group}/ataccorrect"),
        # pooled temp BAM (removed after the job; the dot prefix keeps it out
        # of the deliverable tree)
        tmp_bam=lambda wc: R(f"7.footprint/.{wc.group}_pooled_treat.bam"),
    log:
        R("logs/footprint/{group}_ataccorrect.log"),
    threads: rthreads("footprint_ataccorrect")
    resources:
        mem_mb=rmem("footprint_ataccorrect"),
        runtime_min=rruntime("footprint_ataccorrect"),
        runtime_sec=rruntime_sec("footprint_ataccorrect"),
    shell:
        """
        samtools merge -f -@ {threads} {params.tmp_bam} {input.bams} > {log} 2>&1
        samtools index -@ {threads} {params.tmp_bam}
        {params.tobias} ATACorrect --bam {params.tmp_bam} --peaks {input.peaks} \
            --genome {input.genome} --outdir {params.outdir} --prefix {wildcards.group} \
            --cores {threads} >> {log} 2>&1
        rm -f {params.tmp_bam} {params.tmp_bam}.bai
        """


rule footprint_scorebigwig:
    # TOBIAS 0.8 API: ScoreBigwig scores cutsite tracks over regions
    # (--signal/--regions/--output); it takes no motif input — motif scanning
    # lives in BINDetect (footprint.motifs feeds that rule)
    input:
        bw=R("7.footprint/{group}/ataccorrect/{group}_corrected.bw"),
        peaks=lambda wc: group_final_peak_file(wc.group),
    output:
        R("7.footprint/{group}/scorebigwig/{group}_footprint_scores.bw"),
    wildcard_constraints:
        group=_group_regex(FOOTPRINT_GROUPS),
    params:
        tobias=TOBIAS_BIN,
        outdir=lambda wc, output: os.path.dirname(str(output)),
    log:
        R("logs/footprint/{group}_scorebigwig.log"),
    threads: rthreads("footprint_scorebigwig")  # ScoreBigwig runs single-threaded
    resources:
        mem_mb=rmem("footprint_scorebigwig"),
        runtime_min=rruntime("footprint_scorebigwig"),
        runtime_sec=rruntime_sec("footprint_scorebigwig"),
    shell:
        """
        mkdir -p {params.outdir}
        {params.tobias} ScoreBigwig --signal {input.bw} --regions {input.peaks} \
            --output {output} > {log} 2>&1
        """


rule footprint_bindetect:
    # Gated through the targets (FOOTPRINT_TARGETS only aggregates the
    # bindetect directories when footprint.bindetect is on), so with the
    # switch off the rule simply has no instance in the DAG.
    input:
        bw=R("7.footprint/{group}/ataccorrect/{group}_corrected.bw"),
        peaks=lambda wc: group_final_peak_file(wc.group),
        genome=config["genome_fa"],
        motifs=FOOTPRINT["motifs"],
    output:
        directory(R("7.footprint/{group}/bindetect")),
    wildcard_constraints:
        group=_group_regex(FOOTPRINT_GROUPS),
    params:
        tobias=TOBIAS_BIN,
        pvalue_arg=FOOTPRINT_MOTIF_PVALUE_ARG,
    log:
        R("logs/footprint/{group}_bindetect.log"),
    threads: rthreads("footprint_bindetect")
    resources:
        mem_mb=rmem("footprint_bindetect"),
        runtime_min=rruntime("footprint_bindetect"),
        runtime_sec=rruntime_sec("footprint_bindetect"),
    shell:
        """
        {params.tobias} BINDetect --signals {input.bw} --peaks {input.peaks} \
            --motifs {input.motifs} --genome {input.genome} \
            --cond-names {wildcards.group} --outdir {output} --cores {threads} \
            --skip-excel {params.pvalue_arg} > {log} 2>&1
        """
