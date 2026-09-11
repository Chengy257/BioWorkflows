# SEACR pooled peak calling (peak.caller=seacr; included only when the switch
# selects it): the pooled MACS2 narrow/broad/atac group calls are replaced by
# SEACR (Sparse Enrichment Analysis for CUT&RUN, Yo et al. 2021,
# https://github.com/FredHutch/SEACR), the common CUT&RUN/CUT&Tag alternative
# caller. Per group:
#   1. raw-depth bedGraph coverage of the pooled treat BAMs (samtools merge of
#      the same analysis BAMs the MACS2 rules use, then bamCoverage with
#      --normalizeUsing None and binSize 1 — SEACR integrates the area under
#      the signal curve, so it wants unnormalized, zero-omitted coverage at
#      exact resolution),
#   2. the same for the pooled control BAMs when the group has controls
#      (groups without controls skip the control bedGraph and SEACR falls back
#      to the numeric FDR threshold from the config),
#   3. the SEACR run (external bash script, `bash {CHIP_SEACR} ...`), then the
#      6-column SEACR result (chr, start, end, AUC, max_signal,
#      max_signal_region) is converted to the standard 10-column
#      narrowPeak/broadPeak contract (workflow/scripts/seacr_to_narrowpeak.py)
#      at the SAME {group}_peaks.* path the MACS2 rules would write — FRiP,
#      annotation, blacklist, motif, and DiffBind consumers are unchanged,
#   4. a group coverage bigWig ({group}_FE.bw) from the treat bedGraph via the
#      same bedClip/sort/bedGraphToBigWig recipe as the MACS2 bigwig rule
#      (the bdgcmp fold-enrichment route needs the MACS2 pileup/lambda
#      bedGraphs, which this mode does not produce; the track therefore
#      carries the pooled raw-depth coverage — see the user guide).
# SEACR emits one peak set per group; narrow vs broad naming follows the
# group's assay/peak_type (the same _group_calls_narrow routing the replicate
# stage uses). The replicate/IDR machinery stays MACS2-based:
# caller=seacr + peak.replicate.enabled is rejected at parse time.

_SEACR_GROUPS = list(GROUPS)
_SEACR_CONTROL_GROUPS = [g for g, v in GROUPS.items() if v["control"]]

# SEACR emits one peak set per group; the {peak_suffix} wildcard of the
# seacr_callpeak output matches the group's expected narrowPeak/broadPeak name
# (the same narrow-vs-broad semantics _group_calls_narrow encodes for the
# replicate stage), so group_peak_file()/group_final_peak_file() routing is
# untouched.


rule seacr_bedgraph_treat:
    input:
        bams=lambda wc: group_bams(wc.group, "treat"),
    output:
        R("4.peak/seacr/{group}_treat.bg"),
    wildcard_constraints:
        group=_group_regex(_SEACR_GROUPS),
    params:
        # merged intermediate (removed after coverage); the bamCoverage
        # bedgraph omits zero-coverage bins, as SEACR requires
        tmp_bam=lambda wc: R(f"4.peak/seacr/.{wc.group}_treat.merged.bam"),
    log:
        R("logs/seacr/{group}_treat_bedgraph.log"),
    threads: rthreads("seacr_bedgraph_treat")
    resources:
        mem_mb=rmem("seacr_bedgraph_treat"),
        runtime_min=rruntime("seacr_bedgraph_treat"),
        runtime_sec=rruntime_sec("seacr_bedgraph_treat"),
    shell:
        """
        samtools merge -f -@ {threads} {params.tmp_bam} {input.bams} > {log} 2>&1
        # bamCoverage requires an indexed BAM (WSL real-run finding 2026-09-10)
        samtools index -@ {threads} {params.tmp_bam}
        bamCoverage -b {params.tmp_bam} --outFileFormat bedgraph --normalizeUsing None \
            --binSize 1 -p {threads} -o {output} >> {log} 2>&1
        rm -f {params.tmp_bam} {params.tmp_bam}.bai
        """


rule seacr_bedgraph_control:
    input:
        bams=lambda wc: group_bams(wc.group, "control"),
    output:
        R("4.peak/seacr/{group}_control.bg"),
    wildcard_constraints:
        # only scheduled for groups that actually have controls
        group=_group_regex(_SEACR_CONTROL_GROUPS),
    params:
        tmp_bam=lambda wc: R(f"4.peak/seacr/.{wc.group}_control.merged.bam"),
    log:
        R("logs/seacr/{group}_control_bedgraph.log"),
    threads: rthreads("seacr_bedgraph_control")
    resources:
        mem_mb=rmem("seacr_bedgraph_control"),
        runtime_min=rruntime("seacr_bedgraph_control"),
        runtime_sec=rruntime_sec("seacr_bedgraph_control"),
    shell:
        """
        samtools merge -f -@ {threads} {params.tmp_bam} {input.bams} > {log} 2>&1
        # bamCoverage requires an indexed BAM (WSL real-run finding 2026-09-10)
        samtools index -@ {threads} {params.tmp_bam}
        bamCoverage -b {params.tmp_bam} --outFileFormat bedgraph --normalizeUsing None \
            --binSize 1 -p {threads} -o {output} >> {log} 2>&1
        rm -f {params.tmp_bam} {params.tmp_bam}.bai
        """


rule seacr_callpeak:
    input:
        treat=R("4.peak/seacr/{group}_treat.bg"),
        control=lambda wc: ([R(f"4.peak/seacr/{wc.group}_control.bg")]
                            if GROUPS[wc.group]["control"] else []),
    output:
        peaks=R("4.peak/{group}_peaks.{peak_suffix}"),
    wildcard_constraints:
        group=_group_regex(_SEACR_GROUPS),
        peak_suffix="narrowPeak|broadPeak",
    params:
        seacr=SEACR_BIN,
        # slot 2 of the SEACR CLI takes the control bedgraph, or a numeric FDR
        # threshold when the group has no control; slots 3/4 keep their
        # positions either way (with a numeric threshold the norm/non slot has
        # no control to normalize against and is ignored by SEACR)
        control_arg=lambda wc: (R(f"4.peak/seacr/{wc.group}_control.bg")
                                if GROUPS[wc.group]["control"]
                                else str(PEAK_CALLER["seacr"]["fdr_threshold"])),
        normalize=PEAK_CALLER["seacr"]["normalize"],
        mode=PEAK_CALLER["seacr"]["mode"],
        prefix=lambda wc: R(f"4.peak/seacr/{wc.group}"),
        script=os.path.join(WORKFLOW_DIR, "scripts", "seacr_to_narrowpeak.py"),
        python=PYTHON_BIN,
        name=lambda wc: str(wc.group),
    log:
        # carries {peak_suffix} because every log/output file of a rule must
        # expose the same wildcard set (snakemake requirement)
        R("logs/seacr/{group}_callpeak.{peak_suffix}.log"),
    threads: rthreads("seacr_callpeak")  # SEACR is bash/awk single-threaded
    resources:
        mem_mb=rmem("seacr_callpeak"),
        runtime_min=rruntime("seacr_callpeak"),
        runtime_sec=rruntime_sec("seacr_callpeak"),
    shell:
        """
        bash {params.seacr} {input.treat} {params.control_arg} {params.normalize} \
            {params.mode} {params.prefix} > {log} 2>&1
        # SEACR_1.3.sh renames its merged intermediate to the mode-named bed
        # (stringent/relaxed) and deletes the merge file — the usage text still
        # names the merge file, the code is authoritative (final-review finding
        # 2026-09-10)
        {params.python} {params.script} {params.prefix}.{params.mode}.bed \
            {output.peaks} {params.name} >> {log} 2>&1
        rm -f {params.prefix}.auc.threshold.bed {params.prefix}.{params.mode}.bed
        """


rule seacr_bigwig:
    # Group coverage track under the SEACR caller: the treat bedGraph clipped
    # to chromosome bounds and converted with the exact same
    # clip/sort/bedGraphToBigWig recipe the MACS2 bigwig rule uses (the
    # filename stays {group}_FE.bw to keep the deliverable layout stable; the
    # values are pooled raw depth, not MACS2 fold enrichment).
    input:
        bg=R("4.peak/seacr/{group}_treat.bg"),
        chromsize=config["chromsize"],
    output:
        R("4.peak/{group}_FE.bw"),
    wildcard_constraints:
        group=_group_regex(_SEACR_GROUPS),
    params:
        outdir=lambda wc, output: os.path.dirname(str(output)),
    log:
        R("logs/seacr/{group}_bigwig.log"),
    threads: rthreads("seacr_bigwig")
    resources:
        mem_mb=rmem("seacr_bigwig"),
        runtime_min=rruntime("seacr_bigwig"),
        runtime_sec=rruntime_sec("seacr_bigwig"),
    shell:
        """
        bedtools slop -i {input.bg} -g {input.chromsize} -b 0 \
            | bedClip stdin {input.chromsize} {params.outdir}/{wildcards.group}_FE.clip >> {log} 2>&1
        LC_COLLATE=C sort -k1,1 -k2,2n {params.outdir}/{wildcards.group}_FE.clip \
            > {params.outdir}/{wildcards.group}_FE.clip.sorted
        bedGraphToBigWig {params.outdir}/{wildcards.group}_FE.clip.sorted {input.chromsize} {output} >> {log} 2>&1
        rm -f {params.outdir}/{wildcards.group}_FE.clip {params.outdir}/{wildcards.group}_FE.clip.sorted
        """
