# DiffBind differential binding (v0.5 Phase 3; included only when
# diffbind.enabled is true). Per contrast [groupA, groupB]:
#   1. diffbind_sheet builds the DiffBind sample sheet from the workflow
#      sample table (Condition from the optional condition column, falling
#      back to the group name; optional Batch column carried through; the
#      per-sample peak file is the replicate call when the replicate stage is
#      enabled — recommended — else the pooled group set),
#   2. run_diffbind.R counts reads over the consensus peak set (summit-
#      recentered when diffbind.summit_flank > 0), fits DESeq2 or edgeR
#      (optionally blocked on Batch), and writes the DB tables + MA/volcano/
#      PCA plots into results/6.diffbind/{A}__vs__{B}/.
# Each contrast arm must carry >= 2 treat replicates (parse-time validation
# in common.smk); controls attach as the DiffBind background when
# diffbind.use_controls is true and the group has exactly one control.

rule diffbind_sheet:
    input:
        samples=lambda wc: _resolve_sample_table(config["grouplist"]),
        bams=lambda wc: [sample_bam(s)
                         for g in _contrast_groups(wc.contrast)
                         for s in GROUPS[g]["treat"]],
        peaks=lambda wc: [diffbind_sample_peaks(g, s)
                          for g in _contrast_groups(wc.contrast)
                          for s in GROUPS[g]["treat"]],
    output:
        R("6.diffbind/{contrast}/samplesheet.tsv"),
    wildcard_constraints:
        contrast=_group_regex([slug for _a, _b, slug in DIFFBIND_CONTRASTS]),
    params:
        script=os.path.join(WORKFLOW_DIR, "scripts", "diffbind_sheet.py"),
        python=PYTHON_BIN,
        results_dir=RD.rstrip(os.sep),
        groups=lambda wc: " ".join(_contrast_groups(wc.contrast)),
        replicates="--replicates" if REPLICATE["enabled"] else "",
        use_controls="--use-controls" if DIFFBIND["use_controls"] else "",
    log:
        R("logs/diffbind/{contrast}_sheet.log"),
    threads: rthreads("diffbind_sheet")
    resources:
        mem_mb=rmem("diffbind_sheet"),
        runtime_min=rruntime("diffbind_sheet"),
        runtime_sec=rruntime_sec("diffbind_sheet"),
    shell:
        """
        {params.python} {params.script} \
            --samples {input.samples} --results-dir {params.results_dir} \
            --groups {params.groups} {params.replicates} {params.use_controls} \
            --out {output} > {log} 2>&1
        """


rule diffbind_report:
    input:
        sheet=R("6.diffbind/{contrast}/samplesheet.tsv"),
        bams=lambda wc: [sample_bam(s)
                         for g in _contrast_groups(wc.contrast)
                         for s in GROUPS[g]["treat"]],
        peaks=lambda wc: [diffbind_sample_peaks(g, s)
                          for g in _contrast_groups(wc.contrast)
                          for s in GROUPS[g]["treat"]],
    output:
        tsv=R("6.diffbind/{contrast}/DB_results.tsv"),
        sig=R("6.diffbind/{contrast}/DB_significant.tsv"),
    wildcard_constraints:
        contrast=_group_regex([slug for _a, _b, slug in DIFFBIND_CONTRASTS]),
    params:
        script=os.path.join(WORKFLOW_DIR, "scripts", "run_diffbind.R"),
        rscript=RSCRIPT,
        analysis=DIFFBIND["analysis"],
        summit=DIFFBIND["summit_flank"],
        fdr=DIFFBIND["fdr"],
        fc=DIFFBIND["foldchange"],
        block="T" if DIFFBIND["batch_correction"] else "F",
    log:
        R("logs/diffbind/{contrast}.log"),
    threads: rthreads("diffbind_report")
    resources:
        mem_mb=rmem("diffbind_report"),
        runtime_min=rruntime("diffbind_report"),
        runtime_sec=rruntime_sec("diffbind_report"),
    shell:
        """
        {params.rscript} {params.script} \
            {input.sheet} {params.analysis} {params.summit} \
            {params.fdr} {params.fc} {params.block} > {log} 2>&1
        """
