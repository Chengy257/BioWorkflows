# QC gate summary (qc.gates; included only when the switch is on): one
# PASS/WARN/FAIL row per sample aggregating the QC metrics the workflow
# already computes (mapping rate, duplication rate, FRiP, NSC/RSC, TSS
# enrichment, organelle fraction). Strictly informational — the pipeline
# never hard-fails on a gate; metric sources and semantics live in
# workflow/scripts/gates_summary.py.

# FRiP outputs exist for every treat AND control sample (frip.smk); the
# summary resolves them through the sample=group tokens baked into params.
GATES_FRIP_TSVS = [R(f"5.QC/frip/{g}__{s}.frip.tsv")
                   for g in GROUPS
                   for s in GROUPS[g]["treat"] + GROUPS[g]["control"]]

GATES_THRESHOLD_KEYS = ("mapping_rate_min", "dup_rate_max", "frip_min",
                        "nsc_min", "rsc_min", "tss_min", "organelle_max")

rule gates_flagstat:
    input:
        lambda wc: sample_bam(wc.sample),
    output:
        R("5.QC/gates/{sample}_flagstat.txt"),
    wildcard_constraints:
        sample=_group_regex(SAMPLES),
    log:
        R("logs/gates/{sample}_flagstat.log"),
    threads: rthreads("gates_flagstat")
    resources:
        mem_mb=rmem("gates_flagstat"),
        runtime_min=rruntime("gates_flagstat"),
        runtime_sec=rruntime_sec("gates_flagstat"),
    shell:
        """
        samtools flagstat -@ {threads} {input} > {output} 2> {log}
        """

rule qc_gates:
    input:
        flagstat=expand(R("5.QC/gates/{sample}_flagstat.txt"), sample=SAMPLES),
        # dup metrics exist only for assays that deduplicate (cuttag keeps
        # duplicates, so dup_rate renders NA there)
        dup_metrics=expand(R("3.align/bowtie2/{sample}_dup_metrics.txt"),
                           sample=[s for s in SAMPLES if assay_needs_dedup(SEQTYPE_OF[s])]),
        frip=lambda wc: (GATES_FRIP_TSVS if config["qc"]["frip"] else []),
        spp=lambda wc: (expand(R("5.QC/spp/{sample}_NSC.txt"), sample=SAMPLES)
                        + expand(R("5.QC/spp/{sample}_RSC.txt"), sample=SAMPLES)
                        if config["qc"]["nsc_rsc"] else []),
        tss=lambda wc: (expand(R("5.QC/tss/{sample}_TSSE.txt"), sample=TSS_SAMPLES)
                        if (QC_TSS and TSS_SAMPLES) else []),
        organelle=lambda wc: ([R("5.QC/organelle/Organelle_summary.tsv")]
                              if QC_ORGANELLE else []),
    output:
        tsv=R("5.QC/gates/gate_summary.tsv"),
        mqc=R("5.QC/gates/gate_summary_mqc.tsv"),
    params:
        script=os.path.join(WORKFLOW_DIR, "scripts", "gates_summary.py"),
        python=PYTHON_BIN,
        # sample=group tokens (the group resolves each FRiP filename)
        samples=" ".join(f"{s}={GATES_SAMPLE_GROUPS[s]}" for s in SAMPLES),
        align_dir=R("3.align/bowtie2"),
        qc_dir=R("5.QC"),
        gates_dir=R("5.QC/gates"),
        thresholds=" ".join(f"{k}={float(GATES['thresholds'][k])}"
                            for k in GATES_THRESHOLD_KEYS),
    log:
        R("logs/gates/summary.log"),
    threads: rthreads("qc_gates")
    resources:
        mem_mb=rmem("qc_gates"),
        runtime_min=rruntime("qc_gates"),
        runtime_sec=rruntime_sec("qc_gates"),
    shell:
        """
        {params.python} {params.script} \
            --samples {params.samples} \
            --align-dir {params.align_dir} --qc-dir {params.qc_dir} \
            --gates-dir {params.gates_dir} \
            --out {output.tsv} --mqc {output.mqc} \
            --thresholds {params.thresholds} > {log} 2>&1
        """
