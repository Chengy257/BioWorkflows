# Per-class counting, matrix merging, and cascade read-fate summary.
# Depends on workflow/rules/common.smk for SCRIPTS/_CLASSES/helpers.

rule count_stage:
    input:
        sam=R("3.align/filter/{klass}/{sample}.sam"),
    output:
        R("4.expression/{klass}/{sample}_counts.txt"),
    params:
        script=os.path.join(SCRIPTS, "count_features.py"),
    log:
        R("logs/count/{klass}/{sample}.log"),
    threads: rthreads("count_stage")
    resources:
        mem_mb=rmem("count_stage"),
        runtime_min=rruntime("count_stage"),
        runtime_sec=rruntime_sec("count_stage"),
    shell:
        "python3 {params.script} {input.sam} -o {output} > {log} 2>&1"


rule raw_count:
    input:
        fq=lambda wc: raw_fastq(wc.sample),
    output:
        R("5.QC/raw_counts/{sample}.txt"),
    log:
        R("logs/raw_count/{sample}.log"),
    threads: 1
    resources:
        mem_mb=1024,
        runtime_min=10,
        runtime_sec=600,
    shell:
        "zcat {input.fq} | awk 'END{{print NR/4}}' > {output} 2> {log}"


rule merge_counts:
    input:
        counts=expand(R("4.expression/{klass}/{sample}_counts.txt"),
                      klass=_CLASSES, sample=SAMPLES),
    output:
        all_classes=R("4.expression/all_classes_counts.tsv"),
        matrices=expand(R("4.expression/{klass}/{klass}_counts.tsv"), klass=_CLASSES),
        rpm_matrices=expand(R("4.expression/{klass}/{klass}_RPM.tsv"), klass=_CLASSES),
    params:
        script=os.path.join(SCRIPTS, "merge_counts.py"),
        indir=f"{RD}4.expression",
        classes=" ".join(_CLASSES),
        samples=" ".join(SAMPLES),
    log:
        R("logs/merge_counts.log"),
    threads: rthreads("merge_counts")
    resources:
        mem_mb=rmem("merge_counts"),
        runtime_min=rruntime("merge_counts"),
        runtime_sec=rruntime_sec("merge_counts"),
    shell:
        "python3 {params.script} --indir {params.indir} --classes {params.classes} "
        "--samples {params.samples} > {log} 2>&1"


rule cascade_summary:
    input:
        trim_reports=expand(R("2.cleandata/{sample}_trimming_report.txt"), sample=SAMPLES),
        raw_counts=expand(R("5.QC/raw_counts/{sample}.txt"), sample=SAMPLES),
        counts=expand(R("4.expression/{klass}/{sample}_counts.txt"),
                      klass=_CLASSES, sample=SAMPLES) if _CLASSES else [],
        # Read by path in cascade_summary.py (genome_unmapped column); declares
        # the DAG edge so the summary is ordered after genome_align.
        genome_unmapped=expand(R("3.align/genome/{sample}_unmapped.fq"),
                               sample=SAMPLES) if _GENOME_CONFIGURED else [],
    output:
        summary=R("5.QC/cascade_summary.tsv"),
        mqc=R("5.QC/cascade_summary_mqc.tsv"),
    params:
        script=os.path.join(SCRIPTS, "cascade_summary.py"),
        expr_dir=f"{RD}4.expression",
        qc_dir=f"{RD}5.QC",
        classes=" ".join(_CLASSES),
        samples=" ".join(SAMPLES),
        genome_configured="true" if _GENOME_CONFIGURED else "false",
        raw_count_args=lambda wc: " ".join(
            open(f"{RD}5.QC/raw_counts/{s}.txt").read().strip() for s in SAMPLES),
    log:
        R("logs/cascade_summary.log"),
    threads: rthreads("cascade_summary")
    resources:
        mem_mb=rmem("cascade_summary"),
        runtime_min=rruntime("cascade_summary"),
        runtime_sec=rruntime_sec("cascade_summary"),
    shell:
        "python3 {params.script} --indir {params.expr_dir} --classes {params.classes} "
        "--samples {params.samples} --trim-reports {input.trim_reports} "
        "--raw-counts {params.raw_count_args} --outdir {params.qc_dir} "
        "--genome-configured {params.genome_configured} > {log} 2>&1"
