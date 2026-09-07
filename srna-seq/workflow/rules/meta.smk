# Run metadata: version capture + whole-pipeline MultiQC summary.

rule software_versions:
    output:
        R("5.QC/software_versions.yaml"),
    log:
        R("5.QC/logs/software_versions.log.txt"),
    params:
        script=os.path.join(WORKFLOW_DIR, "scripts", "collect_versions.py"),
        software_config=os.environ.get("SRNA_SOFTWARE_CONFIG",
                                       os.path.join(BASE_DIR, "config", "software.yaml")),
        workflow_dir=WORKFLOW_DIR,
        python=os.environ.get("SRNA_PYTHON", "python3"),
        snakemake_version=__import__("importlib.metadata", fromlist=["version"]).version("snakemake"),
    threads: rthreads("software_versions")
    resources:
        mem_mb=rmem("software_versions"),
        runtime_min=rruntime("software_versions"),
        runtime_sec=rruntime_sec("software_versions"),
    shell:
        "{params.python} {params.script} --software-config {params.software_config} "
        "--workflow {params.workflow_dir} --snakemake-version {params.snakemake_version} "
        "--out {output} >> {log} 2>&1"


rule multiqc:
    input:
        fastqc=expand(R("2.cleandata/fastqc/{sample}_trimmed_fastqc.zip"), sample=SAMPLES),
        trim_reports=expand(R("2.cleandata/{sample}_trimming_report.txt"), sample=SAMPLES),
        cascade_mqc=R("5.QC/cascade_summary_mqc.tsv"),
    output:
        R("5.QC/multiqc/multiqc_report.html"),
    params:
        mqc_config=os.path.join(WORKFLOW_DIR, "multiqc_config.yaml"),
        outdir=lambda wc, output: os.path.dirname(str(output)),
    log:
        R("logs/multiqc.log"),
    threads: rthreads("multiqc")
    resources:
        mem_mb=rmem("multiqc"),
        runtime_min=rruntime("multiqc"),
        runtime_sec=rruntime_sec("multiqc"),
    shell:
        """
        multiqc --force -o {params.outdir} -c {params.mqc_config} \
            --filename multiqc_report.html {input} > {log} 2>&1
        """
