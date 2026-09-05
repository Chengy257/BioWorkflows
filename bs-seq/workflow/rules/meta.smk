# Run metadata: version capture + whole-pipeline MultiQC summary.
#
# software_versions runs workflow/scripts/collect_versions.py through
# params+shell (never the script: directive) with the BSSEQ_ environment
# overrides; multiqc reports on the trim reports + fastqc zips (only when
# trimming is enabled) and the Bismark alignment/dedup/splitting reports,
# using the workflow multiqc_config.yaml.


def _multiqc_inputs():
    """Flat input list for the whole-run MultiQC report: the per-sample
    Trim Galore trimming reports and fastqc zips (layout-aware, and only
    when trim.enabled), plus the Bismark alignment, dedup, and
    methylation-extractor (splitting) reports of every declared sample."""
    files = []
    for sample in SAMPLES:
        if TRIM_ENABLED:
            if layout_of(sample) == "PE":
                files += [
                    R(f"2.cleandata/{sample}_1_val_1.fq.gz_trimming_report.txt"),
                    R(f"2.cleandata/{sample}_2_val_2.fq.gz_trimming_report.txt"),
                    R(f"2.cleandata/fastqc/{sample}_1_val_1_fastqc.zip"),
                    R(f"2.cleandata/fastqc/{sample}_2_val_2_fastqc.zip"),
                ]
            else:
                files += [
                    R(f"2.cleandata/{sample}_trimming_report.txt"),
                    R(f"2.cleandata/fastqc/{sample}_trimmed_fastqc.zip"),
                ]
        files += [
            R(f"3.align/{sample}_report.txt"),
            R(f"4.dedup/{sample}.deduplication_report.txt"),
            R(f"5.methylation/{sample}/{sample}.deduplicated_splitting_report.txt"),
        ]
    return files


rule software_versions:
    output:
        R("5.QC/software_versions.yaml"),
    log:
        R("5.QC/logs/software_versions.log.txt"),
    params:
        script=os.path.join(WORKFLOW_DIR, "scripts", "collect_versions.py"),
        software_config=os.environ.get("BSSEQ_SOFTWARE_CONFIG",
                                       os.path.join(BASE_DIR, "config", "software.yaml")),
        workflow_dir=WORKFLOW_DIR,
        python=os.environ.get("BSSEQ_PYTHON", "python3"),
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
        reports=lambda wc: _multiqc_inputs(),
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
