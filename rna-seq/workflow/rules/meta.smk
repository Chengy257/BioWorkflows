# Workflow metadata and software-version capture.

rule software_versions:
    output:
        R("software_versions.yaml"),
    log:
        R("logs/software_versions.log.txt"),
    params:
        script=os.path.join(SCRIPTS, "collect_versions.py"),
        python=PYTHON,
        software_config=runtime_env("RNASEQ_SOFTWARE_CONFIG", os.path.join(BASE_DIR, "config", "software.yaml")),
        workflow_dir=WORKFLOW_DIR,
        snakemake_version=__import__("importlib.metadata", fromlist=["version"]).version("snakemake"),
    threads:
        rthreads("software_versions")
    resources:
        mem_mb=rmem("software_versions"),
        runtime_min=rruntime("software_versions"),
        runtime_sec=rruntime_sec("software_versions"),
    shell:
        "{params.python} {params.script} --software-config {params.software_config} --workflow {params.workflow_dir} "
        "--snakemake-version {params.snakemake_version} --out {output} >> {log} 2>&1"
