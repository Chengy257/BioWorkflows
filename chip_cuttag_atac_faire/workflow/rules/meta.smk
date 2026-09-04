# Run metadata: record the actually-used software versions.

rule software_versions:
    output:
        R("5.QC/software_versions.yaml"),
    log:
        R("5.QC/logs/software_versions.log.txt"),
    params:
        script=os.path.join(WORKFLOW_DIR, "scripts", "collect_versions.py"),
        software_config=os.environ.get("CHIP_SOFTWARE_CONFIG",
                                       os.path.join(BASE_DIR, "config", "software.yaml")),
        workflow_dir=WORKFLOW_DIR,
        python=os.environ.get("CHIP_PYTHON", "python3"),
    threads: rthreads("software_versions")
    resources:
        mem_mb=rmem("software_versions"),
        runtime_min=rruntime("software_versions"),
        runtime_sec=rruntime_sec("software_versions"),
    shell:
        "{params.python} {params.script} --software-config {params.software_config} "
        "--workflow {params.workflow_dir} --out {output} >> {log} 2>&1"
