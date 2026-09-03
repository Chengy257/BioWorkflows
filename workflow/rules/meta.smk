# 运行元数据：实际软件版本记录。
rule software_versions:
    output:
        "5.QC/software_versions.yaml",
    log:
        "5.QC/logs/software_versions.log.txt",
    params:
        script=os.path.join(WORKFLOW_DIR, "scripts", "collect_versions.py"),
        software_config=os.environ.get("CHIP_SOFTWARE_CONFIG",
                                       os.path.join(BASE_DIR, "config", "software.yaml")),
        workflow_dir=WORKFLOW_DIR,
    threads: 1
    resources:
        mem_mb=1024,
        runtime_min=10,
        runtime_sec=600,
    shell:
        "python3 {params.script} --software-config {params.software_config} "
        "--workflow {params.workflow_dir} --out {output} >> {log} 2>&1"
