# ChIPseeker 峰注释：所有分组的峰文件批量注释 + 分布图

rule peak_annotation:
    input:
        peaks=[group_peak_file(g) for g in GROUPS],
        gtf=config["gtf"],
    output:
        pdf="4.peak/anno_result/Peakanno_PeakDistributions.pdf",
    params:
        script=os.path.join(WORKFLOW_DIR, "scripts", "annoPeak_batch.R"),
        peaklist=lambda wc: ",".join(group_peak_file(g) for g in GROUPS),
        outdir="4.peak/anno_result",
        flank=config["region_flank"],
    threads: 1
    conda:
        os.path.join(ENVS, "r-chipseeker.yaml")
    log:
        "logs/peak_annotation.log",
    shell:
        """
        mkdir -p 4.peak/anno_result logs
        Rscript {params.script} {input.gtf} {params.peaklist} {params.outdir} {params.flank} > {log} 2>&1
        """
