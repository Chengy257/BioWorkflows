# ChIPseeker peak annotation: batch annotation of every group's peak file
# plus the distribution plot.

rule peak_annotation:
    input:
        peaks=[group_peak_file(g) for g in GROUPS],
        gtf=config["gtf"],
    output:
        pdf=R("4.peak/anno_result/Peakanno_PeakDistributions.pdf"),
    params:
        script=os.path.join(WORKFLOW_DIR, "scripts", "annoPeak_batch.R"),
        peaklist=lambda wc: ",".join(group_peak_file(g) for g in GROUPS),
        outdir=lambda wc, output: os.path.dirname(str(output)),
        flank=config["region_flank"],
    threads: rthreads("peak_annotation")
    resources:
        mem_mb=rmem("peak_annotation"),
        runtime_min=rruntime("peak_annotation"),
        runtime_sec=rruntime_sec("peak_annotation"),
    log:
        R("logs/peak_annotation.log"),
    shell:
        """
        Rscript {params.script} {input.gtf} {params.peaklist} {params.outdir} {params.flank} > {log} 2>&1
        """
