# ChIPseeker peak annotation: batch annotation of every group's FINAL peak
# set (pooled by default; the IDR/consensus set when the replicate stage is
# enabled; blacklist-filtered copies when a blacklist is configured) plus
# the distribution plot.

rule peak_annotation:
    input:
        peaks=annot_peak_files(),
        gtf=config["gtf"],
    output:
        pdf=R("4.peak/anno_result/Peakanno_PeakDistributions.pdf"),
    params:
        script=os.path.join(WORKFLOW_DIR, "scripts", "annoPeak_batch.R"),
        rscript=RSCRIPT,
        peaklist=lambda wc: ",".join(annot_peak_files()),
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
        {params.rscript} {params.script} {input.gtf} {params.peaklist} {params.outdir} {params.flank} > {log} 2>&1
        """
