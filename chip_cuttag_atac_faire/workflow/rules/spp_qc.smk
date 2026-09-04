# Optional QC: SPP cross-correlation analysis (NSC / RSC / fragment-length
# estimate), the ENCODE standard. Controlled by config["qc"]["nsc_rsc"];
# recommended off for short-fragment data such as CUT&Tag.
# Thresholds: NSC >= 1.05, RSC >= 0.8 (see the QC table in README).
# Note: run_spp.R -savp plot files derive their names from the input BAM and
# land in the current directory; they are not declared or managed. All
# metrics are contained in the -out text output.

rule spp_crosscorr:
    input:
        bam=lambda wc: sample_bam(wc.sample),
    output:
        txt=R("5.QC/spp/{sample}_spp_crosscorr.txt"),
        fraglen=R("5.QC/spp/{sample}_fragment_len.txt"),
        nsc=R("5.QC/spp/{sample}_NSC.txt"),
        rsc=R("5.QC/spp/{sample}_RSC.txt"),
    log:
        R("logs/spp/{sample}.log"),
    threads: rthreads("spp_crosscorr")
    resources:
        mem_mb=rmem("spp_crosscorr"),
        runtime_min=rruntime("spp_crosscorr"),
        runtime_sec=rruntime_sec("spp_crosscorr"),
    shell:
        """
        run_spp.R -c={input.bam} -p={threads} -out={output.txt} > {log} 2>&1
        sed -r 's/,[^\\t]+//g' {output.txt} | cut -f3 | cut -d, -f1 > {output.fraglen}
        cut -f9 {output.txt} > {output.nsc}
        cut -f10 {output.txt} > {output.rsc}
        """


rule spp_summary:
    input:
        fraglen=expand(R("5.QC/spp/{sample}_fragment_len.txt"), sample=SAMPLES),
        nsc=expand(R("5.QC/spp/{sample}_NSC.txt"), sample=SAMPLES),
        rsc=expand(R("5.QC/spp/{sample}_RSC.txt"), sample=SAMPLES),
    output:
        R("5.QC/spp/NSC_RSC_mqc.tsv"),
    params:
        samples=" ".join(SAMPLES),
        spp_dir=lambda wc, output: os.path.dirname(str(output)),
    log:
        R("logs/spp/summary.log"),
    threads: rthreads("spp_summary")
    resources:
        mem_mb=rmem("spp_summary"),
        runtime_min=rruntime("spp_summary"),
        runtime_sec=rruntime_sec("spp_summary"),
    shell:
        """
        {{
            echo "# id: 'nsc_rsc_table'"
            echo "# section_name: 'SPP cross-correlation (NSC/RSC)'"
            echo "# format: 'tsv'"
            echo "# plot_type: 'table'"
            echo "# pconfig: {{'id': 'nsc_rsc_table', 'title': 'NSC/RSC'}}"
            echo -e "sample\\tfragment_length\\tNSC\\tRSC"
            for s in {params.samples}; do
                printf "%s\\t%s\\t%s\\t%s\\n" "$s" \\
                    "$(cat {params.spp_dir}/${{s}}_fragment_len.txt)" \\
                    "$(cat {params.spp_dir}/${{s}}_NSC.txt)" \\
                    "$(cat {params.spp_dir}/${{s}}_RSC.txt)"
            done
        }} > {output} 2> {log}
        """
