# 可选 QC：SPP 交叉相关分析（NSC / RSC / 片段长度估计），ENCODE 标准
# 由 config["qc"]["nsc_rsc"] 开关控制；CUT&Tag 等短片段数据建议关闭。
# 阈值参考：NSC >= 1.05、RSC >= 0.8（详见 README QC 表）。
# 注：run_spp.R 的 -savp 图件输出文件名依输入 BAM 派生且落在当前目录，
# 不做声明管理，故不启用；全部指标已含于 -out 文本输出。

rule spp_crosscorr:
    input:
        bam=lambda wc: sample_bam(wc.sample),
    output:
        txt="5.QC/spp/{sample}_spp_crosscorr.txt",
        fraglen="5.QC/spp/{sample}_fragment_len.txt",
        nsc="5.QC/spp/{sample}_NSC.txt",
        rsc="5.QC/spp/{sample}_RSC.txt",
    log:
        "logs/spp/{sample}.log",
    threads: config["threads"]
    resources:
        mem_mb=res("spp_crosscorr", 8192),
        runtime_min=res("spp_crosscorr", 180),
        runtime_sec=res("spp_crosscorr", 180) * 60,
    shell:
        """
        mkdir -p 5.QC/spp logs/spp
        run_spp.R -c={input.bam} -p={threads} -out={output.txt} > {log} 2>&1
        sed -r 's/,[^\\t]+//g' {output.txt} | cut -f3 | cut -d, -f1 > {output.fraglen}
        cut -f9 {output.txt} > {output.nsc}
        cut -f10 {output.txt} > {output.rsc}
        """


rule spp_summary:
    input:
        fraglen=expand("5.QC/spp/{sample}_fragment_len.txt", sample=SAMPLES),
        nsc=expand("5.QC/spp/{sample}_NSC.txt", sample=SAMPLES),
        rsc=expand("5.QC/spp/{sample}_RSC.txt", sample=SAMPLES),
    output:
        "5.QC/spp/NSC_RSC_mqc.tsv",
    params:
        samples=" ".join(SAMPLES),
    log:
        "logs/spp/summary.log",
    threads: 1
    resources:
        mem_mb=res("spp_summary", 1024),
        runtime_min=res("spp_summary", 10),
        runtime_sec=res("spp_summary", 10) * 60,
    shell:
        """
        mkdir -p 5.QC/spp logs/spp
        {{
            echo "# id: 'nsc_rsc_table'"
            echo "# section_name: 'SPP cross-correlation (NSC/RSC)'"
            echo "# format: 'tsv'"
            echo "# plot_type: 'table'"
            echo "# pconfig: {{'id': 'nsc_rsc_table', 'title': 'NSC/RSC'}}"
            echo -e "sample\\tfragment_length\\tNSC\\tRSC"
            for s in {params.samples}; do
                printf "%s\\t%s\\t%s\\t%s\\n" "$s" \\
                    "$(cat 5.QC/spp/${{s}}_fragment_len.txt)" \\
                    "$(cat 5.QC/spp/${{s}}_NSC.txt)" \\
                    "$(cat 5.QC/spp/${{s}}_RSC.txt)"
            done
        }} > {output} 2> {log}
        """
