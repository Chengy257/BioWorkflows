# 可选 QC：SPP 交叉相关分析（NSC / RSC / 片段长度估计），ENCODE 标准
# 由 config["qc"]["nsc_rsc"] 开关控制；CUT&Tag 等短片段数据建议关闭。
# 阈值参考：NSC >= 1.05、RSC >= 0.8（详见 README QC 表）。

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
    conda:
        os.path.join(ENVS, "phantompeakqualtools.yaml")
    shell:
        """
        mkdir -p 5.QC/spp logs/spp
        run_spp.R -c={input.bam} -p={threads} -savp -out={output.txt} > {log} 2>&1
        sed -r 's/,[^\\t]+//g' {output.txt} | cut -f3 | cut -d, -f1 > {output.fraglen}
        cut -f9 {output.txt} > {output.nsc}
        cut -f10 {output.txt} > {output.rsc}
        """
