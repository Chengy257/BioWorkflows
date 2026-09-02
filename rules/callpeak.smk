# 峰调用：按分组（treat vs control）并行，规则选择由分组的 assay/peak_type 决定
#   chip / cuttag + narrow  -> callpeak_narrow  (MACS2 -q 0.05 + summits)
#   chip / cuttag + broad   -> callpeak_broad   (MACS2 --broad)
#   atac / faire            -> callpeak_atac    (MACS2 --nomodel --shift -100 --extsize 200)
# 分组无对照时自动省略 -c（MACS2 回退到局部 lambda）。
# bigWig 信号轨道由 bigwig 规则从 treat_pileup/control_lambda 生成。

rule callpeak_narrow:
    input:
        treat=lambda wc: group_bams(wc.group, "treat"),
        control=lambda wc: group_bams(wc.group, "control"),
    output:
        peaks="4.peak/{group}_peaks.narrowPeak",
        summits="4.peak/{group}_summits.bed",
        pileup="4.peak/{group}_treat_pileup.bdg",
        lambda_="4.peak/{group}_control_lambda.bdg",
    wildcard_constraints:
        group=_group_regex(_groups_of("chip", "narrow") + _groups_of("cuttag", "narrow")),
    params:
        control=group_control_arg,
        fmt=lambda wc: "BAMPE" if GROUPS[wc.group]["layout"] == "PE" else "BAMSOLO",
        gsize=config["genome_size"],
        keepdup=config["peak"]["keepdup"],
        qvalue=config["peak"]["qvalue"],
    log:
        "logs/callpeak/{group}_narrow.log",
    threads: config["threads"]
    conda:
        os.path.join(ENVS, "macs2.yaml")
    shell:
        """
        mkdir -p 4.peak logs/callpeak
        macs2 callpeak \
            -t {input.treat} {params.control} \
            -f {params.fmt} -g {params.gsize} \
            --keep-dup {params.keepdup} -q {params.qvalue} \
            -B --SPMR --call-summits \
            --outdir 4.peak -n {wildcards.group} > {log} 2>&1
        """


rule callpeak_broad:
    input:
        treat=lambda wc: group_bams(wc.group, "treat"),
        control=lambda wc: group_bams(wc.group, "control"),
    output:
        peaks="4.peak/{group}_peaks.broadPeak",
        pileup="4.peak/{group}_treat_pileup.bdg",
        lambda_="4.peak/{group}_control_lambda.bdg",
    wildcard_constraints:
        group=_group_regex(_groups_of("chip", "broad") + _groups_of("cuttag", "broad")),
    params:
        control=group_control_arg,
        fmt=lambda wc: "BAMPE" if GROUPS[wc.group]["layout"] == "PE" else "BAMSOLO",
        gsize=config["genome_size"],
        keepdup=config["peak"]["keepdup"],
        broad_cutoff=config["peak"]["broad_cutoff"],
    log:
        "logs/callpeak/{group}_broad.log",
    threads: config["threads"]
    conda:
        os.path.join(ENVS, "macs2.yaml")
    shell:
        """
        mkdir -p 4.peak logs/callpeak
        macs2 callpeak \
            -t {input.treat} {params.control} \
            -f {params.fmt} -g {params.gsize} \
            --keep-dup {params.keepdup} --broad --broad-cutoff {params.broad_cutoff} \
            -B --SPMR \
            --outdir 4.peak -n {wildcards.group} > {log} 2>&1
        """


rule callpeak_atac:
    input:
        treat=lambda wc: group_bams(wc.group, "treat"),
        control=lambda wc: group_bams(wc.group, "control"),
    output:
        peaks="4.peak/{group}_peaks.narrowPeak",
        summits="4.peak/{group}_summits.bed",
        pileup="4.peak/{group}_treat_pileup.bdg",
        lambda_="4.peak/{group}_control_lambda.bdg",
    wildcard_constraints:
        group=_group_regex(_groups_of("atac") + _groups_of("faire")),
    params:
        control=group_control_arg,
        fmt=lambda wc: "BAMPE" if GROUPS[wc.group]["layout"] == "PE" else "BAMSOLO",
        gsize=config["genome_size"],
        keepdup=config["peak"]["keepdup"],
        qvalue=config["peak"]["qvalue"],
        nomodel_shift=config["peak"]["atac"]["shift"],
        nomodel_extsize=config["peak"]["atac"]["extsize"],
    log:
        "logs/callpeak/{group}_atac_faire.log",
    threads: config["threads"]
    conda:
        os.path.join(ENVS, "macs2.yaml")
    shell:
        """
        mkdir -p 4.peak logs/callpeak
        macs2 callpeak \
            -t {input.treat} {params.control} \
            -f {params.fmt} -g {params.gsize} \
            --keep-dup {params.keepdup} -q {params.qvalue} \
            --nomodel --shift {params.nomodel_shift} --extsize {params.nomodel_extsize} \
            -B --SPMR --call-summits \
            --outdir 4.peak -n {wildcards.group} > {log} 2>&1
        """


rule bigwig:
    input:
        pileup="4.peak/{group}_treat_pileup.bdg",
        lambda_="4.peak/{group}_control_lambda.bdg",
    output:
        "4.peak/{group}_FE.bw",
    wildcard_constraints:
        group=_group_regex(list(GROUPS)),
    params:
        chromsize=config["chromsize"],
    log:
        "logs/bigwig/{group}.log",
    threads: 1
    conda:
        os.path.join(ENVS, "bigwig.yaml")
    shell:
        """
        mkdir -p 4.peak logs/bigwig
        macs2 bdgcmp \
            -t {input.pileup} -c {input.lambda_} \
            -o 4.peak/{wildcards.group}_FE.bdg -m FE -p 0.00001 > {log} 2>&1
        bedtools slop -i 4.peak/{wildcards.group}_FE.bdg -g {params.chromsize} -b 0 \
            | bedClip stdin {params.chromsize} 4.peak/{wildcards.group}_FE.clip >> {log} 2>&1
        LC_COLLATE=C sort -k1,1 -k2,2n 4.peak/{wildcards.group}_FE.clip \
            > 4.peak/{wildcards.group}_FE.clip.sorted
        bedGraphToBigWig 4.peak/{wildcards.group}_FE.clip.sorted {params.chromsize} {output} >> {log} 2>&1
        rm -f 4.peak/{wildcards.group}_FE.bdg \
              4.peak/{wildcards.group}_FE.clip 4.peak/{wildcards.group}_FE.clip.sorted
        """
