# 峰调用：按分组（treat vs control）并行，规则选择由分组的 assay/peak_type 决定
#   chip / cuttag + narrow  -> callpeak_narrow  (MACS2 -q 0.05 + summits)
#   chip / cuttag + broad   -> callpeak_broad   (MACS2 --broad)
#   atac / faire            -> callpeak_atac    (见 config peak.atac.mode)
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
        fmt=lambda wc: "BAMPE" if GROUPS[wc.group]["layout"] == "PE" else "BAM",
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
        fmt=lambda wc: "BAMPE" if GROUPS[wc.group]["layout"] == "PE" else "BAM",
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
        # mode=bampe（默认）: ENCODE ATAC v2 管线做法，按真实片段长度堆叠，
        #   此时 MACS2 忽略 shift/extsize，故不传；
        # mode=shifted: 经典 Tn5 偏移校正配方（-f BAM --nomodel --shift -100 --extsize 200）
        fmt=lambda wc: ("BAMPE" if GROUPS[wc.group]["layout"] == "PE"
                        else "BAM") if config["peak"]["atac"]["mode"] == "bampe" else "BAM",
        nomodel=lambda wc: "" if config["peak"]["atac"]["mode"] == "bampe" else
            "--nomodel --shift {} --extsize {}".format(
                config["peak"]["atac"]["shift"], config["peak"]["atac"]["extsize"]),
        gsize=config["genome_size"],
        keepdup=config["peak"]["keepdup"],
        qvalue=config["peak"]["qvalue"],
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
            {params.nomodel} \
            -B --SPMR --call-summits \
            --outdir 4.peak -n {wildcards.group} > {log} 2>&1
        """


rule bigwig:
    input:
        pileup="4.peak/{group}_treat_pileup.bdg",
        lambda_="4.peak/{group}_control_lambda.bdg",
        chromsize=config["chromsize"],
    output:
        "4.peak/{group}_FE.bw",
    wildcard_constraints:
        group=_group_regex(list(GROUPS)),
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
        bedtools slop -i 4.peak/{wildcards.group}_FE.bdg -g {input.chromsize} -b 0 \
            | bedClip stdin {input.chromsize} 4.peak/{wildcards.group}_FE.clip >> {log} 2>&1
        # bedGraphToBigWig 要求染色体顺序与 chrom.sizes 一致，用 bedtools sort -g 按
        # genome 文件顺序排序（字典序 sort 在 Chr10/Chr2 这类命名下顺序会错）
        bedtools sort -g {input.chromsize} -i 4.peak/{wildcards.group}_FE.clip \
            > 4.peak/{wildcards.group}_FE.clip.sorted
        bedGraphToBigWig 4.peak/{wildcards.group}_FE.clip.sorted {input.chromsize} {output} >> {log} 2>&1
        rm -f 4.peak/{wildcards.group}_FE.bdg \
              4.peak/{wildcards.group}_FE.clip 4.peak/{wildcards.group}_FE.clip.sorted
        """
