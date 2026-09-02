# deeptools QC：样本相关性热图 / PCA / 指纹图 / 片段长分布 / 基因区信号 profile
# 输入为各样本的"分析用 BAM"（去重 assay 用 rmdup，否则 sorted），信号轨道用分组 FE bigWig。

ANALYSIS_BAMS = [sample_bam(s) for s in SAMPLES]

rule deeptools_multibamsummary:
    input:
        ANALYSIS_BAMS,
    output:
        npz="5.QC_deeptools/multiBamSummary.npz",
        counts="5.QC_deeptools/readCounts.tab",
    params:
        labels=" ".join(SAMPLES),
    log:
        "logs/deeptools/multiBamSummary.log",
    threads: config["threads"]
    conda:
        os.path.join(ENVS, "deeptools.yaml")
    shell:
        """
        mkdir -p 5.QC_deeptools logs/deeptools
        multiBamSummary bins -p {threads} --bamfiles {input} \
            --minMappingQuality {config[min_mapq]} --labels {params.labels} \
            -out {output.npz} --outRawCounts {output.counts} > {log} 2>&1
        """


rule deeptools_correlation:
    input:
        "5.QC_deeptools/multiBamSummary.npz",
    output:
        png="5.QC_deeptools/heatmap_SpearmanCorr_readCounts.png",
        tab="5.QC_deeptools/SpearmanCorr_readCounts.tab",
    log:
        "logs/deeptools/correlation.log",
    conda:
        os.path.join(ENVS, "deeptools.yaml")
    shell:
        """
        plotCorrelation -in {input} --corMethod spearman --skipZeros \
            --plotTitle "Spearman Correlation of Read Counts" \
            --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
            -o {output.png} --outFileCorMatrix {output.tab} > {log} 2>&1
        """


rule deeptools_pca:
    input:
        "5.QC_deeptools/multiBamSummary.npz",
    output:
        "5.QC_deeptools/PCA_readCounts.png",
    log:
        "logs/deeptools/pca.log",
    conda:
        os.path.join(ENVS, "deeptools.yaml")
    shell:
        """
        plotPCA -in {input} -o {output} -T "PCA of read counts" > {log} 2>&1
        """


rule deeptools_fingerprint:
    input:
        ANALYSIS_BAMS,
    output:
        png="5.QC_deeptools/fingerprints.png",
        tab="5.QC_deeptools/fingerprints.tab",
    params:
        labels=" ".join(SAMPLES),
    log:
        "logs/deeptools/fingerprint.log",
    threads: config["threads"]
    conda:
        os.path.join(ENVS, "deeptools.yaml")
    shell:
        """
        plotFingerprint -b {input} --labels {params.labels} \
            --minMappingQuality {config[min_mapq]} --skipZeros \
            -T "Fingerprints of different samples" \
            --plotFile {output.png} --outRawCounts {output.tab} > {log} 2>&1
        """


rule deeptools_fragmentsize:
    input:
        ANALYSIS_BAMS,
    output:
        "5.QC_deeptools/fragmentsize.png",
    params:
        labels=" ".join(SAMPLES),
    log:
        "logs/deeptools/fragmentsize.log",
    threads: config["threads"]
    conda:
        os.path.join(ENVS, "deeptools.yaml")
    shell:
        """
        bamPEFragmentSize -p {threads} -hist {output} \
            -T "Fragment size of PE data" --maxFragmentLength 1000 \
            -b {input} --samplesLabel {params.labels} > {log} 2>&1
        """


rule deeptools_profile:
    input:
        bed=config["bed"],
        bws=expand("4.peak/{group}_FE.bw", group=GROUPS),
    output:
        matrix="5.QC_deeptools/matrix_scaled.gz",
        profile="5.QC_deeptools/profile_scaled.png",
        heatmap="5.QC_deeptools/heatmap_scaled.png",
    log:
        "logs/deeptools/profile.log",
    threads: config["threads"]
    conda:
        os.path.join(ENVS, "deeptools.yaml")
    shell:
        """
        computeMatrix scale-regions -R {input.bed} -S {input.bws} \
            -b 3000 -a 3000 -o {output.matrix} -p {threads} > {log} 2>&1
        plotProfile -m {output.matrix} -out {output.profile} \
            --plotTitle "Mean signal over genes" >> {log} 2>&1
        plotHeatmap -m {output.matrix} -out {output.heatmap} >> {log} 2>&1
        """
