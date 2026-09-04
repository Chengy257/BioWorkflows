# deeptools QC: sample-correlation heatmap / PCA / fingerprints / fragment
# size distribution / gene-region signal profiles.
# Inputs are each sample's "analysis BAM" (rmdup for dedup assays, otherwise
# sorted); signal tracks use the per-group FE bigWigs.

ANALYSIS_BAMS = [sample_bam(s) for s in SAMPLES]

rule deeptools_multibamsummary:
    input:
        ANALYSIS_BAMS,
    output:
        npz=R("5.QC/deeptools/multiBamSummary.npz"),
        counts=R("5.QC/deeptools/readCounts.tab"),
    params:
        labels=" ".join(SAMPLES),
    log:
        R("logs/deeptools/multiBamSummary.log"),
    threads: rthreads("deeptools_multibamsummary")
    resources:
        mem_mb=rmem("deeptools_multibamsummary"),
        runtime_min=rruntime("deeptools_multibamsummary"),
        runtime_sec=rruntime_sec("deeptools_multibamsummary"),
    shell:
        """
        multiBamSummary bins -p {threads} --bamfiles {input} \
            --minMappingQuality {config[min_mapq]} --labels {params.labels} \
            -out {output.npz} --outRawCounts {output.counts} > {log} 2>&1
        """


rule deeptools_correlation:
    input:
        R("5.QC/deeptools/multiBamSummary.npz"),
    output:
        png=R("5.QC/deeptools/heatmap_SpearmanCorr_readCounts.png"),
        tab=R("5.QC/deeptools/SpearmanCorr_readCounts.tab"),
    log:
        R("logs/deeptools/correlation.log"),
    threads: rthreads("deeptools_correlation")
    resources:
        mem_mb=rmem("deeptools_correlation"),
        runtime_min=rruntime("deeptools_correlation"),
        runtime_sec=rruntime_sec("deeptools_correlation"),
    shell:
        """
        plotCorrelation -in {input} --corMethod spearman --skipZeros \
            --plotTitle "Spearman Correlation of Read Counts" \
            --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
            -o {output.png} --outFileCorMatrix {output.tab} > {log} 2>&1
        """


rule deeptools_pca:
    input:
        R("5.QC/deeptools/multiBamSummary.npz"),
    output:
        R("5.QC/deeptools/PCA_readCounts.png"),
    log:
        R("logs/deeptools/pca.log"),
    threads: rthreads("deeptools_pca")
    resources:
        mem_mb=rmem("deeptools_pca"),
        runtime_min=rruntime("deeptools_pca"),
        runtime_sec=rruntime_sec("deeptools_pca"),
    shell:
        """
        plotPCA -in {input} -o {output} -T "PCA of read counts" > {log} 2>&1
        """


rule deeptools_fingerprint:
    input:
        ANALYSIS_BAMS,
    output:
        png=R("5.QC/deeptools/fingerprints.png"),
        tab=R("5.QC/deeptools/fingerprints.tab"),
    params:
        labels=" ".join(SAMPLES),
    log:
        R("logs/deeptools/fingerprint.log"),
    threads: rthreads("deeptools_fingerprint")
    resources:
        mem_mb=rmem("deeptools_fingerprint"),
        runtime_min=rruntime("deeptools_fingerprint"),
        runtime_sec=rruntime_sec("deeptools_fingerprint"),
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
        R("5.QC/deeptools/fragmentsize.png"),
    params:
        labels=" ".join(SAMPLES),
    log:
        R("logs/deeptools/fragmentsize.log"),
    threads: rthreads("deeptools_fragmentsize")
    resources:
        mem_mb=rmem("deeptools_fragmentsize"),
        runtime_min=rruntime("deeptools_fragmentsize"),
        runtime_sec=rruntime_sec("deeptools_fragmentsize"),
    shell:
        """
        bamPEFragmentSize -p {threads} -hist {output} \
            -T "Fragment size of PE data" --maxFragmentLength 1000 \
            -b {input} --samplesLabel {params.labels} > {log} 2>&1
        """


rule deeptools_profile:
    input:
        bed=config["bed"],
        bws=expand(R("4.peak/{group}_FE.bw"), group=GROUPS),
    output:
        matrix=R("5.QC/deeptools/matrix_scaled.gz"),
        profile=R("5.QC/deeptools/profile_scaled.png"),
        heatmap=R("5.QC/deeptools/heatmap_scaled.png"),
    log:
        R("logs/deeptools/profile.log"),
    threads: rthreads("deeptools_profile")
    resources:
        mem_mb=rmem("deeptools_profile"),
        runtime_min=rruntime("deeptools_profile"),
        runtime_sec=rruntime_sec("deeptools_profile"),
    params:
        flank=config["region_flank"],
    shell:
        """
        computeMatrix scale-regions -R {input.bed} -S {input.bws} \
            -b {params.flank} -a {params.flank} -o {output.matrix} -p {threads} > {log} 2>&1
        plotProfile -m {output.matrix} -out {output.profile} \
            --plotTitle "Mean signal over genes" >> {log} 2>&1
        plotHeatmap -m {output.matrix} -out {output.heatmap} >> {log} 2>&1
        """
