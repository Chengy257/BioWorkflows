# Peak calling: one job per group (treat vs control); the concrete rule is
# selected by the group's assay/peak_type:
#   chip / cuttag + narrow  -> callpeak_narrow  (MACS2 -q + summits)
#   chip / cuttag + broad   -> callpeak_broad   (MACS2 --broad)
#   atac / faire            -> callpeak_atac    (see config peak.atac.mode)
# Groups without a control automatically drop -c (MACS2 falls back to a
# local lambda). bigWig signal tracks come from the bigwig rule, built from
# treat_pileup/control_lambda.

rule callpeak_narrow:
    input:
        treat=lambda wc: group_bams(wc.group, "treat"),
        control=lambda wc: group_bams(wc.group, "control"),
    output:
        peaks=R("4.peak/{group}_peaks.narrowPeak"),
        summits=R("4.peak/{group}_summits.bed"),
        pileup=R("4.peak/{group}_treat_pileup.bdg"),
        lambda_=R("4.peak/{group}_control_lambda.bdg"),
    wildcard_constraints:
        group=_group_regex(_groups_of("chip", "narrow") + _groups_of("cuttag", "narrow")),
    params:
        control=group_control_arg,
        fmt=lambda wc: "BAMPE" if GROUPS[wc.group]["layout"] == "PE" else "BAM",
        gsize=config["genome_size"],
        keepdup=config["peak"]["keepdup"],
        qvalue=config["peak"]["qvalue"],
        outdir=lambda wc, output: os.path.dirname(output.peaks),
    log:
        R("logs/callpeak/{group}_narrow.log"),
    threads: rthreads("callpeak_narrow")  # MACS2 is single-threaded; request accordingly
    resources:
        mem_mb=rmem("callpeak_narrow"),
        runtime_min=rruntime("callpeak_narrow"),
        runtime_sec=rruntime_sec("callpeak_narrow"),
    shell:
        """
        macs2 callpeak \
            -t {input.treat} {params.control} \
            -f {params.fmt} -g {params.gsize} \
            --keep-dup {params.keepdup} -q {params.qvalue} \
            -B --SPMR --call-summits \
            --outdir {params.outdir} -n {wildcards.group} > {log} 2>&1
        """


rule callpeak_broad:
    input:
        treat=lambda wc: group_bams(wc.group, "treat"),
        control=lambda wc: group_bams(wc.group, "control"),
    output:
        peaks=R("4.peak/{group}_peaks.broadPeak"),
        pileup=R("4.peak/{group}_treat_pileup.bdg"),
        lambda_=R("4.peak/{group}_control_lambda.bdg"),
    wildcard_constraints:
        group=_group_regex(_groups_of("chip", "broad") + _groups_of("cuttag", "broad")),
    params:
        control=group_control_arg,
        fmt=lambda wc: "BAMPE" if GROUPS[wc.group]["layout"] == "PE" else "BAM",
        gsize=config["genome_size"],
        keepdup=config["peak"]["keepdup"],
        broad_cutoff=config["peak"]["broad_cutoff"],
        outdir=lambda wc, output: os.path.dirname(output.peaks),
    log:
        R("logs/callpeak/{group}_broad.log"),
    threads: rthreads("callpeak_broad")  # MACS2 is single-threaded; request accordingly
    resources:
        mem_mb=rmem("callpeak_broad"),
        runtime_min=rruntime("callpeak_broad"),
        runtime_sec=rruntime_sec("callpeak_broad"),
    shell:
        """
        macs2 callpeak \
            -t {input.treat} {params.control} \
            -f {params.fmt} -g {params.gsize} \
            --keep-dup {params.keepdup} --broad --broad-cutoff {params.broad_cutoff} \
            -B --SPMR \
            --outdir {params.outdir} -n {wildcards.group} > {log} 2>&1
        """


rule callpeak_atac:
    input:
        treat=lambda wc: group_bams(wc.group, "treat"),
        control=lambda wc: group_bams(wc.group, "control"),
    output:
        peaks=R("4.peak/{group}_peaks.narrowPeak"),
        summits=R("4.peak/{group}_summits.bed"),
        pileup=R("4.peak/{group}_treat_pileup.bdg"),
        lambda_=R("4.peak/{group}_control_lambda.bdg"),
    wildcard_constraints:
        group=_group_regex(_groups_of("atac") + _groups_of("faire")),
    params:
        control=group_control_arg,
        # mode=bampe (default): the ENCODE ATAC v2 recipe piling up real
        #   fragment lengths; MACS2 then ignores shift/extsize, so they are
        #   not passed at all;
        # mode=shifted: the classic Tn5 offset correction recipe
        #   (-f BAM --nomodel --shift -100 --extsize 200)
        fmt=lambda wc: ("BAMPE" if GROUPS[wc.group]["layout"] == "PE"
                        else "BAM") if config["peak"]["atac"]["mode"] == "bampe" else "BAM",
        nomodel=lambda wc: "" if config["peak"]["atac"]["mode"] == "bampe" else
            "--nomodel --shift {} --extsize {}".format(
                config["peak"]["atac"]["shift"], config["peak"]["atac"]["extsize"]),
        gsize=config["genome_size"],
        keepdup=config["peak"]["keepdup"],
        qvalue=config["peak"]["qvalue"],
        outdir=lambda wc, output: os.path.dirname(output.peaks),
    log:
        R("logs/callpeak/{group}_atac_faire.log"),
    threads: rthreads("callpeak_atac")  # MACS2 is single-threaded; request accordingly
    resources:
        mem_mb=rmem("callpeak_atac"),
        runtime_min=rruntime("callpeak_atac"),
        runtime_sec=rruntime_sec("callpeak_atac"),
    shell:
        """
        macs2 callpeak \
            -t {input.treat} {params.control} \
            -f {params.fmt} -g {params.gsize} \
            --keep-dup {params.keepdup} -q {params.qvalue} \
            {params.nomodel} \
            -B --SPMR --call-summits \
            --outdir {params.outdir} -n {wildcards.group} > {log} 2>&1
        """


rule bigwig:
    input:
        pileup=R("4.peak/{group}_treat_pileup.bdg"),
        lambda_=R("4.peak/{group}_control_lambda.bdg"),
        chromsize=config["chromsize"],
    output:
        R("4.peak/{group}_FE.bw"),
    wildcard_constraints:
        group=_group_regex(list(GROUPS)),
    params:
        measure=PEAK_MEASURE,   # FE (fold enrichment, default) | logFE
        outdir=lambda wc, output: os.path.dirname(str(output)),
    log:
        R("logs/bigwig/{group}.log"),
    threads: rthreads("bigwig")
    resources:
        mem_mb=rmem("bigwig"),
        runtime_min=rruntime("bigwig"),
        runtime_sec=rruntime_sec("bigwig"),
    shell:
        """
        macs2 bdgcmp \
            -t {input.pileup} -c {input.lambda_} \
            -o {params.outdir}/{wildcards.group}_FE.bdg -m {params.measure} -p 0.00001 > {log} 2>&1
        bedtools slop -i {params.outdir}/{wildcards.group}_FE.bdg -g {input.chromsize} -b 0 \
            | bedClip stdin {input.chromsize} {params.outdir}/{wildcards.group}_FE.clip >> {log} 2>&1
        # bedGraphToBigWig validates C-collation order (Chr1 < Chr10 < Chr11 <
        # Chr12 < Chr2) plus numeric starts, regardless of the chrom.sizes line
        # order; bedtools sort -g follows the chrom.sizes line order instead
        # and fails the check whenever that order differs from C collation.
        LC_COLLATE=C sort -k1,1 -k2,2n {params.outdir}/{wildcards.group}_FE.clip \
            > {params.outdir}/{wildcards.group}_FE.clip.sorted
        bedGraphToBigWig {params.outdir}/{wildcards.group}_FE.clip.sorted {input.chromsize} {output} >> {log} 2>&1
        rm -f {params.outdir}/{wildcards.group}_FE.bdg \
              {params.outdir}/{wildcards.group}_FE.clip {params.outdir}/{wildcards.group}_FE.clip.sorted
        """


rule bigwig_sample:
    # Per-sample normalized coverage track (v0.5; requested only when
    # bigwig.per_sample is true): browser-comparable replicate tracks built
    # from the same analysis BAM every other signal consumer uses. The group
    # FE track filename stays {group}_FE.bw in logFE mode (the measure config
    # changes the track values, not the layout).
    input:
        lambda wc: sample_bam(wc.sample),
    output:
        R("4.peak/samples/{sample}.bw"),
    wildcard_constraints:
        sample=_group_regex(SAMPLES),
    params:
        normalize=BIGWIG["normalize"],
        binsize=BIGWIG["bin"],
        # RPGC needs an integer effective genome size; MACS2-style strings
        # ("3.7e8") convert through float
        gsize=int(float(str(config["genome_size"]))),
    log:
        R("logs/bigwig/{sample}_sample.log"),
    threads: rthreads("bigwig_sample")
    resources:
        mem_mb=rmem("bigwig_sample"),
        runtime_min=rruntime("bigwig_sample"),
        runtime_sec=rruntime_sec("bigwig_sample"),
    shell:
        """
        bamCoverage -b {input} --normalizeUsing {params.normalize} \
            --effectiveGenomeSize {params.gsize} --binSize {params.binsize} \
            -p {threads} -o {output} > {log} 2>&1
        """
