# Repeats/sncRNA pre-filter + unique genome alignment.
# Depends on workflow/rules/common.smk for FILTER_REPEATS, _align_input,
# and resource helpers (the input-routing function lives in common.smk to
# keep this module rules-only).

if FILTER_REPEATS:
    rule star_filter_repeats:
        input:
            fq=R("2.cleandata/{sample}_clean.fqTrTr.sorted.fq.gz"),
            index=R("0.index/repeats_STARindex/SA"),
        output:
            R("3.align/repeats/{sample}_Unmapped.out.mate1"),
        params:
            index=lambda wc, input: os.path.dirname(str(input.index)),
            prefix=lambda wc: f"{RD}3.align/repeats/{wc.sample}_",
            multimap=config["star"]["filter_multimap_nmax"],
        log:
            R("logs/star_filter_repeats/{sample}.log"),
        threads: rthreads("star_filter_repeats")
        resources:
            mem_mb=rmem("star_filter_repeats"),
            runtime_min=rruntime("star_filter_repeats"),
            runtime_sec=rruntime_sec("star_filter_repeats"),
        shell:
            """
            ## STAR spawns FIFO files for threaded input; results/ may live on a
            ## non-FIFO filesystem (e.g. WSL /mnt/* NTFS mounts), so keep the temp
            ## dir on the scheduler-provided tmpdir (Linux /tmp by default).
            ## STAR refuses an existing --outTmpDir — clean, then let it create.
            rm -rf {resources.tmpdir}/STAR_{wildcards.sample}_filter_repeats
            STAR --runThreadN {threads} --runMode alignReads --alignEndsType EndToEnd \
                --genomeDir {params.index} --readFilesCommand zcat \
                --outTmpDir {resources.tmpdir}/STAR_{wildcards.sample}_filter_repeats \
                --genomeLoad NoSharedMemory --outBAMcompression 10 \
                --outFileNamePrefix {params.prefix} \
                --outFilterMultimapNmax {params.multimap} --outFilterMultimapScoreRange 1 \
                --outFilterScoreMin 10 --outFilterType BySJout --outReadsUnmapped Fastx \
                --outSAMattrRGline ID:{wildcards.sample} --outSAMattributes All \
                --outSAMmode Full --outSAMtype BAM Unsorted --outSAMunmapped Within \
                --outStd Log --readFilesIn {input.fq} >> {log} 2>&1
            """


rule star_align:
    input:
        fq=_align_input,
        index=R("0.index/genome_STARindex/SA"),
    output:
        bam=R("3.align/genome/{sample}_Aligned.out.bam"),
        final_log=R("3.align/genome/{sample}_Log.final.out"),
    params:
        index=lambda wc, input: os.path.dirname(str(input.index)),
        prefix=lambda wc: f"{RD}3.align/genome/{wc.sample}_",
        multimap=config["star"]["align_multimap_nmax"],
        # The repeats filter emits a PLAIN (uncompressed) unmapped Fastx;
        # zcat would reject it and STAR would silently see zero reads.
        read_cmd="cat" if FILTER_REPEATS else "zcat",
    log:
        R("logs/star_align/{sample}.log"),
    threads: rthreads("star_align")
    resources:
        mem_mb=rmem("star_align"),
        runtime_min=rruntime("star_align"),
        runtime_sec=rruntime_sec("star_align"),
    shell:
        """
        ## See star_filter_repeats: STAR FIFOs need a Linux tmpdir.
        rm -rf {resources.tmpdir}/STAR_{wildcards.sample}_align
        STAR --runThreadN {threads} --runMode alignReads --alignEndsType EndToEnd \
            --genomeDir {params.index} --readFilesCommand {params.read_cmd} \
            --outTmpDir {resources.tmpdir}/STAR_{wildcards.sample}_align \
            --genomeLoad NoSharedMemory --outBAMcompression 10 \
            --outFileNamePrefix {params.prefix} \
            --outFilterMultimapNmax {params.multimap} --outFilterMultimapScoreRange 1 \
            --outFilterScoreMin 10 --outFilterType BySJout --outReadsUnmapped Fastx \
            --outSAMattrRGline ID:{wildcards.sample} --outSAMattributes All \
            --outSAMmode Full --outSAMtype BAM Unsorted --outSAMunmapped Within \
            --outStd Log --readFilesIn {input.fq} >> {log} 2>&1
        """
