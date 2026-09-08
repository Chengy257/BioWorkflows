# UMI deduplication, read counts, and peak calling.
# Depends on workflow/rules/common.smk for CLIPPER/CLIPPER_ENABLED/helpers.
# Target-set note: the rule definitions here stay sample-agnostic (upstream
# trim/align/dedup run for every declared sample). Which samples join the
# PureCLIP target set is decided in common.smk (PEAKCALL_SAMPLES): with the
# condition/role sample table and reproducible_peaks enabled, the ip-role
# samples only — unless reproducible_peaks.input_control is enabled, in which
# case role=input samples run PureCLIP too so their beds can form the
# per-condition background consensus (rules/consensus.smk).

rule umi_dedup:
    input:
        bam=R("3.align/genome/{sample}_Aligned.out.bam"),
    output:
        bam=R("4.rmdup/{sample}.rmDupSo.bam"),
        bai=R("4.rmdup/{sample}.rmDupSo.bam.bai"),
        stats=R("4.rmdup/{sample}_stats/{sample}_edit_distance.tsv"),
    params:
        tmp_bam=lambda wc: f"{RD}4.rmdup/{wc.sample}.rmDup.bam",
        sorted_bam=lambda wc: f"{RD}4.rmdup/{wc.sample}.Aligned.sorted.bam",
        stats_dir=lambda wc: f"{RD}4.rmdup/{wc.sample}_stats",
        stats_prefix=lambda wc: f"{RD}4.rmdup/{wc.sample}_stats/{wc.sample}",
    log:
        R("logs/umi_dedup/{sample}.log"),
    threads: rthreads("umi_dedup")
    resources:
        mem_mb=rmem("umi_dedup"),
        runtime_min=rruntime("umi_dedup"),
        runtime_sec=rruntime_sec("umi_dedup"),
    shell:
        """
        ## umi_tools dedup requires a coordinate-sorted, INDEXED input BAM
        ## (pysam fetch); STAR writes unsorted BAMs, so sort+index first.
        ## Dedup preserves input order, so its output is already sorted.
        mkdir -p {params.stats_dir}
        samtools sort -@ {threads} -o {params.sorted_bam} {input.bam} >> {log} 2>&1
        samtools index -@ {threads} {params.sorted_bam} >> {log} 2>&1
        umi_tools dedup --random-seed 1 -I {params.sorted_bam} --method unique \
            --output-stats {params.stats_prefix} -S {output.bam} >> {log} 2>&1
        samtools index -@ {threads} {output.bam} >> {log} 2>&1
        rm -f {params.sorted_bam} {params.sorted_bam}.bai
        """


rule read_count:
    input:
        R("4.rmdup/{sample}.rmDupSo.bam"),
    output:
        R("4.rmdup/{sample}_readnum.txt"),
    log:
        R("logs/read_count/{sample}.log"),
    threads: rthreads("read_count")
    resources:
        mem_mb=rmem("read_count"),
        runtime_min=rruntime("read_count"),
        runtime_sec=rruntime_sec("read_count"),
    shell:
        """
        samtools view -c -F 4 {input} > {output} 2> {log}
        """


if CLIPPER_ENABLED:
    rule callpeak_clipper:
        input:
            bam=R("4.rmdup/{sample}.rmDupSo.bam"),
        output:
            R("5.callpeak/{sample}.clipper.peakClusters.bed"),
        params:
            exe=CLIPPER,
            species=config["callpeak"]["clipper_species"],
        log:
            R("logs/callpeak_clipper/{sample}.log"),
        threads: rthreads("callpeak_clipper")
        resources:
            mem_mb=rmem("callpeak_clipper"),
            runtime_min=rruntime("callpeak_clipper"),
            runtime_sec=rruntime_sec("callpeak_clipper"),
        shell:
            """
            {params.exe} --species {params.species} --bam {input.bam} \
                --outfile {output} > {log} 2>&1
            rm -f {output}.tsv
            """


rule callpeak_pureclip:
    input:
        bam=R("4.rmdup/{sample}.rmDupSo.bam"),
        bai=R("4.rmdup/{sample}.rmDupSo.bam.bai"),
        genome=config["genome"],
    output:
        R("5.callpeak/{sample}.pureclip.bed"),
    log:
        R("logs/callpeak_pureclip/{sample}.log"),
    threads: rthreads("callpeak_pureclip")
    resources:
        mem_mb=rmem("callpeak_pureclip"),
        runtime_min=rruntime("callpeak_pureclip"),
        runtime_sec=rruntime_sec("callpeak_pureclip"),
    shell:
        """
        pureclip -i {input.bam} -bai {input.bai} -g {input.genome} \
            -nt {threads} -o {output} > {log} 2>&1
        """
