# rules/align.smk — bismark alignment, deduplication, nucleotide stats.
#
# Alignment core of the bs-seq pipeline: Bismark (bowtie2 backend) aligns
# trim-aware reads against the bisulfite genome prepared in rules/index.smk
# (the GA-conversion bowtie2 index is the sentinel input; the genome folder
# is three dirname levels above it), then deduplicate_bismark removes PCR
# duplicates and bam2nuc collects nucleotide composition — genome-wide once
# and per sample on the dedup BAM (the per-sample definition therefore also
# blocks on the genome-wide totals).
#
# Derived-name convention these definitions emit (keep in sync with TARGETS
# in common.smk and the C5 methylation module):
# - bismark --basename {sample} writes {sample}.bam plus
#   {sample}_report.txt into the --od folder (3.align/);
# - deduplicate_bismark inherits the full alignment basename and writes
#   into --od, so the dedup BAM is {sample}.bam.deduplicated.bam; the
#   declared dedup report is {sample}.dedup_report.txt (4.dedup/) while
#   the tool may itself write {sample}.bam.dedup_report.txt (Phase D
#   reconciliation);
# - bam2nuc names the per-sample output after the input BAM; the declared
#   stats path is {sample}.nucleotide_stats.txt under
#   5.methylation/{sample}/.
#
# Phase D real-run validation against the pinned bismark 0.24.0 (tracked
# in docs/TODO.md) must confirm: bismark accepts --basename in the pinned
# version (fallback = post-run rename, documented in docs/TODO.md), the
# legacy "--pe"/-p paired flag, the exact dedup report name, and the exact
# bam2nuc output name — adjust the declared names if the tools differ.
# Keep --phred33-quals in bismark.align_extra.

rule bismark_align:
    input:
        fq=lambda wc: align_input(wc.sample),
        genome=R("0.index/bismark_genome/Bisulfite_Genome/GA_conversion/BS_conv.1.bt2.lf"),
    output:
        bam=R("3.align/{sample}.bam"),
        report=R("3.align/{sample}_report.txt"),
    params:
        genome_dir=lambda wc, input: os.path.dirname(os.path.dirname(
            os.path.dirname(str(input.genome)))),
        outdir=lambda wc, output: os.path.dirname(str(output.bam)),
        reads=lambda wc: align_input_arg(wc.sample),
        pe_flag=lambda wc: "--pe" if layout_of(wc.sample) == "PE" else "",
        extra=config["bismark"]["align_extra"],
    log: R("logs/bismark_align/{sample}.log"),
    threads: rthreads("bismark_align")
    resources:
        mem_mb=rmem("bismark_align"), runtime_min=rruntime("bismark_align"),
        runtime_sec=rruntime_sec("bismark_align"),
    shell:
        """
        bismark --parallel {threads} --genome {params.genome_dir} \
            {params.extra} {params.pe_flag} --basename {wildcards.sample} \
            --bam --od {params.outdir} {params.reads} > {log} 2>&1
        """


rule deduplicate:
    input:
        bam=R("3.align/{sample}.bam"),
    output:
        # deduplicate_bismark derives the output name from the input file
        # name and writes it into --od: {sample}.bam -> {sample}.bam.deduplicated.bam
        bam=R("4.dedup/{sample}.bam.deduplicated.bam"),
        report=R("4.dedup/{sample}.dedup_report.txt"),
    params:
        outdir=lambda wc, output: os.path.dirname(str(output.bam)),
        pe_flag=lambda wc: "--pe" if layout_of(wc.sample) == "PE" else "",
    log: R("logs/deduplicate/{sample}.log"),
    threads: rthreads("deduplicate")
    resources:
        mem_mb=rmem("deduplicate"), runtime_min=rruntime("deduplicate"),
        runtime_sec=rruntime_sec("deduplicate"),
    shell:
        """
        deduplicate_bismark --bam {params.pe_flag} --od {params.outdir} \
            {input.bam} > {log} 2>&1
        """


rule bam2nuc_genome:
    input:
        genome=R("0.index/bismark_genome/Bisulfite_Genome/GA_conversion/BS_conv.1.bt2.lf"),
    output:
        R("0.index/bismark_genome/genomic_nucleotide_totals.txt"),
    params:
        genome_dir=lambda wc, input: os.path.dirname(os.path.dirname(
            os.path.dirname(str(input.genome)))),
    log: R("logs/bam2nuc_genome.log"),
    threads: rthreads("bam2nuc_genome")
    resources:
        mem_mb=rmem("bam2nuc_genome"), runtime_min=rruntime("bam2nuc_genome"),
        runtime_sec=rruntime_sec("bam2nuc_genome"),
    shell:
        """
        bam2nuc --genome_folder {params.genome_dir} \
            --genomic_composition_only > {log} 2>&1
        """


rule bam2nuc_sample:
    input:
        bam=R("4.dedup/{sample}.bam.deduplicated.bam"),
        totals=R("0.index/bismark_genome/genomic_nucleotide_totals.txt"),
    output:
        R("5.methylation/{sample}/{sample}.nucleotide_stats.txt"),
    params:
        genome_dir=lambda wc, input: os.path.dirname(os.path.dirname(
            os.path.dirname(str(input.totals)))),
        outdir=lambda wc, output: os.path.dirname(str(output)),
    log: R("logs/bam2nuc/{sample}.log"),
    threads: rthreads("bam2nuc_sample")
    resources:
        mem_mb=rmem("bam2nuc_sample"), runtime_min=rruntime("bam2nuc_sample"),
        runtime_sec=rruntime_sec("bam2nuc_sample"),
    shell:
        """
        bam2nuc --genome_folder {params.genome_dir} \
            --genomic_composition {input.totals} --dir {params.outdir} \
            {input.bam} > {log} 2>&1
        """
