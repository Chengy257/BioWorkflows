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
# Derived-name convention these definitions emit (reconciled against the
# pinned bismark 0.24.0 in WSL real-run validation, 2026-09-05/06; keep in
# sync with TARGETS in common.smk and the methylation module):
# - bismark --basename {sample} writes {sample}_pe.bam/_se.bam plus
#   {sample}_PE_report.txt/_SE_report.txt into --output_dir; the rule
#   renames/copies them to the layout-neutral {sample}.bam +
#   {sample}_report.txt and keeps the native copies + a {sample}_pe.bam
#   symlink for bismark2summary discovery (see the rule body);
# - deduplicate_bismark writes {sample}.deduplicated.bam plus
#   {sample}.deduplication_report.txt into --output_dir (4.dedup/);
# - bam2nuc names the per-sample output after the input BAM:
#   {sample}.deduplicated.nucleotide_stats.txt under 5.methylation/{sample}/
#   (composition auto-detected from the genome folder — never pass
#   --genomic_composition, see the bam2nuc_sample rule body).
# No per-job --parallel: bismark 0.24 rejects --basename together with
# --multicore (docs/TODO.md §1); keep --phred33-quals in bismark.align_extra.

rule bismark_align:
    input:
        fq=lambda wc: align_input(wc.sample),
        genome=R("0.index/bismark_genome/Bisulfite_Genome/GA_conversion/BS_GA.1.bt2"),
    output:
        bam=R("3.align/{sample}.bam"),
        report=R("3.align/{sample}_report.txt"),
    params:
        genome_dir=lambda wc, input: os.path.dirname(os.path.dirname(
            os.path.dirname(str(input.genome)))),
        outdir=lambda wc, output: os.path.dirname(str(output.bam)),
        reads=lambda wc: align_input_arg(wc.sample),
        extra=config["bismark"]["align_extra"],
        layout_tag=lambda wc: "pe" if layout_of(wc.sample) == "PE" else "se",
    log: R("logs/bismark_align/{sample}.log"),
    threads: rthreads("bismark_align")
    resources:
        mem_mb=rmem("bismark_align"), runtime_min=rruntime("bismark_align"),
        runtime_sec=rruntime_sec("bismark_align"),
    shell:
        """
        ## No --parallel: bismark 0.24 rejects --basename together with
        ## --multicore; per-job multicore is deferred to v0.2 (docs/TODO.md)
        ## and parallelism comes from Snakemake scheduling samples
        ## concurrently (the legacy script parallelized the same way via
        ## ParaFly across samples).
        bismark --genome_folder {params.genome_dir} \
            {params.extra} --basename {wildcards.sample} \
            --bam --output_dir {params.outdir} {params.reads} > {log} 2>&1
        ## Bismark 0.24 with --basename X writes X_pe.bam / X_se.bam
        ## (lowercase layout tag) and X_PE_report.txt / X_SE_report.txt
        ## (uppercase). The BAM is renamed to the layout-neutral contract
        ## path; the report is COPIED (native name kept) because
        ## bismark2summary classifies alignment reports by the native
        ## _PE_report.txt / _SE_report.txt suffix.
        mv -f {params.outdir}/{wildcards.sample}_pe.bam {output.bam} 2>/dev/null || true
        mv -f {params.outdir}/{wildcards.sample}_se.bam {output.bam} 2>/dev/null || true
        cp -f {params.outdir}/{wildcards.sample}_PE_report.txt {output.report} 2>/dev/null || true
        cp -f {params.outdir}/{wildcards.sample}_SE_report.txt {output.report} 2>/dev/null || true
        ## bismark2summary requires BAM basenames ending in _pe/_se to
        ## locate the native report; keep a zero-copy symlink under the
        ## layout tag pointing at the contract BAM (same directory).
        ln -sfn $(basename "{output.bam}") {params.outdir}/{wildcards.sample}_{params.layout_tag}.bam
        """


rule deduplicate:
    input:
        bam=R("3.align/{sample}.bam"),
    output:
        # deduplicate_bismark derives the output name from the input file
        # name and writes it into --od: {sample}.bam -> {sample}.deduplicated.bam
        bam=R("4.dedup/{sample}.deduplicated.bam"),
        report=R("4.dedup/{sample}.deduplication_report.txt"),
    params:
        outdir=lambda wc, output: os.path.dirname(str(output.bam)),
        pe_flag=lambda wc: "-p" if layout_of(wc.sample) == "PE" else "",
    log: R("logs/deduplicate/{sample}.log"),
    threads: rthreads("deduplicate")
    resources:
        mem_mb=rmem("deduplicate"), runtime_min=rruntime("deduplicate"),
        runtime_sec=rruntime_sec("deduplicate"),
    shell:
        """
        deduplicate_bismark --bam {params.pe_flag} --output_dir {params.outdir} \
            {input.bam} > {log} 2>&1
        """


rule bam2nuc_genome:
    input:
        genome=R("0.index/bismark_genome/Bisulfite_Genome/GA_conversion/BS_GA.1.bt2"),
    output:
        R("0.index/bismark_genome/genomic_nucleotide_frequencies.txt"),
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
        bam=R("4.dedup/{sample}.deduplicated.bam"),
        totals=R("0.index/bismark_genome/genomic_nucleotide_frequencies.txt"),
    output:
        R("5.methylation/{sample}/{sample}.deduplicated.nucleotide_stats.txt"),
    params:
        # The totals file sits directly inside the genome folder (one
        # dirname up), unlike the deep Bisulfite_Genome sentinel.
        genome_dir=lambda wc, input: os.path.dirname(str(input.totals)),
        outdir=lambda wc, output: os.path.dirname(str(output)),
    log: R("logs/bam2nuc/{sample}.log"),
    threads: rthreads("bam2nuc_sample")
    resources:
        mem_mb=rmem("bam2nuc_sample"), runtime_min=rruntime("bam2nuc_sample"),
        runtime_sec=rruntime_sec("bam2nuc_sample"),
    shell:
        """
        ## Do NOT pass --genomic_composition: bam2nuc has no such option,
        ## and Perl Getopt abbreviation silently resolves it to
        ## --genomic_composition_only (genome composition, exit, BAM never
        ## processed). The composition file is auto-detected from the
        ## genome folder; the totals input stays as the DAG ordering edge.
        bam2nuc --genome_folder {params.genome_dir} \
            --dir {params.outdir} \
            {input.bam} > {log} 2>&1
        """
