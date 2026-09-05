# rules/methylation.smk — methylation extraction, CpG merge, per-sample
# reports, and the run-level Bismark summary.
#
# Downstream of rules/align.smk: the extractor runs on the deduplicated BAM
# ({sample}.deduplicated.bam) and inherits its basename, so every extractor
# output carries the ".deduplicated" infix (the derived-name convention
# locked in TARGETS, common.smk, and the align.smk header):
#   5.methylation/{sample}/{sample}.deduplicated.bismark.cov.gz
#   5.methylation/{sample}/{sample}.deduplicated.bedGraph.gz
#   5.methylation/{sample}/{sample}.deduplicated_splitting_report.txt
#   5.methylation/{sample}/{sample}.deduplicated.M-bias.txt
# The extractor may also drop a side-effect splitter report next to them
# (not a declared output). coverage2cytosine --merge_CpG writes
# 5.methylation/{sample}/{sample}.CpG_merged.tsv.gz; the
# definition is unconditional but only requested when
# methylation_extractor.merge_cpg is true (via TARGETS in common.smk).
#
# Phase D real-run validation against the pinned bismark 0.24.0 (tracked in
# docs/TODO.md) must confirm: the exact extractor output names, the exact
# coverage2cytosine -o / ".gz" naming (the trailing `gzip -f` line is the
# safe fallback and stays until confirmed), the bismark2report html name
# ({sample}.html: bismark2report -o X writes X.html),
# that bismark2summary (-o bismark2summary) writes into the working
# directory (hence the `cd` into 5.QC/ below; a --basename alternative is
# documented in docs/TODO.md).

rule methylation_extractor:
    input:
        bam=R("4.dedup/{sample}.deduplicated.bam"),
        genome=R("0.index/bismark_genome/Bisulfite_Genome/GA_conversion/BS_GA.1.bt2"),
    output:
        # Derived names: the extractor inherits the full dedup basename.
        cov=R("5.methylation/{sample}/{sample}.deduplicated.bismark.cov.gz"),
        bedgraph=R("5.methylation/{sample}/{sample}.deduplicated.bedGraph.gz"),
        splitting=R("5.methylation/{sample}/{sample}.deduplicated_splitting_report.txt"),
        mbias=R("5.methylation/{sample}/{sample}.deduplicated.M-bias.txt"),
    params:
        genome_dir=lambda wc, input: os.path.dirname(os.path.dirname(
            os.path.dirname(str(input.genome)))),
        outdir=lambda wc, output: os.path.dirname(str(output.cov)),
        pe_flag=lambda wc: "-p" if layout_of(wc.sample) == "PE" else "",
        buffer_size=lambda wc: f"{max(1, rmem('methylation_extractor') // 1024 // config['methylation_extractor']['buffer_frac'])}G",
        cx_flag=lambda wc: "--CX --cytosine_report" if config["methylation_extractor"]["cx_report"] else "",
    log: R("logs/methylation_extractor/{sample}.log"),
    threads: rthreads("methylation_extractor")
    resources:
        mem_mb=rmem("methylation_extractor"), runtime_min=rruntime("methylation_extractor"),
        runtime_sec=rruntime_sec("methylation_extractor"),
    shell:
        """
        bismark_methylation_extractor {params.pe_flag} \
            --genome {params.genome_dir} --gzip --bedGraph \
            --buffer_size {params.buffer_size} {params.cx_flag} \
            --output_dir {params.outdir} {input.bam} > {log} 2>&1
        """


rule coverage2cytosine:
    input:
        cov=R("5.methylation/{sample}/{sample}.deduplicated.bismark.cov.gz"),
        genome=R("0.index/bismark_genome/Bisulfite_Genome/GA_conversion/BS_GA.1.bt2"),
    output:
        R("5.methylation/{sample}/{sample}.CpG_merged.CpG_report.merged_CpG_evidence.cov.gz"),
    params:
        genome_dir=lambda wc, input: os.path.dirname(os.path.dirname(
            os.path.dirname(str(input.genome)))),
        outdir=lambda wc, output: os.path.dirname(str(output)),
        prefix=lambda wc: f"{wc.sample}.CpG_merged",
    log: R("logs/coverage2cytosine/{sample}.log"),
    threads: rthreads("coverage2cytosine")
    resources:
        mem_mb=rmem("coverage2cytosine"), runtime_min=rruntime("coverage2cytosine"),
        runtime_sec=rruntime_sec("coverage2cytosine"),
    shell:
        """
        coverage2cytosine --genome_folder {params.genome_dir} --merge_CpG \
            -o {params.prefix} --dir {params.outdir} {input.cov} > {log} 2>&1
        ## coverage2cytosine writes {{prefix}}.CpG_report.merged_CpG_evidence.cov
        ## (--merge_CpG); gzip it to the declared .cov.gz output. The main
        ## {{prefix}}.CpG_report.txt cytosine report stays beside it (side effect).
        gzip -f {params.outdir}/{params.prefix}.CpG_report.merged_CpG_evidence.cov
        """


rule bismark2report:
    input:
        align_report=R("3.align/{sample}_report.txt"),
        dedup_report=R("4.dedup/{sample}.deduplication_report.txt"),
        splitting=R("5.methylation/{sample}/{sample}.deduplicated_splitting_report.txt"),
        mbias=R("5.methylation/{sample}/{sample}.deduplicated.M-bias.txt"),
        nuc=R("5.methylation/{sample}/{sample}.deduplicated.nucleotide_stats.txt"),
    output:
        R("5.methylation/{sample}/{sample}.html"),
    params:
        outdir=lambda wc, output: os.path.dirname(str(output)),
    log: R("logs/bismark2report/{sample}.log"),
    threads: rthreads("bismark2report")
    resources:
        mem_mb=rmem("bismark2report"), runtime_min=rruntime("bismark2report"),
        runtime_sec=rruntime_sec("bismark2report"),
    shell:
        """
        bismark2report --dir {params.outdir} --output {wildcards.sample} \
            --alignment_report {input.align_report} \
            --dedup_report {input.dedup_report} \
            --splitting_report {input.splitting} \
            --mbias_report {input.mbias} \
            --nucleotide_report {input.nuc} > {log} 2>&1
        """


rule bismark2summary:
    input:
        # Declared inputs are DAG dependency edges only: bismark2summary
        # discovers the per-sample reports itself when run inside 5.QC/;
        # the dedup entry exists purely to order the job after dedup.
        reports=expand(R("3.align/{sample}_report.txt"), sample=SAMPLES),
        dedup=expand(R("4.dedup/{sample}.deduplication_report.txt"), sample=SAMPLES),
    output:
        R("5.QC/bismark2summary.html"),
    params:
        outdir=lambda wc, output: os.path.dirname(str(output)),
        # bismark2summary takes ALIGNMENT BAM files whose basenames end in
        # _pe/_se (absolute paths; after the `cd` below the relative paths
        # would no longer resolve) and reads each native report
        # <base>_PE_report.txt / _SE_report.txt beside the BAM — the symlink
        # + native-name copy kept by rules/align.smk. Dedup/splitting stats
        # are skipped unless reports sit beside the BAM under the tool's
        # own names (v0.1 limitation, docs/TODO.md).
        abs_bams=" ".join(
            os.path.abspath(R(f"3.align/{s}_{'pe' if layout_of(s) == 'PE' else 'se'}.bam"))
            for s in SAMPLES),
        abs_log=lambda wc: os.path.abspath(R("logs/bismark2summary.log")),
    log: R("logs/bismark2summary.log"),
    threads: rthreads("bismark2summary")
    resources:
        mem_mb=rmem("bismark2summary"), runtime_min=rruntime("bismark2summary"),
        runtime_sec=rruntime_sec("bismark2summary"),
    shell:
        """
        cd {params.outdir} && bismark2summary -o bismark2summary {params.abs_bams} > {params.abs_log} 2>&1
        """
