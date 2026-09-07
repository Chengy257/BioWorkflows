# Optional novel-miRNA discovery stage (docs/TODO.md item 2): miRDeep-P2
# (bioconda mirdeep-p2=1.1.4) run per sample over the trimmed reads.
# Candidate hairpins are excised from the reference genome and scored by
# mapping signature; the tool-native result tree lands under
# results/6.novel_mirna/{sample}/ (main table {sample}_filter_P_prediction,
# plus the .bed export, precursor fasta/structures, signatures, and the
# {sample}_predictions raw scores -- file names are the tool's own and
# version-dependent). Parse-time guarded by _NOVEL_RUN from rules/common.smk:
# with the default novel_mirna.enabled=false this module declares nothing
# and the DAG is unchanged. Depends on rules/common.smk for R(), _NOVEL_RUN,
# and the rthreads/rmem/rruntime helpers.
#
# VERIFIED COMMAND (recorded per the implementation contract): the bioconda
# package mirdeep-p2=1.1.4 ships miRDP2-v1.1.4_pipeline.bash; the flag
# surface below was verified against the package sources (the exact
# miRDP2-v1.1.4.tar.gz artifact, sha256-matched to the bioconda recipe pin
# in bioconda-recipes/recipes/mirdeep-p2) and cross-checked with the
# miRDeep-P2 v1.1.4 usage documented at
# https://github.com/TF-Chan-Lab/miRDeep-P2_pipeline:
#   miRDP2-v1.1.4_pipeline.bash -f -g GENOME.fa -x INDEX_PREFIX -i READS.fa \
#       -o OUT_DIR -p THREADS
# Tool facts that shape this rule (all verified in the 1.1.4 sources):
# - -x requires a PRE-BUILT bowtie index prefix; the tool never builds a
#   genome index itself (it builds only its precursor index inside its own
#   output tree), so this rule consumes the workflow genome index
#   0.index/genome (bowtie1 small index, exactly the tool's default
#   expectation; only --large-index/.ebwtl users deviate).
# - There is NO mature-reference flag: known-miRNA filtering runs against
#   the plant mature-miRNA index bundled in the package, so the configured
#   novel_mirna.mature_fasta is declared as an input (existence +
#   provenance) but is not passed to the tool.
# - The -q fastq mode is broken in 1.1.4 (the script reassigns the
#   converted-fasta path and drops the per-sample directory component), and
#   scripts/preprocess_reads.pl requires collapsed
#   ">readNNNNNNNN_x<count>" headers for its RPM arithmetic, so this rule
#   converts the trimmed fq.gz into exactly that collapsed-fasta form and
#   feeds -f.
# - The rfam ncRNA index the tool references (scripts/index/rfam_index) is
#   not shipped in the package: that filter silently no-ops (its stderr
#   lands in 6.novel_mirna/{sample}/script_err); known-miRNA filtering
#   works. rRNA/tRNA reads therefore reach the scoring step.

if _NOVEL_RUN:
    rule novel_mirna:
        input:
            fq=R("2.cleandata/{sample}_trimmed.fq.gz"),
            genome=str(config["genome"]["fasta"]),
            genome_index=R("0.index/genome/genome.1.ebwt"),
            mature=str(_NOVEL_MATURE),
        output:
            fa=R("6.novel_mirna/{sample}.fa"),
            flag=R("6.novel_mirna/{sample}/flag.log"),
        log:
            R("logs/novel_mirna/{sample}.log"),
        params:
            mirdp2="miRDP2-v1.1.4_pipeline.bash",
            outdir=R("6.novel_mirna"),
            index_prefix=lambda wc, input: str(input.genome_index)[:-len(".1.ebwt")],
            final=lambda wc: R(f"6.novel_mirna/{wc.sample}/{wc.sample}_filter_P_prediction"),
        threads:
            rthreads("novel_mirna")
        resources:
            mem_mb=rmem("novel_mirna"),
            runtime_min=rruntime("novel_mirna"),
            runtime_sec=rruntime_sec("novel_mirna"),
        shell:
            """
            ## Collapsed unique-read fasta in the exact form miRDP2 expects
            ## (">readNNNNNNNN_x<count>"): sequence lines only, sorted, with
            ## duplicate counts encoded in the header (mirrors the tool's own
            ## fastq-to-fasta converter).
            gzip -cd {input.fq} | awk 'NR%4==2' | LC_ALL=C sort | uniq -c \\
                | awk '{{printf ">read%08d_x%d\\n%s\\n", NR, $1, $2}}' > {output.fa}
            {params.mirdp2} -f -g {input.genome} -x {params.index_prefix} \\
                -i {output.fa} -o {params.outdir} -p {threads} >> {log} 2>&1
            ## The tool ignores internal step failures (no set -e), so require
            ## the final plant-criteria prediction table to exist (it is
            ## created even when nothing passes the filters) before declaring
            ## the stage done.
            test -e {params.final}
            echo `date` " : novel miRNA prediction done" > {output.flag}
            """
