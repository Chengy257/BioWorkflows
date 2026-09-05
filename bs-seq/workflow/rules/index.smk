# rules/index.smk — bisulfite genome preparation + genome nucleotide composition.
#
# bismark_genome_preparation requires the fasta INSIDE the genome folder,
# so the species fasta named by config["genome"] is copied into
# 0.index/bismark_genome/ first (a plain copy is cheap and cluster-safe)
# and the preparation runs on that folder. The folder may also hold
# multiple contig fastas: bismark_genome_preparation picks up every fasta
# present, so extra contig files can simply be dropped next to the copy.
# The GA-conversion bowtie2 index below doubles as the alignment-side
# sentinel consumed by the bismark alignment definition in rules/align.smk.
rule bismark_genome_prep:
    input:
        fasta=config["genome"],
    output:
        R("0.index/bismark_genome/Bisulfite_Genome/GA_conversion/BS_conv.1.bt2.lf"),
    params:
        genome_dir=lambda wc, output: os.path.dirname(os.path.dirname(os.path.dirname(str(output)))),
    log: R("logs/bismark_genome_prep.log"),
    threads: rthreads("bismark_genome_prep")
    resources:
        mem_mb=rmem("bismark_genome_prep"), runtime_min=rruntime("bismark_genome_prep"),
        runtime_sec=rruntime_sec("bismark_genome_prep"),
    shell:
        """
        mkdir -p {params.genome_dir}
        cp {input.fasta} {params.genome_dir}/
        bismark_genome_preparation --parallel {threads} --verbose \
            {params.genome_dir} > {log} 2>&1
        """
