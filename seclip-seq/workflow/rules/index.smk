# STAR genome / repeats index generation.
# Depends on workflow/rules/common.smk for R()/resource helpers and config.

rule star_index_genome:
    input:
        fasta=config["genome"],
        gtf=config["gtf"],
    output:
        R("0.index/genome_STARindex/SA"),
    params:
        outdir=lambda wc, output: os.path.dirname(str(output)),
        sjdb_overhang=config["star"]["sjdb_overhang"],
        sa_nbases=config["star"]["genome_sa_index_nbases"],
    log:
        R("logs/star_index_genome.log"),
    threads: rthreads("star_index_genome")
    resources:
        mem_mb=rmem("star_index_genome"),
        runtime_min=rruntime("star_index_genome"),
        runtime_sec=rruntime_sec("star_index_genome"),
    shell:
        """
        STAR --runThreadN {threads} --runMode genomeGenerate \
            --genomeDir {params.outdir} --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} --sjdbOverhang {params.sjdb_overhang} \
            --genomeSAindexNbases {params.sa_nbases} > {log} 2>&1
        """


rule star_index_repeats:
    input:
        fasta=config["repeats_fa"],
    output:
        R("0.index/repeats_STARindex/SA"),
    params:
        outdir=lambda wc, output: os.path.dirname(str(output)),
        sa_nbases=config["star"]["repeats_sa_index_nbases"],
        limit_ram=config["star"]["repeats_limit_ram"],
    log:
        R("logs/star_index_repeats.log"),
    threads: rthreads("star_index_repeats")
    resources:
        mem_mb=rmem("star_index_repeats"),
        runtime_min=rruntime("star_index_repeats"),
        runtime_sec=rruntime_sec("star_index_repeats"),
    shell:
        """
        STAR --runThreadN {threads} --runMode genomeGenerate \
            --genomeDir {params.outdir} --genomeFastaFiles {input.fasta} \
            --limitGenomeGenerateRAM {params.limit_ram} \
            --genomeSAindexNbases {params.sa_nbases} > {log} 2>&1
        """
