# Upstream shared steps: trim_galore SE trimming (with --fastqc) and
# bowtie1 index building for the running cascade classes and the genome.
# Depends on workflow/rules/common.smk for helpers.

rule trim:
    input:
        fq=lambda wc: raw_fastq(wc.sample),
    output:
        fq=R("2.cleandata/{sample}_trimmed.fq.gz"),
        report=R("2.cleandata/{sample}_trimming_report.txt"),
        fastqc_zip=R("2.cleandata/fastqc/{sample}_trimmed_fastqc.zip"),
    params:
        quality=config["trim"]["quality"],
        min_len=config["trim"]["min_len"],
        adapter=config["trim"]["adapter"],
        stringency=config["trim"]["stringency"],
        error_rate=config["trim"]["error_rate"],
        extra=config["trim"]["extra"],
        outdir=lambda wc, output: os.path.dirname(str(output.fq)),
        fastqc_dir=lambda wc, output: os.path.dirname(str(output.fastqc_zip)),
    log:
        R("logs/trim/{sample}.log"),
    threads: rthreads("trim")
    resources:
        mem_mb=rmem("trim"),
        runtime_min=rruntime("trim"),
        runtime_sec=rruntime_sec("trim"),
    shell:
        """
        trim_galore -q {params.quality} --stringency {params.stringency} \
            -e {params.error_rate} --length {params.min_len} -a {params.adapter} \
            {params.extra} --gzip -j {threads} -o {params.outdir}/ {input.fq} \
            --fastqc --fastqc_args "--outdir {params.fastqc_dir}" > {log} 2>&1
        ## Trim Galore (0.6.x) names the report after the INPUT file full basename; rename it to the contract path above.
        mv -f "{params.outdir}/$(basename "{input.fq}")_trimming_report.txt" {output.report} 2>/dev/null || true
        """


rule bowtie_index:
    input:
        lambda wc: config["genome"]["fasta"] if wc.klass == "genome" else cascade_fasta(wc.klass),
    output:
        R("0.index/{klass}/{klass}.1.ebwt"),
    params:
        prefix=lambda wc, output: str(output)[:-len(".1.ebwt")],
    log:
        R("logs/bowtie_index_{klass}.log"),
    threads: rthreads("bowtie_index")
    resources:
        mem_mb=rmem("bowtie_index"),
        runtime_min=rruntime("bowtie_index"),
        runtime_sec=rruntime_sec("bowtie_index"),
    wildcard_constraints:
        klass=_INDEX_WILDCARD,
    shell:
        """
        bowtie-build --threads {threads} {input} {params.prefix} > {log} 2>&1
        """
