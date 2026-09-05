# rules/upstream.smk — optional adapter/quality trimming with Trim Galore
# (+ per-mate FastQC, driven by trim_galore's own --fastqc integration).
#
# A PE and an SE trim definition are both declared unconditionally and are
# selected per sample by the detected library layout: the PE definition
# takes the raw_pe(wc) pair via unpack and writes flat
# 2.cleandata/{sample}_1_val_1.fq.gz / {sample}_2_val_2.fq.gz plus
# per-mate trimming reports and per-mate fastqc zips; the SE definition
# takes raw_se(wc) and writes 2.cleandata/{sample}_trimmed.fq.gz plus
# report and fastqc zip. Output paths are FLAT under 2.cleandata/ (no
# trim/ subdirectory) — this is the locked align_input() contract in
# common.smk that the bismark alignment module consumes and that the QC
# aggregation in rules/meta.smk reports on. (rna-seq nests the same
# outputs under 2.cleandata/trim/ on purpose; do not align the two.)
#
# TRIM_ENABLED (common.smk) never gates these definitions: the plan note
# "when TRIM_ENABLED is false neither is defined" is realized
# DAG-equivalently by align_input(), which returns the raw fastqs when
# trim.enabled is false — the trimmed outputs are then never requested,
# so neither trim job is ever instantiated (raw_pe/raw_se also raise on a
# layout mismatch).
#
# The resources "fastqc" entry stays unused here by design: trim_galore's
# --fastqc runs inside the trim job envelope.

rule trim_pe:
    input:
        unpack(raw_pe),
    output:
        fq1=R("2.cleandata/{sample}_1_val_1.fq.gz"),
        fq2=R("2.cleandata/{sample}_2_val_2.fq.gz"),
        report1=R("2.cleandata/{sample}_1_val_1.fq.gz_trimming_report.txt"),
        report2=R("2.cleandata/{sample}_2_val_2.fq.gz_trimming_report.txt"),
        fastqc1=R("2.cleandata/fastqc/{sample}_1_val_1_fastqc.zip"),
        fastqc2=R("2.cleandata/fastqc/{sample}_2_val_2_fastqc.zip"),
    params:
        quality=config["trim"]["quality"],
        min_len=config["trim"]["min_len"],
        adapter=config["trim"]["adapter"],
        stringency=config["trim"]["stringency"],
        error_rate=config["trim"]["error_rate"],
        extra=config["trim"]["extra"],
        outdir=lambda wc, output: os.path.dirname(str(output.fq1)),
        fastqc_dir=lambda wc, output: os.path.dirname(str(output.fastqc1)),
    log:
        R("logs/trim_pe/{sample}.log"),
    threads: rthreads("trim")
    resources:
        mem_mb=rmem("trim"),
        runtime_min=rruntime("trim"),
        runtime_sec=rruntime_sec("trim"),
    shell:
        """
        mkdir -p {params.fastqc_dir}
        trim_galore --paired -q {params.quality} --stringency {params.stringency} \
            -e {params.error_rate} --length {params.min_len} -a {params.adapter} \
            -A {params.adapter} {params.extra} --gzip -j {threads} \
            -o {params.outdir}/ {input.fq1} {input.fq2} \
            --fastqc --fastqc_args "--outdir {params.fastqc_dir}" > {log} 2>&1
        ## Trim Galore (0.6.x) names the per-mate reports after the INPUT
        ## files; rename them to the sample-keyed contract paths above.
        mv -f "{params.outdir}/$(basename "{input.fq1}")_trimming_report.txt" {output.report1} 2>/dev/null || true
        mv -f "{params.outdir}/$(basename "{input.fq2}")_trimming_report.txt" {output.report2} 2>/dev/null || true
        """


rule trim_se:
    input:
        fq=raw_se,
    output:
        fq=R("2.cleandata/{sample}_trimmed.fq.gz"),
        report=R("2.cleandata/{sample}_trimming_report.txt"),
        fastqc=R("2.cleandata/fastqc/{sample}_trimmed_fastqc.zip"),
    params:
        quality=config["trim"]["quality"],
        min_len=config["trim"]["min_len"],
        adapter=config["trim"]["adapter"],
        stringency=config["trim"]["stringency"],
        error_rate=config["trim"]["error_rate"],
        extra=config["trim"]["extra"],
        outdir=lambda wc, output: os.path.dirname(str(output.fq)),
        fastqc_dir=lambda wc, output: os.path.dirname(str(output.fastqc)),
    log:
        R("logs/trim_se/{sample}.log"),
    threads: rthreads("trim")
    resources:
        mem_mb=rmem("trim"),
        runtime_min=rruntime("trim"),
        runtime_sec=rruntime_sec("trim"),
    shell:
        """
        mkdir -p {params.fastqc_dir}
        trim_galore -q {params.quality} --stringency {params.stringency} \
            -e {params.error_rate} --length {params.min_len} -a {params.adapter} \
            {params.extra} --gzip -j {threads} -o {params.outdir}/ {input.fq} \
            --fastqc --fastqc_args "--outdir {params.fastqc_dir}" > {log} 2>&1
        ## Trim Galore (0.6.x) names the report after the INPUT file;
        ## rename it to the sample-keyed contract path above (a no-op under
        ## tool versions that already write the canonical name).
        mv -f "{params.outdir}/$(basename "{input.fq}")_trimming_report.txt" {output.report} 2>/dev/null || true
        """
