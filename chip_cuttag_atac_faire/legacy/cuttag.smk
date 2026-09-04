configfile: "config.yaml"

def get_samples():
    ids = []
    gps = []
    with open(config["grouplist"], "r") as samples_list:
        next(samples_list)
        for line in samples_list:
            line = line.strip().split(",")
            # ids.append(line[0])
            if line[1] not in ids:
                ids.append(line[0])
            if line[3] not in gps:
                gps.append(line[3])
    samples_list.close()
    return ids,gps
SAMPLES,GROUPS = get_samples()

rule all:
    input:
        "2.cleandata/fastqc/multiqc/multiqc_report.html",
        expand("3.align/bowtie2/{sample}_sorted.bam",sample=SAMPLES),
        # expand("3.align/bowtie2/{sample}_rmdup.bam",sample=SAMPLES),
        # expand("4.peak/{sample}_peaks.narrowPeak",sample=SAMPLES),
        # "4.peak/anno_result/Peakanno_PeakDistributions.pdf",
        # "5.QC_deeptools/matrix_scaled.gz",
        # "6.motif/flag.txt",

include: "/home/chengyu/workflows/snakemake/chip_cuttag_atac_faire-workflows/rules/upsteam.smk"


