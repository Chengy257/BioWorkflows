#!/usr/bin/env python3
"""Create the synthetic fixture tree for the rnaseq-ext `test` profiles.

The fixture layout mirrors the PINED nf-core/rnaseq 3.27.0 publish layout
(conf/modules/*.config): markdup BAMs under {aligner}/, the merged Salmon
matrix under salmon/. File contents are minimal but parseable so the
extension helper scripts can also run for real outside containers.

Fixture scientific contract (fixed seed):
- 4 samples: WT_1/WT_2 (control) + drugA_1/drugA_2 (treatment)
- 40 genes; genes 10-19 up in drugA, 20-29 down (planted DE signal so the
  DA adapter / enrichment split helpers run on real numbers)
- the DA fixture table carries DA's deseq2.results.tsv column layout
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures")


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


SAMPLES = ["WT_1", "WT_2", "drugA_1", "drugA_2"]
GENES = [f"gene{i:03d}" for i in range(1, 41)]

write(os.path.join(FIX, "samplesheet.csv"),
      "sample,fastq_1,fastq_2,strandedness\n" +
      "".join(f"{s},1.rawdata/{s}_1.fastq.gz,1.rawdata/{s}_2.fastq.gz,auto\n" for s in SAMPLES))

write(os.path.join(FIX, "metadata.csv"),
      "sample,group,batch\n" +
      "".join(f"{s},{'control' if s.startswith('WT') else 'drugA'},b1\n" for s in SAMPLES))

# Placeholder BAMs: file existence is what dry-runs validate (real BAMs for
# real-container smokes come from tests/make_real_bams.sh when needed).
for s in SAMPLES:
    for suffix in (".markdup.sorted.bam", ".markdup.sorted.bam.bai"):
        write(os.path.join(FIX, "results/star_salmon", s + suffix), "BAM\n")

# Native Salmon matrix (authoritative matrix on the route=salmon side)
rng = __import__("random")
rng.seed(42)
salmon_cols = "\t".join(SAMPLES)
lines = ["gene_id\tgene_name\t" + salmon_cols]
for i, g in enumerate(GENES):
    row = [g, g]
    for s in SAMPLES:
        base = 10
        if 9 <= i <= 18 and s.startswith("drugA"):
            base *= 4          # up in treatment
        if 19 <= i <= 28 and s.startswith("drugA"):
            base = max(1, base // 4)   # down in treatment
        row.append(str(base + rng.randint(0, 3)))
    lines.append("\t".join(row))
write(os.path.join(FIX, "results/salmon/salmon.merged.gene_counts.tsv"), "\n".join(lines) + "\n")
write(os.path.join(FIX, "results/salmon/salmon.merged.gene_tpm.tsv"), "\n".join(lines) + "\n")

# Reference: a tiny GTF (exon rows for all 40 genes) + matching BED12
gtf_rows = []
bed_rows = []
for i, g in enumerate(GENES):
    start = 1000 * (i + 1)
    end = start + 900
    gtf_rows.append(f"1\tfix\texon\t{start}\t{end}\t.\t+\t.\tgene_id \"{g}\"; transcript_id \"{g}.1\";")
    bed_rows.append(f"1\t{start - 1}\t{end}\t{g}\t0\t+")
write(os.path.join(FIX, "reference/genes.gtf"), "\n".join(gtf_rows) + "\n")
write(os.path.join(FIX, "reference/genes.bed"), "\n".join(bed_rows) + "\n")

# DA differential table for one contrast (flat under tables/differential,
# matching the pinned differentialabundance 2.0.0 publish: prefix {contrast id}
# + the module's *.deseq2.results.tsv naming)
rng.seed(7)
de_rows = ["gene_id\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj"]
for i, g in enumerate(GENES):
    if 9 <= i <= 18:
        lfc, padj = 2.5, 0.001          # significant UP
    elif 19 <= i <= 28:
        lfc, padj = -2.0, 0.003         # significant DOWN
    else:
        lfc, padj = round(rng.uniform(-0.5, 0.5), 3), round(rng.uniform(0.1, 1.0), 3)
    stat = lfc / 0.5
    de_rows.append(f"{g}\t{lfc}\t0.5\t{stat}\t{padj / 10:.6f}\t{padj}")
write(os.path.join(FIX, "da_out/tables/differential/drugA_vs_control.deseq2.results.tsv"),
      "\n".join(de_rows) + "\n")

print(f"fixtures written to {FIX}")
