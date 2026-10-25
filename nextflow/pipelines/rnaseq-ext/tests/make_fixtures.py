#!/usr/bin/env python3
"""Create the synthetic fixture tree for the rnaseq-ext `test` profiles.

The fixture layout mirrors the PINED nf-core/rnaseq 3.27.0 publish layout
(conf/modules/*.config): markdup BAMs under {aligner}/, the merged Salmon
matrix under salmon/. File contents are minimal but parseable so the
extension helper scripts can also run for real outside containers.

Fixture scientific contract (fixed seed):
- 6 samples: WT_1/WT_2 (control) + drugA_1/drugA_2 + drugB_1/drugB_2
- 40 genes; genes 10-19 up in drugA, 20-29 down (planted DE signal so the
  DA adapter / enrichment split helpers run on real numbers); genes 30-34 up
  in drugB only - the SECOND DE table pins the per-contrast enrichment
  pairing (>= 2 contrasts must not cross-stage tables)
- the DA fixture table carries DA's deseq2.results.tsv column layout
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures")


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


SAMPLES = ["WT_1", "WT_2", "drugA_1", "drugA_2", "drugB_1", "drugB_2"]
GENES = [f"gene{i:03d}" for i in range(1, 41)]

write(os.path.join(FIX, "samplesheet.csv"),
      "sample,fastq_1,fastq_2,strandedness\n" +
      "".join(f"{s},1.rawdata/{s}_1.fastq.gz,1.rawdata/{s}_2.fastq.gz,auto\n" for s in SAMPLES))


def group_of(s):
    if s.startswith("WT"):
        return "control"
    return s.rsplit("_", 1)[0]          # drugA / drugB


write(os.path.join(FIX, "metadata.csv"),
      "sample,group,batch\n" +
      "".join(f"{s},{group_of(s)},b1\n" for s in SAMPLES))

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
        if 29 <= i <= 33 and s.startswith("drugB"):
            base *= 4          # up in drugB only (second-contrast signal)
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

# DA differential tables (flat under tables/differential, matching the
# pinned differentialabundance 2.0.0 publish: prefix {contrast id} + the
# module's *.deseq2.results.tsv naming). TWO contrasts so the enrichment
# wiring is pinned at >= 2 tables (the cross-staging failure mode).
rng.seed(7)


def de_table(contrast, up_rows, down_rows):
    rows = ["gene_id\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj"]
    for i, g in enumerate(GENES):
        if i in up_rows:
            lfc, padj = 2.5, 0.001          # significant UP
        elif i in down_rows:
            lfc, padj = -2.0, 0.003         # significant DOWN
        else:
            lfc, padj = round(rng.uniform(-0.5, 0.5), 3), round(rng.uniform(0.1, 1.0), 3)
        stat = lfc / 0.5
        rows.append(f"{g}\t{lfc}\t0.5\t{stat}\t{padj / 10:.6f}\t{padj}")
    write(os.path.join(FIX, f"da_out/tables/differential/{contrast}.deseq2.results.tsv"),
          "\n".join(rows) + "\n")


de_table("drugA_vs_control", set(range(9, 19)), set(range(19, 29)))
de_table("drugB_vs_control", set(range(29, 34)), set())

print(f"fixtures written to {FIX}")
