#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate the miniature RNA-seq regression test dataset.

Pure standard library implementation with a fixed random seed; data is fully reproducible
for a given version (the data is not committed, it is generated on the fly on the test
machine). Outputs (--outdir, default tests/data/):

    reference/genome.fa           2 x 100 kb chromosomes
    reference/genes.gtf           60 simulated genes (gene/transcript/exon)
    reference/genes.bed           BED12 gene models (used by RSeQC infer_experiment)
    reference/annotation_full.tsv gene_id functional annotation (used for DEG annotation)
    rawdata/{sample}_1.fastq.gz   simulated PE reads (default 50000 pairs per sample;
    rawdata/{sample}_2.fastq.gz   sampled from transcripts + sequencing errors + 10% adapter + 5% noise)
    truth_degenes.tsv             expected differential gene directions (drugA: 6 up x8 / 6 down x0.125)
    samples.csv                   sample table (2 groups x 2 replicates)

Samples: ctrl_1/ctrl_2 (control group) and drugA_1/drugA_2 (treatment group).

Usage:
    python3 make_testdata.py [--outdir tests/data] [--reads 50000] [--seed 20260903]
"""
import argparse
import gzip
import os
import random

CHROM_LEN = 100000
N_CHROM = 2
N_GENES = 60
READ_LEN = 150
FRAG_MIN, FRAG_MAX = 280, 420
ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ERROR_RATE = 0.005
ADAPTER_FRAC = 0.10
NOISE_FRAC = 0.05
N_UP, N_DOWN = 6, 6
SAMPLES = [("ctrl_1", "control", 1.00),
           ("ctrl_2", "control", 1.05),
           ("drugA_1", "drugA", 1.00),
           ("drugA_2", "drugA", 1.10)]
FOLD_UP, FOLD_DOWN = 8.0, 0.125

BASES = "ACGT"
_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def mutate(seq, rng):
    """Introduce substitution-type sequencing errors at ERROR_RATE."""
    if ERROR_RATE <= 0:
        return seq
    out = []
    for b in seq:
        if rng.random() < ERROR_RATE:
            out.append(rng.choice(BASES))
        else:
            out.append(b)
    return "".join(out)


def build_reference():
    """Return (chromosome sequence dict, gene structure list)."""
    rng = random.Random(777)
    chroms = {f"chr{i + 1}": "".join(rng.choices(BASES, k=CHROM_LEN))
              for i in range(N_CHROM)}
    genes = []
    per_chrom = N_GENES // N_CHROM
    for i in range(N_GENES):
        chrom = f"chr{i // per_chrom + 1}"
        idx = i % per_chrom
        strand = "+" if i % 2 == 0 else "-"
        gid = f"gene{i + 1:04d}"
        tx_start = 500 + idx * 3000          # 0-based
        exons = []                            # (0-based start, end)
        for j in range(5):
            s = tx_start + j * 550
            exons.append((s, s + 300))
        genes.append({"id": gid, "chrom": chrom, "strand": strand,
                      "tx_start": exons[0][0], "tx_end": exons[-1][1],
                      "exons": exons})
    return chroms, genes


def transcript_seq(chroms, gene):
    """Concatenate exon sequences in transcript orientation."""
    seq = chroms[gene["chrom"]]
    exons = gene["exons"]
    if gene["strand"] == "+":
        return "".join(seq[s:e] for s, e in exons)
    return "".join(revcomp(seq[s:e]) for s, e in reversed(exons))


def write_reference(outdir, chroms, genes):
    ref_dir = os.path.join(outdir, "reference")
    os.makedirs(ref_dir, exist_ok=True)

    with open(os.path.join(ref_dir, "genome.fa"), "w") as fh:
        for c in sorted(chroms):
            fh.write(f">{c}\n")
            s = chroms[c]
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")

    with open(os.path.join(ref_dir, "genes.gtf"), "w") as fh:
        for g in genes:
            attr_g = f'gene_id "{g["id"]}";'
            fh.write(f'{g["chrom"]}\tsim\tgene\t{g["tx_start"] + 1}\t{g["tx_end"]}\t.\t{g["strand"]}\t.\t{attr_g};\n')
            fh.write(f'{g["chrom"]}\tsim\ttranscript\t{g["tx_start"] + 1}\t{g["tx_end"]}\t.\t{g["strand"]}\t.\t{attr_g} transcript_id "{g["id"]}.1";\n')
            for s, e in g["exons"]:
                fh.write(f'{g["chrom"]}\tsim\texon\t{s + 1}\t{e}\t.\t{g["strand"]}\t.\t{attr_g} transcript_id "{g["id"]}.1";\n')

    with open(os.path.join(ref_dir, "genes.bed"), "w") as fh:
        for g in genes:
            blocks = sorted(g["exons"])           # BED blocks in ascending genomic coordinate order
            sizes = [e - s for s, e in blocks]
            starts = [s - g["tx_start"] for s, e in blocks]
            fh.write("\t".join(str(x) for x in [
                g["chrom"], g["tx_start"], g["tx_end"], g["id"], "0", g["strand"],
                g["tx_start"], g["tx_end"], "0", len(blocks),
                ",".join(map(str, sizes)) + ",", ",".join(map(str, starts)) + ","]) + "\n")

    with open(os.path.join(ref_dir, "annotation_full.tsv"), "w") as fh:
        for i, g in enumerate(genes):
            fh.write(f'{g["id"]}\tsimulated gene {g["id"]}, function class {(i % 5) + 1}\n')


def write_fastqs(outdir, chroms, genes, reads_per_sample, seed):
    """Sample simulated PE reads according to expression levels."""
    raw_dir = os.path.join(outdir, "rawdata")
    os.makedirs(raw_dir, exist_ok=True)
    txs = [transcript_seq(chroms, g) for g in genes]
    tx_len = [len(t) for t in txs]

    up_idx = set(range(10, 10 + N_UP))
    down_idx = set(range(30, 30 + N_DOWN))

    truth = os.path.join(outdir, "truth_degenes.tsv")
    with open(truth, "w") as fh:
        fh.write("gene_id\texpected_drugA_vs_control\n")
        for i in range(len(genes)):
            if i in up_idx:
                fh.write(f"{genes[i]['id']}\tup\n")
            elif i in down_idx:
                fh.write(f"{genes[i]['id']}\tdown\n")

    for k, (sample, group, factor) in enumerate(SAMPLES):
        rng = random.Random(seed + k)
        weights = []
        for i in range(len(genes)):
            w = rng.uniform(0.6, 1.5)
            if group == "drugA":
                if i in up_idx:
                    w *= FOLD_UP
                elif i in down_idx:
                    w *= FOLD_DOWN
            weights.append(w * factor)

        fq1 = gzip.open(os.path.join(raw_dir, f"{sample}_1.fastq.gz"), "wt", compresslevel=6)
        fq2 = gzip.open(os.path.join(raw_dir, f"{sample}_2.fastq.gz"), "wt", compresslevel=6)
        try:
            for n in range(reads_per_sample):
                q = "I" * READ_LEN
                if rng.random() < NOISE_FRAC:
                    r1 = "".join(rng.choices(BASES, k=READ_LEN))
                    r2 = revcomp("".join(rng.choices(BASES, k=READ_LEN)))
                else:
                    gi = rng.choices(range(len(genes)), weights=weights)[0]
                    tx = txs[gi]
                    frag_len = rng.randint(FRAG_MIN, min(FRAG_MAX, tx_len[gi]))
                    start = rng.randint(0, tx_len[gi] - frag_len)
                    frag = tx[start:start + frag_len]
                    r1 = mutate(frag[:READ_LEN], rng)
                    r2 = revcomp(mutate(frag[-READ_LEN:], rng))
                    if rng.random() < ADAPTER_FRAC:
                        r1 += ADAPTER[:rng.randint(15, 25)]
                    if rng.random() < ADAPTER_FRAC:
                        r2 += ADAPTER[:rng.randint(15, 25)]
                fq1.write(f"@r{n:07d} 1:N:0:1\n{r1}\n+\n{q}\n")
                fq2.write(f"@r{n:07d} 2:N:0:1\n{r2}\n+\n{q}\n")
        finally:
            fq1.close()
            fq2.close()

    return truth


def write_samples(outdir):
    with open(os.path.join(outdir, "samples.csv"), "w") as fh:
        fh.write("id,group,layout\n")
        for sample, group, _ in SAMPLES:
            fh.write(f"{sample},{group},auto\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=os.path.join(os.path.dirname(__file__), "data"))
    ap.add_argument("--reads", type=int, default=50000, help="PE read pairs per sample")
    ap.add_argument("--seed", type=int, default=20260903)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    chroms, genes = build_reference()
    write_reference(args.outdir, chroms, genes)
    truth = write_fastqs(args.outdir, chroms, genes, args.reads, args.seed)
    write_samples(args.outdir)

    with open(truth) as fh:
        n_truth = sum(1 for _ in fh) - 1
    print(f"[make_testdata] genes {len(genes)}, chromosomes {N_CHROM} x {CHROM_LEN}bp, "
          f"samples {len(SAMPLES)} x {args.reads} PE read pairs")
    print(f"[make_testdata] truth differential genes {n_truth} -> {truth}")
    print(f"[make_testdata] output root: {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()
