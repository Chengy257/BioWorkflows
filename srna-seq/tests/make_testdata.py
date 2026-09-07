#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate the miniature srna-seq regression test dataset.

Pure standard library implementation with a fixed random seed and the gzip
mtime pinned to 0, so outputs are byte-reproducible for a given version (the
data is not committed, it is generated on the fly on the test machine). A
single master RNG stream (seeded 42) builds the references first and then
draws every read for every sample in order. Outputs (--outdir, default
tests/data/):

    ref/genome.fa             2 x 20 kb chromosomes (chr1/chr2, random ACGT)
    ref/rRNA.fa               2 sequences x 200 bp
    ref/tRNA.fa               3 sequences x 75 bp
    ref/miRNA.fa              4 sequences x 21 nt (osa-MIRmock1..4)
    1.rawdata/{sample}.fq.gz  simulated SE reads of 21-26 nt: 55% miRNA-derived
                              (full 21 nt sequence, 0-1 substitution
                              mismatches), 15% rRNA-derived, 10% tRNA-derived,
                              20% intergenic genome substrings; 35% of the
                              reads carry a shifted 3' adapter tail
                              (AGATCGGAAGAGC + shifts) (default 2000 reads per
                              sample, --reads N)
    samples.csv               sample table (sample_id column)
    config.yaml               ready-to-run miniature configuration (relative
                              ref/ paths; cascade rRNA/tRNA/miRNA filled,
                              snoRNA/snRNA/mRNA/rhizo empty -- the dry-run
                              regression therefore also exercises the
                              empty-fasta skip path)

Samples: s1, s2.

Usage:
    python3 make_testdata.py [--outdir tests/data] [--reads 2000] [--seed 42]
"""
import argparse
import gzip
import io
import os
import random

CHROM_LEN = 20000
N_CHROM = 2
N_RRNA = 2
RRNA_LEN = 200
N_TRNA = 3
TRNA_LEN = 75
N_MIRNA = 4
MIRNA_LEN = 21
MIN_READ = 21            # simulated mature sRNA read length range
MAX_READ = 26
ADAPTER = "AGATCGGAAGAGC"   # trim.adapter from config.yaml
MIRNA_FRAC = 0.55        # miRNA-derived reads (0-1 substitution mismatches)
RRNA_FRAC = 0.15         # rRNA-derived reads
TRNA_FRAC = 0.10         # tRNA-derived reads (intergenic genome takes the rest)
ADAPTER_FRAC = 0.35      # fraction of reads carrying a 3' adapter tail
SEED = 42
SAMPLES = ["s1", "s2"]

BASES = "ACGT"

# Miniature configuration written next to the generated data. Relative ref/
# paths; the cascade lists every class in workflow order with rRNA/tRNA/miRNA
# filled and snoRNA/snRNA/mRNA/rhizo left empty (the empty-fasta skip path is
# exercised by the dry-run regression).
CONFIG_TEMPLATE = """\
SampleListFile: "samples.csv"
results_dir: "results"
threads: 2
species: "none"
trim:
  quality: 25
  min_len: 15
  adapter: "AGATCGGAAGAGC"
  stringency: 3
  error_rate: 0.1
  extra: ""
cascade:
  - name: rRNA
    fasta: "ref/rRNA.fa"
  - name: snoRNA
    fasta: ""
  - name: snRNA
    fasta: ""
  - name: tRNA
    fasta: "ref/tRNA.fa"
  - name: miRNA
    fasta: "ref/miRNA.fa"
  - name: mRNA
    fasta: ""
  - name: rhizo
    fasta: ""
genome:
  fasta: "ref/genome.fa"
bowtie:
  extra: ""
"""


def build_reference(rng):
    """Draw chromosomes and the rRNA/tRNA/miRNA class fastas from the master
    RNG stream."""
    chroms = {f"chr{i + 1}": "".join(rng.choices(BASES, k=CHROM_LEN))
              for i in range(N_CHROM)}
    rrnas = [("rRNA_mock%d" % (i + 1), "".join(rng.choices(BASES, k=RRNA_LEN)))
             for i in range(N_RRNA)]
    trnas = [("tRNA_mock%d" % (i + 1), "".join(rng.choices(BASES, k=TRNA_LEN)))
             for i in range(N_TRNA)]
    mirnas = [("osa-MIRmock%d" % (i + 1),
               "".join(rng.choices(BASES, k=MIRNA_LEN)))
              for i in range(N_MIRNA)]
    return chroms, rrnas, trnas, mirnas


def make_read(chroms, rrnas, trnas, mirnas, rng):
    """One simulated SE small-RNA read (21-26 nt, optional 3' adapter tail).

    Body source mix: 55% miRNA (the full 21 nt sequence carrying 0-1
    substitution mismatches), 15% rRNA substring, 10% tRNA substring, 20%
    intergenic genome substring. 35% of the reads get a shifted 3' adapter
    tail so Trim Galore has realistic input."""
    u = rng.random()
    if u < MIRNA_FRAC:
        body = mirnas[rng.randrange(N_MIRNA)][1]
        for _ in range(rng.randrange(2)):          # 0-1 mismatches
            pos = rng.randrange(len(body))
            repl = rng.choice([b for b in BASES if b != body[pos]])
            body = body[:pos] + repl + body[pos + 1:]
    else:
        if u < MIRNA_FRAC + RRNA_FRAC:
            source = rrnas[rng.randrange(N_RRNA)][1]
        elif u < MIRNA_FRAC + RRNA_FRAC + TRNA_FRAC:
            source = trnas[rng.randrange(N_TRNA)][1]
        else:
            source = chroms["chr%d" % (rng.randrange(N_CHROM) + 1)]
        length = rng.randint(MIN_READ, MAX_READ)
        start = rng.randrange(len(source) - length + 1)
        body = source[start:start + length]
    if rng.random() < ADAPTER_FRAC:
        shift = rng.randrange(min(6, len(ADAPTER)))
        return body + ADAPTER[shift:]
    return body


def _gzip_text(path):
    """Open a deterministic gzip text writer (mtime pinned to 0, no stored
    file name -- byte-reproducible output)."""
    raw = open(path, "wb")
    gz = gzip.GzipFile(filename="", mode="wb", compresslevel=6,
                       fileobj=raw, mtime=0)
    return io.TextIOWrapper(gz, encoding="ascii", newline="\n")


def write_reference(outdir, chroms, rrnas, trnas, mirnas):
    ref_dir = os.path.join(outdir, "ref")
    os.makedirs(ref_dir, exist_ok=True)

    def write_fasta(name, seqs):
        with open(os.path.join(ref_dir, name), "w") as fh:
            for seq_name, seq in seqs:
                fh.write(">%s\n" % seq_name)
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + "\n")

    write_fasta("genome.fa", sorted(chroms.items()))
    write_fasta("rRNA.fa", rrnas)
    write_fasta("tRNA.fa", trnas)
    write_fasta("miRNA.fa", mirnas)


def write_fastqs(outdir, chroms, rrnas, trnas, mirnas, reads_per_sample, rng):
    """Draw the simulated SE reads for every sample from the shared stream."""
    raw_dir = os.path.join(outdir, "1.rawdata")
    os.makedirs(raw_dir, exist_ok=True)
    for sample in SAMPLES:
        fq = _gzip_text(os.path.join(raw_dir, "%s.fq.gz" % sample))
        try:
            for n in range(reads_per_sample):
                read = make_read(chroms, rrnas, trnas, mirnas, rng)
                fq.write("@r%07d 1:N:0:1\n%s\n+\n%s\n"
                         % (n, read, "I" * len(read)))
        finally:
            fq.close()


def write_samples(outdir):
    with open(os.path.join(outdir, "samples.csv"), "w") as fh:
        fh.write("sample_id\n")
        for sample in SAMPLES:
            fh.write(sample + "\n")


def write_config(outdir):
    with open(os.path.join(outdir, "config.yaml"), "w") as fh:
        fh.write(CONFIG_TEMPLATE)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir",
                    default=os.path.join(os.path.dirname(__file__), "data"))
    ap.add_argument("--reads", type=int, default=2000,
                    help="SE reads per sample")
    ap.add_argument("--seed", type=int, default=SEED)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(args.seed)
    chroms, rrnas, trnas, mirnas = build_reference(rng)
    write_reference(args.outdir, chroms, rrnas, trnas, mirnas)
    write_fastqs(args.outdir, chroms, rrnas, trnas, mirnas, args.reads, rng)
    write_samples(args.outdir)
    write_config(args.outdir)

    print("[make_testdata] chromosomes %d x %dbp, rRNA %d x %dbp, "
          "tRNA %d x %dbp, miRNA %d x %dnt"
          % (N_CHROM, CHROM_LEN, N_RRNA, RRNA_LEN,
             N_TRNA, TRNA_LEN, N_MIRNA, MIRNA_LEN))
    print("[make_testdata] samples %s x %d SE reads (21-26 nt, %d%% with a 3' "
          "adapter tail)"
          % (", ".join(SAMPLES), args.reads, round(ADAPTER_FRAC * 100)))
    print("[make_testdata] output root: %s" % os.path.abspath(args.outdir))


if __name__ == "__main__":
    main()
