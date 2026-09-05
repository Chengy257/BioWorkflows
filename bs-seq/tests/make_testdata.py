#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate the miniature bs-seq regression test dataset.

Pure standard library implementation with a fixed random seed and the gzip
mtime pinned to 0, so outputs are byte-reproducible for a given version (the
data is not committed, it is generated on the fly on the test machine). A
single master RNG stream (seeded 42) builds the genome first and then draws
every fragment for every sample in order. Outputs (--outdir, default
tests/data/):

    ref/genome.fa                    2 x 20 kb chromosomes (chr1/chr2, random
                                     ACGT)
    1.rawdata/{sample}_1.fastq.gz    simulated PE 100 bp reads, post-bisulfite
    1.rawdata/{sample}_2.fastq.gz    (default 3000 pairs per sample, --reads N):
                                     150-300 bp fragments are drawn from the
                                     genome; read 1 sequences the fragment
                                     start on the forward strand with C->T
                                     conversion, read 2 sequences the fragment
                                     end from the opposite strand with G->A
                                     conversion (the directional Bismark read
                                     layout); 1% sequencing errors on top of
                                     the converted read; 30% of the pairs carry
                                     a 3' adapter tail (AGATCGGAAGAGC plus a
                                     few extra bases) appended past the 100 bp
                                     genomic body, so adapter trimming
                                     recovers an alignable read
    samples.csv                      sample table (sample_id column)
    config.yaml                      ready-to-run miniature configuration
                                     (relative ref/ path, full repository key
                                     set)

Samples: s1, s2.

Usage:
    python3 make_testdata.py [--outdir tests/data] [--reads 3000] [--seed 42]
"""
import argparse
import gzip
import io
import os
import random

CHROM_LEN = 20000
N_CHROM = 2
READ_LEN = 100            # simulated PE read length
MIN_FRAG = 150            # simulated insert size range
MAX_FRAG = 300
ERROR_RATE = 0.01         # sequencing error rate on the converted read
ADAPTER_FRAC = 0.30       # fraction of pairs carrying a 3' adapter tail
ADAPTER = "AGATCGGAAGAGC"   # trim.adapter from config.yaml
TAIL_MIN = 2              # extra bases beyond the fixed 13-mer adapter core
TAIL_MAX = 8
SEED = 42
SAMPLES = ["s1", "s2"]

BASES = "ACGT"
_REVCOMP = str.maketrans("ACGT", "TGCA")

# Miniature configuration written next to the generated data. Relative ref/
# path; the full repository key set at reduced scale (threads capped at 2)
# so the dry-run regression exercises the same schema as config/config.yaml.
CONFIG_TEMPLATE = """\
SampleListFile: "samples.csv"
results_dir: "results"
threads: 2
species: "none"
trim:
  enabled: true
  quality: 20
  min_len: 20
  adapter: "AGATCGGAAGAGC"
  stringency: 3
  error_rate: 0.1
  extra: ""
bismark:
  align_extra: "--phred33-quals"
methylation_extractor:
  cx_report: false
  merge_cpg: true
  buffer_frac: 4
genome: "ref/genome.fa"
"""


def revcomp(seq):
    """Reverse complement of an ACGT string."""
    return seq.translate(_REVCOMP)[::-1]


def convert(seq, src, dst):
    """Bisulfite conversion of one strand: src -> dst in place."""
    return seq.replace(src, dst)


def add_errors(seq, rng):
    """Apply 1% substitution errors to the converted read."""
    out = []
    for base in seq:
        if rng.random() < ERROR_RATE:
            out.append(rng.choice([b for b in BASES if b != base]))
        else:
            out.append(base)
    return "".join(out)


def build_reference(rng):
    """Draw the chromosomes from the master RNG stream."""
    chroms = {f"chr{i + 1}": "".join(rng.choices(BASES, k=CHROM_LEN))
              for i in range(N_CHROM)}
    return chroms


def make_pair(chroms, rng):
    """One simulated PE read pair (100 bp, post-bisulfite).

    A 150-300 bp fragment is drawn from the genome. Read 1 is the C->T-
    converted forward substrate of the fragment start; read 2 is the G->A-
    converted reverse-complement substrate of the fragment end. Sequencing
    errors (1%) are applied on top of the converted read. 30% of the pairs
    get the same adapter tail appended to both reads, always past the 100 bp
    genomic body, so trimming restores the full alignable insert ends."""
    chrom = chroms["chr%d" % (rng.randrange(N_CHROM) + 1)]
    frag_len = rng.randint(MIN_FRAG, MAX_FRAG)
    start = rng.randrange(CHROM_LEN - frag_len + 1)
    fragment = chrom[start:start + frag_len]
    r1 = convert(fragment[:READ_LEN], "C", "T")
    r2 = convert(revcomp(fragment[-READ_LEN:]), "G", "A")
    r1 = add_errors(r1, rng)
    r2 = add_errors(r2, rng)
    if rng.random() < ADAPTER_FRAC:
        tail = ADAPTER + "".join(
            rng.choices(BASES, k=rng.randint(TAIL_MIN, TAIL_MAX)))
        r1 += tail
        r2 += tail
    return r1, r2


def _gzip_text(path):
    """Open a deterministic gzip text writer (mtime pinned to 0, no stored
    file name -- byte-reproducible output)."""
    raw = open(path, "wb")
    gz = gzip.GzipFile(filename="", mode="wb", compresslevel=6,
                       fileobj=raw, mtime=0)
    return io.TextIOWrapper(gz, encoding="ascii", newline="\n")


def write_reference(outdir, chroms):
    ref_dir = os.path.join(outdir, "ref")
    os.makedirs(ref_dir, exist_ok=True)
    with open(os.path.join(ref_dir, "genome.fa"), "w") as fh:
        for seq_name, seq in sorted(chroms.items()):
            fh.write(">%s\n" % seq_name)
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def write_fastqs(outdir, chroms, pairs_per_sample, rng):
    """Draw the simulated PE reads for every sample from the shared stream."""
    raw_dir = os.path.join(outdir, "1.rawdata")
    os.makedirs(raw_dir, exist_ok=True)
    for sample in SAMPLES:
        fq1 = _gzip_text(os.path.join(raw_dir, "%s_1.fastq.gz" % sample))
        fq2 = _gzip_text(os.path.join(raw_dir, "%s_2.fastq.gz" % sample))
        try:
            for n in range(pairs_per_sample):
                r1, r2 = make_pair(chroms, rng)
                fq1.write("@r%07d 1:N:0:1\n%s\n+\n%s\n"
                          % (n, r1, "I" * len(r1)))
                fq2.write("@r%07d 2:N:0:1\n%s\n+\n%s\n"
                          % (n, r2, "I" * len(r2)))
        finally:
            fq1.close()
            fq2.close()


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
    ap.add_argument("--reads", type=int, default=3000,
                    help="PE read pairs per sample")
    ap.add_argument("--seed", type=int, default=SEED)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(args.seed)
    chroms = build_reference(rng)
    write_reference(args.outdir, chroms)
    write_fastqs(args.outdir, chroms, args.reads, rng)
    write_samples(args.outdir)
    write_config(args.outdir)

    print("[make_testdata] chromosomes %d x %dbp, PE %dbp reads, fragments %d-%dbp"
          % (N_CHROM, CHROM_LEN, READ_LEN, MIN_FRAG, MAX_FRAG))
    print("[make_testdata] samples %s x %d PE pairs (1%% errors, %d%% with a 3' "
          "adapter tail)"
          % (", ".join(SAMPLES), args.reads, round(ADAPTER_FRAC * 100)))
    print("[make_testdata] output root: %s" % os.path.abspath(args.outdir))


if __name__ == "__main__":
    main()
