#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate the miniature seclip-seq regression test dataset.

Pure standard library implementation with a fixed random seed and the gzip
mtime pinned to 0, so outputs are byte-reproducible for a given version (the
data is not committed, it is generated on the fly on the test machine). A
single master RNG stream (seeded 42) builds the reference first and then draws
every read for every sample in order. Outputs (--outdir, default tests/data/):

    ref/genome.fa                2 x 20 kb chromosomes (chr1/chr2, random ACGT)
    ref/genes.gtf                10 genes per chromosome (gene + exon records,
                                 gene_id/transcript_id set)
    ref/repeats.fa               3 snRNA-like sequences x 300 bp
    1.rawdata/{sample}_R1.fq.gz  simulated SE reads: 10 N UMI + 30 bp
                                 genomic/repeat substring + shifted 3' adapter
                                 tail (AGATCGGAAGAGCAC + shifts); 60%
                                 genome-derived, 25% repeats-derived, 15%
                                 random; 1% substitution errors (default 2000
                                 reads per sample, --reads N)
    samples.csv                  sample table (sample_id column)
    config.yaml                  ready-to-run miniature configuration (relative
                                 ref/ paths, tiny STAR indices, CLIPper
                                 disabled -- the test environment has no
                                 CLIPper, so this also exercises the skip path)

Samples: FC_rep1, FC_rep2.

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
GENES_PER_CHROM = 10
N_REPEATS = 3
REPEAT_LEN = 300
READ_BODY = 30          # genomic/repeat-derived insert length
UMI_LEN = 10
ADAPTER = "AGATCGGAAGAGCAC"   # first of the 20 shifted 3' adapter variants
ERROR_RATE = 0.01       # substitution errors in genomic/repeat-derived bodies
GENOME_FRAC = 0.60      # genome-derived reads
REPEAT_FRAC = 0.25      # repeats-derived reads (random noise takes the rest)
SEED = 42
SAMPLES = ["FC_rep1", "FC_rep2"]

BASES = "ACGT"

# Miniature configuration written next to the generated data. Relative ref/
# paths; STAR indices shrunk for the 2 x 20 kb genome; CLIPper off because the
# test environment has no CLIPper install.
CONFIG_TEMPLATE = """\
SampleListFile: "samples.csv"
results_dir: "results"
threads: 2
species: "hsa"
filter_repeats: true
umi:
  pattern: "NNNNNNNNNN"
cutadapt:
  min_len: 18
  quality_cutoff: 6
  error_rate: 0.1
  adapters:
    - AGATCGGAAGAGCAC
    - GATCGGAAGAGCACA
    - ATCGGAAGAGCACAC
    - TCGGAAGAGCACACG
    - CGGAAGAGCACACGT
    - GGAAGAGCACACGTC
    - GAAGAGCACACGTCT
    - AAGAGCACACGTCTG
    - AGAGCACACGTCTGA
    - GAGCACACGTCTGAA
    - AGCACACGTCTGAAC
    - GCACACGTCTGAACT
    - CACACGTCTGAACTC
    - ACACGTCTGAACTCC
    - CACGTCTGAACTCCA
    - ACGTCTGAACTCCAG
    - CGTCTGAACTCCAGT
    - GTCTGAACTCCAGTC
    - TCTGAACTCCAGTCA
    - CTGAACTCCAGTCAC
star:
  sjdb_overhang: 49
  genome_sa_index_nbases: 4
  repeats_sa_index_nbases: 2
  repeats_limit_ram: 500000000000
  filter_multimap_nmax: 30
  align_multimap_nmax: 1
callpeak:
  pureclip: true
  clipper: false
  clipper_species: "GRCh38_v40"
genome: "ref/genome.fa"
gtf: "ref/genes.gtf"
repeats_fa: "ref/repeats.fa"
"""


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


def build_reference(rng):
    """Draw chromosomes and snRNA-like repeats from the master RNG stream."""
    chroms = {f"chr{i + 1}": "".join(rng.choices(BASES, k=CHROM_LEN))
              for i in range(N_CHROM)}
    repeats = [("snRNA_mock%d" % (i + 1), "".join(rng.choices(BASES, k=REPEAT_LEN)))
               for i in range(N_REPEATS)]
    genes = []
    for i in range(N_CHROM * GENES_PER_CHROM):
        chrom = "chr%d" % (i // GENES_PER_CHROM + 1)
        idx = i % GENES_PER_CHROM
        strand = "+" if i % 2 == 0 else "-"
        gid = "gene%03d" % (i + 1)
        tx_start = 500 + idx * 1800           # 0-based
        exons = [(tx_start, tx_start + 400),
                 (tx_start + 700, tx_start + 1100),
                 (tx_start + 1400, tx_start + 1700)]
        genes.append({"id": gid, "chrom": chrom, "strand": strand,
                      "tx_start": tx_start, "tx_end": exons[-1][1],
                      "exons": exons})
    return chroms, genes, repeats


def make_read(chroms, repeats, rng):
    """One simulated SE read: UMI + body(+errors) + shifted adapter prefix.

    Body source mix: 60% genome substring, 25% repeats substring, 15% random
    noise. Every read keeps the same 10 N UMI + 30 bp body + adapter tail
    structure so umi_tools extract and cutadapt both see realistic input."""
    umi = "".join(rng.choice(BASES) for _ in range(UMI_LEN))
    u = rng.random()
    if u < GENOME_FRAC:
        source = chroms["chr%d" % (rng.randrange(N_CHROM) + 1)]
        start = rng.randrange(len(source) - READ_BODY)
        body = mutate(source[start:start + READ_BODY], rng)
    elif u < GENOME_FRAC + REPEAT_FRAC:
        source = repeats[rng.randrange(N_REPEATS)][1]
        start = rng.randrange(len(source) - READ_BODY)
        body = mutate(source[start:start + READ_BODY], rng)
    else:
        body = "".join(rng.choice(BASES) for _ in range(READ_BODY))
    shift = rng.randrange(min(6, len(ADAPTER)))
    return umi + body + ADAPTER[shift:]


def _gzip_text(path):
    """Open a deterministic gzip text writer (mtime pinned to 0, no stored
    file name -- byte-reproducible output)."""
    raw = open(path, "wb")
    gz = gzip.GzipFile(filename="", mode="wb", compresslevel=6,
                       fileobj=raw, mtime=0)
    return io.TextIOWrapper(gz, encoding="ascii", newline="\n")


def write_reference(outdir, chroms, genes, repeats):
    ref_dir = os.path.join(outdir, "ref")
    os.makedirs(ref_dir, exist_ok=True)

    with open(os.path.join(ref_dir, "genome.fa"), "w") as fh:
        for c in sorted(chroms):
            fh.write(">%s\n" % c)
            s = chroms[c]
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")

    with open(os.path.join(ref_dir, "genes.gtf"), "w") as fh:
        for g in genes:
            attr = 'gene_id "%s"; transcript_id "%s.1"' % (g["id"], g["id"])
            fh.write("%s\tsim\tgene\t%d\t%d\t.\t%s\t.\t%s;\n"
                     % (g["chrom"], g["tx_start"] + 1, g["tx_end"],
                        g["strand"], attr))
            for s, e in g["exons"]:
                fh.write("%s\tsim\texon\t%d\t%d\t.\t%s\t.\t%s;\n"
                         % (g["chrom"], s + 1, e, g["strand"], attr))

    with open(os.path.join(ref_dir, "repeats.fa"), "w") as fh:
        for name, seq in repeats:
            fh.write(">%s\n" % name)
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def write_fastqs(outdir, chroms, repeats, reads_per_sample, rng):
    """Draw the simulated SE reads for every sample from the shared stream."""
    raw_dir = os.path.join(outdir, "1.rawdata")
    os.makedirs(raw_dir, exist_ok=True)
    for sample in SAMPLES:
        fq = _gzip_text(os.path.join(raw_dir, "%s_R1.fq.gz" % sample))
        try:
            for n in range(reads_per_sample):
                read = make_read(chroms, repeats, rng)
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
    chroms, genes, repeats = build_reference(rng)
    write_reference(args.outdir, chroms, genes, repeats)
    write_fastqs(args.outdir, chroms, repeats, args.reads, rng)
    write_samples(args.outdir)
    write_config(args.outdir)

    print("[make_testdata] chromosomes %d x %dbp, genes %d, repeats %d x %dbp"
          % (N_CHROM, CHROM_LEN, len(genes), N_REPEATS, REPEAT_LEN))
    print("[make_testdata] samples %s x %d SE reads"
          % (", ".join(SAMPLES), args.reads))
    print("[make_testdata] output root: %s" % os.path.abspath(args.outdir))


if __name__ == "__main__":
    main()
