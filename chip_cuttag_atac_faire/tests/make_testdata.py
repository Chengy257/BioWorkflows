#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate the miniature ChIP/CUT&Tag/ATAC regression test dataset (shared by the
dry-run regression and --real-run validation).

Pure standard library, fixed random seed (default 42); two runs with the same
arguments are byte-identical (gzip header mtime pinned to 0). Data is not
committed to the repository; it is generated on the test machine on demand.
Outputs (written to the --outdir directory; paths follow the data working-directory layout):

    ref/genome.fa             2 x 100 kb chromosomes (chr1/chr2)
    ref/genes.gtf             30 mock genes per chromosome (gene + exon levels;
                              exons carry gene_id/transcript_id so ChIPseeker accepts them)
    ref/genes.bed             BED6 gene models (chrom/start/end/name/score/strand)
    ref/chrom.sizes           chromosome length table (for bigWig)
    1.rawdata/{sample}_1.fq.gz and _2.fq.gz
                              simulated PE reads (length 50, fragments 150-300bp,
                              sequences are real substrings of the reference genome,
                              so --real-run passes bowtie2)
    samples.csv               5-sample, 6-column sample table (sample_id,role,group,seqtype,layout,peak_type)
    config.yaml               test config (relative paths, threads: 2, all three QC switches on)

Samples: chip_treat_rep1/rep2 + chip_control (narrow group g1),
         atac_treat_rep1/rep2 (atac group g2).
Enrichment design: three 2kb peak regions pre-seeded on chr1; treat samples draw
70% of fragments from peak regions, control samples sample uniformly genome-wide.

Usage:
    python tests/make_testdata.py --outdir <dir> [--reads 50000] [--seed 42]
        [--replicate] [--qc-full] [--motif] [--diffbind] [--gates] [--spike-in]
"""
import argparse
import gzip
import io
import os
import random

# ---------------------------------------------------------------------
# Fixed parameters (matches implementation plan Task 4.2; update the
# run_tests.py assertions in lockstep before changing)
# ---------------------------------------------------------------------
CHROM_LEN = 100000        # per-chromosome length (bp)
N_CHROM = 2               # number of chromosomes (chr1/chr2)
SPIKE_CONTIGS = 2         # spike-in contigs (--spike-in scenario)
SPIKE_CONTIG_LEN = 3000   # spike-in per-contig length (bp)
GENES_PER_CHROM = 30      # genes per chromosome
READ_LEN = 50             # read length
FRAG_MIN, FRAG_MAX = 150, 300   # fragment length range (bp)
GENE_START = 1000         # first gene start (0-based)
GENE_SPACING = 3000       # gene spacing
GENE_LEN = 2000           # gene span
PEAK_REGIONS = [(20000, 22000), (50000, 52000), (80000, 82000)]  # chr1 enrichment peaks (0-based half-open)
PEAK_PROB = 0.7           # probability a treat fragment comes from a peak region (control samples uniformly)
ERROR_RATE = 0.005        # substitution-type sequencing error rate (a few errors avoid pathological exact duplicates)
SEED = 42                 # default random seed

# Sample set: (sample_id, role, group, seqtype, peak_type)
# Names and group names all satisfy the workflow's _NAME_RE (alphanumerics plus
# . _ - , no consecutive underscores __)
SAMPLES = [
    ("chip_treat_rep1", "treat",   "g1", "chip", "narrow"),
    ("chip_treat_rep2", "treat",   "g1", "chip", "narrow"),
    ("chip_control",    "control", "g1", "chip", "narrow"),
    ("atac_treat_rep1", "treat",   "g2", "atac", "none"),
    ("atac_treat_rep2", "treat",   "g2", "atac", "none"),
]

# Extra group for the --replicate scenario: a broad chip group with two
# treats (exercises callpeak_broad_replicate + multiinter consensus on top
# of the narrow IDR already covered by g1/g2 being two-treat groups).
BROAD_SAMPLES = [
    ("hist_treat_rep1", "treat",   "g3", "chip", "broad"),
    ("hist_treat_rep2", "treat",   "g3", "chip", "broad"),
    ("hist_control",    "control", "g3", "chip", "broad"),
]

# Extra group for the --diffbind scenario: a second two-treat narrow chip
# group so one contrast ([g1, g4]) exists with replicate structure on both
# arms; carries the optional condition/batch columns.
DIFFBIND_SAMPLES = [
    ("db_treat_rep1", "treat",   "g4", "chip", "narrow", "mutant", "b1"),
    ("db_treat_rep2", "treat",   "g4", "chip", "narrow", "mutant", "b2"),
    ("db_control",    "control", "g4", "chip", "narrow", "",       ""),
]

# Sample-table header variants: the base 6 columns, and the extended schema
# with the optional condition/batch columns (written only when the
# --diffbind group set is present).
BASE_HEADER = "sample_id,role,group,seqtype,layout,peak_type"
EXT_HEADER = BASE_HEADER + ",condition,batch"

# Config insertions for the scenario flags (see write_config): valid YAML
# sub-blocks spliced into the base CONFIG_YAML (YAML forbids duplicate
# top-level keys, so the blocks are inserted under peak:/qc: instead of
# redefining them).
PEAK_REPLICATE_TAIL = """\
  replicate:
    enabled: true
    qvalue: 0.01
    idr_threshold: 0.05
    idr_rank: "p.value"
    consensus_min_replicates: 2
    frip_on: "pooled"
"""

QC_FULL_TAIL = """\
  tss: true
  organelle: true
  organelle_patterns: [chrc, chrm]
"""

QC_FULL_BLACKLIST = '\n# ---------- Extended QC (--qc-full scenario) ----------\nblacklist: "ref/blacklist.bed"\n'

GATES_TAIL = """\
  gates:
    enabled: true
    thresholds:
      mapping_rate_min: 0.70
      dup_rate_max: 0.50
      frip_min: 0.01
      nsc_min: 1.05
      rsc_min: 0.8
      tss_min: 6.0
      organelle_max: 0.20
"""

MOTIF_CONFIG_BLOCK = """
# ---------- Motif enrichment (--motif scenario) ----------
motif:
  enabled: true
  homer_genome: "test_genome"
  size: "given"
  background: ""
  extra: ""
"""

DIFFBIND_CONFIG_BLOCK = """
# ---------- Differential binding (--diffbind scenario) ----------
diffbind:
  enabled: true
  contrasts: [["g1", "g4"]]
  analysis: "DESeq2"
  summit_flank: 250
  use_controls: false
  fdr: 0.05
  foldchange: 1.0
  batch_correction: true
"""

SPIKE_CONFIG_BLOCK = """
# ---------- Spike-in normalization (--spike-in scenario) ----------
spike_in:
  enabled: true
  fasta: "ref/spike.fa"
  name: "lambda"
  scale_bigwigs: true
"""

# Synthetic blacklist for --qc-full: overlaps the first pre-seeded peak
# region so a real run can observe peaks being removed.
BLACKLIST_BED = "chr1\t19800\t20600\n"

BASES = "ACGT"
_COMP = str.maketrans("ACGT", "TGCA")

# Test config: the key set covers every key required by validate_config; paths
# are relative to the data working directory.
# qc.nsc_rsc=true only builds the DAG in the CI dry-run (spp is not executed), so it is safe.
CONFIG_YAML = """\
# =====================================================================
# Config dedicated to the synthetic test dataset (generated by
# tests/make_testdata.py; do not edit by hand).
# All paths are relative to the data working directory.
# Usage: bash run.sh -P <workdir> -c <this file> -n
# =====================================================================

# ---------- Reference genome (synthetic 2 x 100kb) ----------
genome_fa: "ref/genome.fa"     # bowtie2 index input
gtf: "ref/genes.gtf"           # for peak annotation (parseable by ChIPseeker makeTxDbFromGFF)
bed: "ref/genes.bed"           # reserved (used by deeptools)
chromsize: "ref/chrom.sizes"   # for bigWig generation
genome_size: "180000"          # MACS2 effective genome size: two 100kb chromosomes

# ---------- Sample table and resources ----------
grouplist: "samples.csv"       # sample table written next to the generator outputs
threads: 2

# ---------- Quality trimming (trim_galore) ----------
trim:
  quality: 25
  stringency: 3
  error_rate: 0.1
  extra: ""

# ---------- Alignment ----------
bowtie2_extra: "--end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700"
min_mapq: 30

# ---------- Deduplication ----------
dedup:
  chip: true
  cuttag: false
  atac: true
  faire: true

# ---------- Peak calling ----------
peak:
  keepdup: all
  qvalue: 0.05
  broad_cutoff: 0.05
  atac:
    mode: bampe
    shift: -100
    extsize: 200

# ---------- Annotation and signal windows ----------
region_flank: 3000

# ---------- QC (all three switches on, so the dry-run covers every QC branch) ----------
qc:
  nsc_rsc: true    # safe: dry-run does not execute spp; --real-run needs phantompeakqualtools
  frip: true
  deeptools: true
"""


def revcomp(s):
    """Reverse complement (R2 is taken from the other end of the fragment)."""
    return s.translate(_COMP)[::-1]


def mutate(seq, rng):
    """Introduce substitution-type sequencing errors at ERROR_RATE (keeps sequences mappable to the reference)."""
    if ERROR_RATE <= 0:
        return seq
    out = []
    for b in seq:
        if rng.random() < ERROR_RATE:
            out.append(rng.choice(BASES))
        else:
            out.append(b)
    return "".join(out)


def gzip_text(path):
    """gzip text handle with a fixed mtime (guarantees byte-identical output across runs with the same args; LF line endings)."""
    raw = gzip.GzipFile(path, "wb", compresslevel=6, mtime=0)
    return io.TextIOWrapper(raw, encoding="ascii", newline="\n")


def build_reference(seed):
    """Return (chromosome-sequence dict, gene-structure list).

    Chromosome sequences come from the fixed seed; gene positions use a
    deterministic layout (equidistant, alternating strands) that consumes no
    random numbers — so genome.fa is independent of --reads and reproducible on its own.
    """
    rng = random.Random(seed)
    chroms = {f"chr{i + 1}": "".join(rng.choices(BASES, k=CHROM_LEN))
              for i in range(N_CHROM)}
    genes = []
    for i in range(N_CHROM * GENES_PER_CHROM):
        chrom = f"chr{i // GENES_PER_CHROM + 1}"
        idx = i % GENES_PER_CHROM
        strand = "+" if i % 2 == 0 else "-"
        gid = f"gene{i + 1}"
        s0 = GENE_START + idx * GENE_SPACING          # gene start (0-based)
        exons = [(s0, s0 + 800), (s0 + 1200, s0 + GENE_LEN)]  # two exons (0-based half-open)
        genes.append({"id": gid, "chrom": chrom, "strand": strand,
                      "start": s0, "end": s0 + GENE_LEN, "exons": exons})
    return chroms, genes


def write_reference(outdir, chroms, genes):
    """Write the reference file set ref/{genome.fa, genes.gtf, genes.bed, chrom.sizes} (LF line endings)."""
    ref_dir = os.path.join(outdir, "ref")
    os.makedirs(ref_dir, exist_ok=True)

    with open(os.path.join(ref_dir, "genome.fa"), "w",
              encoding="ascii", newline="\n") as fh:
        for c in sorted(chroms):
            fh.write(f">{c}\n")
            s = chroms[c]
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")

    # GTF: minimal gene + exon two-level records; exons carry transcript_id for
    # makeTxDbFromGFF grouping
    with open(os.path.join(ref_dir, "genes.gtf"), "w",
              encoding="ascii", newline="\n") as fh:
        for g in genes:
            attr = f'gene_id "{g["id"]}";'
            fh.write(f'{g["chrom"]}\ttest\tgene\t{g["start"] + 1}\t{g["end"]}\t'
                     f'.\t{g["strand"]}\t.\t{attr}\n')
            exon_attr = attr + f' transcript_id "{g["id"]}.1";'
            for s, e in g["exons"]:
                fh.write(f'{g["chrom"]}\ttest\texon\t{s + 1}\t{e}\t'
                         f'.\t{g["strand"]}\t.\t{exon_attr}\n')

    # BED6: chrom/start/end/name/score/strand (ChIPseeker/rtracklayer compatible)
    with open(os.path.join(ref_dir, "genes.bed"), "w",
              encoding="ascii", newline="\n") as fh:
        for g in genes:
            fh.write(f'{g["chrom"]}\t{g["start"]}\t{g["end"]}\t'
                     f'{g["id"]}\t0\t{g["strand"]}\n')

    with open(os.path.join(ref_dir, "chrom.sizes"), "w",
              encoding="ascii", newline="\n") as fh:
        for c in sorted(chroms):
            fh.write(f"{c}\t{len(chroms[c])}\n")


def sample_fragment(chroms, rng, role):
    """Sample one genomic fragment by role: treat draws from a chr1 peak region
    with probability PEAK_PROB; everything else (including control) samples the
    genome uniformly. The sequence is always a real substring of the reference."""
    frag_len = rng.randint(FRAG_MIN, FRAG_MAX)
    if role == "treat" and rng.random() < PEAK_PROB:
        chrom = "chr1"
        rs, re_ = rng.choice(PEAK_REGIONS)
        start = rng.randint(rs, re_ - frag_len)   # fragment fully inside the peak region
    else:
        chrom = f"chr{rng.randint(1, N_CHROM)}"
        start = rng.randint(0, CHROM_LEN - frag_len)
    return chroms[chrom][start:start + frag_len]


def write_fastqs(outdir, chroms, genes, reads_per_sample, seed, samples):
    """Generate PE reads per sample (1.rawdata/{sample}_1.fq.gz and _2.fq.gz).

    Each sample uses its own rng (seed + sample index), decoupled from sample
    generation order; R2 is the reverse complement of the fragment's other end,
    quality lines are a fixed repeat of 'I'.
    """
    raw_dir = os.path.join(outdir, "1.rawdata")
    os.makedirs(raw_dir, exist_ok=True)

    for k, (sid, role, *_rest) in enumerate(samples):
        rng = random.Random(seed + k)
        fq1 = gzip_text(os.path.join(raw_dir, f"{sid}_1.fq.gz"))
        fq2 = gzip_text(os.path.join(raw_dir, f"{sid}_2.fq.gz"))
        try:
            for n in range(reads_per_sample):
                frag = sample_fragment(chroms, rng, role)
                r1 = mutate(frag[:READ_LEN], rng)
                r2 = revcomp(mutate(frag[-READ_LEN:], rng))
                q = "I" * READ_LEN
                fq1.write(f"@r{n:07d} 1:N:0:1\n{r1}\n+\n{q}\n")
                fq2.write(f"@r{n:07d} 2:N:0:1\n{r2}\n+\n{q}\n")
        finally:
            fq1.close()
            fq2.close()


def _normalize(rows):
    """Extend 5-tuples to the 7-field form (condition/batch default empty)."""
    out = []
    for row in rows:
        sid, role, grp, seqtype, pt = row[:5]
        cond = row[5] if len(row) > 5 else ""
        bat = row[6] if len(row) > 6 else ""
        out.append((sid, role, grp, seqtype, pt, cond, bat))
    return out


def write_samples(outdir, samples, extended):
    """Write the sample table samples.csv: the base 6 columns, plus the
    optional condition/batch columns when the diffbind group set is present
    (header matches the workflow's REQUIRED_COLUMNS + optional columns)."""
    path = os.path.join(outdir, "samples.csv")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write((EXT_HEADER if extended else BASE_HEADER) + "\n")
        for sid, role, grp, seqtype, pt, cond, bat in samples:
            if extended:
                fh.write(f"{sid},{role},{grp},{seqtype},PE,{pt},{cond},{bat}\n")
            else:
                fh.write(f"{sid},{role},{grp},{seqtype},PE,{pt}\n")


def write_blacklist(outdir):
    """Write the --qc-full synthetic blacklist (overlaps the first peak region)."""
    path = os.path.join(outdir, "ref", "blacklist.bed")
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write(BLACKLIST_BED)


def write_spike_reference(outdir, seed):
    """Write the --spike-in synthetic spike-in genome ref/spike.fa (2 x 3kb
    contigs, generated from its own seeded rng so the sequences are distinct
    from the main reference)."""
    rng = random.Random(seed + 101)
    path = os.path.join(outdir, "ref", "spike.fa")
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        for i in range(SPIKE_CONTIGS):
            fh.write(f">spike_ctg{i + 1}\n")
            s = "".join(rng.choices(BASES, k=SPIKE_CONTIG_LEN))
            for j in range(0, len(s), 60):
                fh.write(s[j:j + 60] + "\n")


def write_config(outdir, replicate=False, qc_full=False, motif=False, diffbind=False,
                 gates=False, spike=False):
    """Write the test config.yaml (relative paths, consumed via run.sh -c).

    Scenario flags splice the matching sub-blocks into the base config."""
    text = CONFIG_YAML
    if replicate:
        anchor = "    extsize: 200\n"
        assert anchor in text, "base config anchor for the replicate block moved"
        text = text.replace(anchor, anchor + PEAK_REPLICATE_TAIL, 1)
    if qc_full:
        anchor = "  deeptools: true\n"
        assert anchor in text, "base config anchor for the qc-full block moved"
        text = text.replace(anchor, anchor + QC_FULL_TAIL, 1)
        text += QC_FULL_BLACKLIST
    if gates:
        anchor = "  deeptools: true\n"
        assert anchor in text, "base config anchor for the gates block moved"
        text = text.replace(anchor, anchor + GATES_TAIL, 1)
    if motif:
        text += MOTIF_CONFIG_BLOCK
    if diffbind:
        text += DIFFBIND_CONFIG_BLOCK
    if spike:
        text += SPIKE_CONFIG_BLOCK
    path = os.path.join(outdir, "config.yaml")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True,
                    help="output root directory (required; reused if it already exists)")
    ap.add_argument("--reads", type=int, default=50000,
                    help="PE read pairs per sample (default 50000; CI regression uses 2000)")
    ap.add_argument("--seed", type=int, default=SEED,
                    help=f"random seed (default {SEED}, for reproducibility)")
    ap.add_argument("--replicate", action="store_true",
                    help="add a 2-treat broad chip group (g3) and enable the "
                         "replicate-aware peak stage (IDR + consensus scenario)")
    ap.add_argument("--qc-full", action="store_true", dest="qc_full",
                    help="enable tss/organelle QC and the synthetic blacklist "
                         "(extended QC scenario)")
    ap.add_argument("--motif", action="store_true",
                    help="enable the HOMER motif stage with a dummy genome tag "
                         "(dry-run only executes the DAG)")
    ap.add_argument("--diffbind", action="store_true",
                    help="add a 2-treat narrow chip group (g4) with condition/"
                         "batch columns and enable one DiffBind contrast")
    ap.add_argument("--gates", action="store_true",
                    help="enable the QC gate summary stage (qc.gates) with the "
                         "default thresholds")
    ap.add_argument("--spike-in", action="store_true", dest="spike_in",
                    help="write a synthetic spike-in reference (ref/spike.fa) and "
                         "enable the spike_in stage with scaled bigWigs")
    args = ap.parse_args()
    if args.reads < 1:
        ap.error("--reads must be a positive integer")

    samples = _normalize(SAMPLES)
    if args.replicate:
        samples += _normalize(BROAD_SAMPLES)
    if args.diffbind:
        samples += _normalize(DIFFBIND_SAMPLES)
    os.makedirs(args.outdir, exist_ok=True)
    chroms, genes = build_reference(args.seed)
    write_reference(args.outdir, chroms, genes)
    write_fastqs(args.outdir, chroms, genes, args.reads, args.seed, samples)
    write_samples(args.outdir, samples, extended=args.diffbind)
    write_config(args.outdir, replicate=args.replicate, qc_full=args.qc_full,
                 motif=args.motif, diffbind=args.diffbind, gates=args.gates,
                 spike=args.spike_in)
    if args.qc_full:
        write_blacklist(args.outdir)
    if args.spike_in:
        write_spike_reference(args.outdir, args.seed)

    tags = [t for t, on in (("+replicate", args.replicate),
                            ("+qc-full", args.qc_full),
                            ("+motif", args.motif),
                            ("+diffbind", args.diffbind),
                            ("+gates", args.gates),
                            ("+spike-in", args.spike_in)) if on]
    print(f"[make_testdata] chromosomes {N_CHROM} x {CHROM_LEN}bp, {len(genes)} genes, "
          f"{len(samples)} samples x {args.reads} PE read pairs (seed={args.seed}"
          + (", " + ", ".join(tags) if tags else "") + ")")
    print(f"[make_testdata] output root: {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()
