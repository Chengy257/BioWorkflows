#!/usr/bin/env python3
"""Create the synthetic fixture tree consumed by the `test` / `test_full`
profiles. Fixed seed; the layout mirrors the Stage-2 rnaseq handoff
({aligner}/ with markdup BAMs) plus the reference set the chain consumes.

The fixtures are REAL-RUN capable (the deferred real-container round runs
the whole chain on them; see tests/build_fixture_dbs.sh for the container
side and tests/gate_planted_novel.py for the scientific gate):

- Read-carrying mini-BAMs (built with samtools from deterministic SAM
  text): spliced + exonic single-end reads over every planted locus and
  every reference gene, so StringTie really assembles transcripts and
  gffcompare really assigns classcodes.
- Designed genome content, not plain random sequence:
  * novA / novB (planted NONCODING): stop-codon-rich sequence in every
    frame (max ORF < 150 nt) - must survive every coding-potential filter;
  * novC (planted CODING-like): exon bodies are back-translated stop-free
    codons (a ~1366-codon ORF on the spliced transcript) - must be
    excluded by the coding-potential filters;
  * three extra ORF islands and three extra stop-rich blocks exist only as
    CPAT/CPC2 training material (no reads planted over them).
- novel.expected.gtf carries the planted loci with a planted_class
  attribute; the gate script consumes it (locus-based matching - StringTie
  assigns its own MSTRG ids).
- toy_proteins.fa / pfam_seed.sto are the raw material for the mini
  DIAMOND / Pfam databases (built in containers by build_fixture_dbs.sh).

Usage: SAMTOOLS_BIN=<samtools> python3 tests/make_fixtures.py
"""
import os
import random
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FX = HERE / "fixtures"
SAMPLES = ["w1", "w2", "m1", "m2"]
SAMPLE_SEEDS = {"w1": 101, "w2": 202, "m1": 303, "m2": 404}
CHROM, GENOME_LEN = "chr1", 300000
READ_LEN = 100

# Planted novel loci: two real exons + a real intron each, in intergenic
# windows (reference genes occupy [2000,5900], [62000,65900], [122000,
# 125900], [182000,185900], [242000,245900]).
NOVEL_LOCI = {
    # noncoding exon lengths are multiples of 3 so the spliced transcript
    # keeps one unbroken frame layout (the cassette design relies on it)
    "novA": {"exons": [(30000, 30800), (31400, 32200)], "strand": "+",
             "planted_class": "noncoding"},
    "novB": {"exons": [(80000, 81199), (82000, 83400)], "strand": "-",
             "planted_class": "noncoding"},
    "novC": {"exons": [(100000, 101999), (102600, 104699)], "strand": "+",
             "planted_class": "coding"},
}
# Training-only sequence islands (intergenic, read-free)
ORF_ISLANDS = [(140000, 144599), (150000, 154599), (160000, 164599)]   # coding
STOP_BLOCKS = [(40000, 40999), (110000, 110999), (130000, 130999)]     # noncoding

STOPS = ("TAA", "TAG", "TGA")
AA_CODON_TABLE = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT",
    "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAT",
    "P": "CCT", "Q": "CAA", "R": "CGT", "S": "TCT", "T": "ACT", "V": "GTT",
    "W": "TGG", "Y": "TAT",
}
CODON_TO_AA = {cod: aa for aa, cod in AA_CODON_TABLE.items()}
PROTEIN_REPEAT = "ACDEFGHIKLMNPQRSTVWY" * 70   # 1400 stop-free codons


def revcomp(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def coding_seq(length, rng):
    """Stop-free codon sequence in frame 1 (canonical codons, no wobble)."""
    n_codons = length // 3
    seq = "".join(AA_CODON_TABLE[aa]
                  for aa in PROTEIN_REPEAT[:n_codons] + "A" * (length % 3))
    return seq[:length]


# CPAT/CPC2 scan all six frames (both strands): a stop on the forward
# strand disappears on the reverse complement. The 66-nt cassette below
# carries a stop codon in each FORWARD frame (codons at 3/21/45) AND a
# reverse-complement stop (TTA/CTA/TCA) at codons 18/26/34, which maps to
# each REVERSE frame of the cassette - so no ORF in any of the six frames
# runs past ~66 nt, on either strand.
CASSETTE_LEN = 66
# forward frames 0/1/2 die at codons 3/22/44; the reverse-complement
# stops at 18/26/34 map to reverse frames 0/1/2 (a reverse-strand codon
# starting at fwd position p lives in reverse frame (63 - p) % 3)
CASSETTE_STOPS = [(3, "TAA"), (22, "TAG"), (44, "TGA"),   # forward frames
                  (18, "TTA"), (26, "CTA"), (34, "TCA")]  # reverse frames


def stop_rich_seq(length, rng):
    out = []
    made = 0
    while made < length:
        cassette = [None] * CASSETTE_LEN
        for pos, codon in CASSETTE_STOPS:
            cassette[pos:pos + 3] = list(codon)
        for i in range(CASSETTE_LEN):
            if cassette[i] is None:
                cassette[i] = rng.choice("ACGT")
        out.append("".join(cassette))
        made += CASSETTE_LEN
    return "".join(out)[:length]


def max_orf_len(seq):
    best = 0
    for strand_seq in (seq, revcomp(seq)):
        for frame in range(3):
            run = 0
            for i in range(frame, len(strand_seq) - 2, 3):
                if strand_seq[i:i + 3] in STOPS:
                    best = max(best, run)
                    run = 0
                else:
                    run += 3
            best = max(best, run)
    return best


def splice(locus, genome_seq):
    return "".join(genome_seq[s - 1:e] for s, e in locus["exons"])


def translate(seq):
    return "".join(CODON_TO_AA.get(seq[i:i + 3], "X")
                   for i in range(0, len(seq) - len(seq) % 3, 3))


def write_reference(ref, rng):
    """Genome: seeded random background, then designed islands written in.
    Returns the full genome sequence string."""
    genome = [rng.choice("ACGT") for _ in range(GENOME_LEN)]

    def implant(start, end, seq):
        genome[start - 1:end] = list(seq)

    for locus in NOVEL_LOCI.values():
        (s1, e1), (s2, e2) = locus["exons"]
        if locus["planted_class"] == "coding":
            need = sum(e - s + 1 for s, e in locus["exons"])
            cds = coding_seq(need, rng)
            pos = 0
            for s, e in locus["exons"]:
                implant(s, e, cds[pos:pos + e - s + 1])
                pos += e - s + 1
        else:
            implant(s1, e1, stop_rich_seq(e1 - s1 + 1, rng))
            implant(s2, e2, stop_rich_seq(e2 - s2 + 1, rng))
        implant(e1 + 1, s2 - 1, stop_rich_seq(s2 - 1 - e1, rng))
    for start, end in ORF_ISLANDS:
        implant(start, end, coding_seq(end - start + 1, rng))
    for start, end in STOP_BLOCKS:
        implant(start, end, stop_rich_seq(end - start + 1, rng))

    genome_seq = "".join(genome)
    with open(ref / "genome.fa", "w") as fh:
        fh.write(f">{CHROM}\n")
        for i in range(0, len(genome_seq), 60):
            fh.write(genome_seq[i:i + 60] + "\n")

    # reference annotation: 5 genes with 2 exons each (random bodies; they
    # are 'known' via gtf_pcg coordinate matching, not sequence)
    with open(ref / "genes.gtf", "w") as gtf, open(ref / "genes_pcg.gtf", "w") as pcg:
        gtf.write("#!genome-build test\n")
        pcg.write("#!genome-build test (protein-coding CDS only)\n")
        for k in range(5):
            start = 2000 + k * 60000
            gene = f"gene{k + 1:02d}"
            tx = f"tx{k + 1:02d}"
            for feat, s, e in (("exon", start, start + 900),
                               ("exon", start + 3000, start + 3900)):
                line = (f"{CHROM}\ttest\t{feat}\t{s}\t{e}\t.\t+\t.\t"
                        f"gene_id \"{gene}\"; transcript_id \"{tx}\";\n")
                gtf.write(line)
                pcg.write(line)
            tline = (f"{CHROM}\ttest\ttranscript\t{start}\t{start + 3900}\t.\t+\t.\t"
                     f"gene_id \"{gene}\"; transcript_id \"{tx}\";\n")
            gtf.write(tline)
            pcg.write(tline)

    # planted-loci manifest (the gate's ground truth)
    with open(ref / "novel.expected.gtf", "w") as gtf:
        for tid, locus in sorted(NOVEL_LOCI.items()):
            s = locus["exons"][0][0]
            e = locus["exons"][-1][1]
            for ex_s, ex_e in locus["exons"]:
                gtf.write(f"{CHROM}\tplanted\texon\t{ex_s}\t{ex_e}\t.\t"
                          f"{locus['strand']}\t.\t"
                          f"gene_id \"{tid}\"; transcript_id \"{tid}\";\n")
            gtf.write(f"{CHROM}\tplanted\ttranscript\t{s}\t{e}\t.\t"
                      f"{locus['strand']}\t.\t"
                      f"gene_id \"{tid}\"; transcript_id \"{tid}\"; "
                      f"planted_class \"{locus['planted_class']}\";\n")

    # CPAT training material (transcript-space sequences)
    with open(ref / "cpat_coding.fa", "w") as fh:
        seqs = [splice(NOVEL_LOCI["novC"], genome_seq)]
        seqs += [genome_seq[s - 1:e] for s, e in ORF_ISLANDS]
        for i, seq in enumerate(seqs):
            fh.write(f">coding{i}\n")
            for j in range(0, len(seq), 60):
                fh.write(seq[j:j + 60] + "\n")
    with open(ref / "cpat_noncoding.fa", "w") as fh:
        seqs = [splice(NOVEL_LOCI["novA"], genome_seq),
                splice(NOVEL_LOCI["novB"], genome_seq)]
        seqs += [genome_seq[s - 1:e] for s, e in STOP_BLOCKS]
        for i, seq in enumerate(seqs):
            fh.write(f">noncoding{i}\n")
            for j in range(0, len(seq), 60):
                fh.write(seq[j:j + 60] + "\n")

    # toy NR protein set: the novC ORF translation (a real homology hit for
    # the coding locus) + three unrelated proteins
    with open(ref / "toy_proteins.fa", "w") as fh:
        prot = translate(splice(NOVEL_LOCI["novC"], genome_seq))
        fh.write(">novC_orf\n")
        for j in range(0, len(prot), 60):
            fh.write(prot[j:j + 60] + "\n")
        rr = random.Random(99991)
        for k in range(3):
            p = "".join(rr.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(120))
            fh.write(f">prot{k}\n{p}\n")

    # Pfam seed alignment: a 60-aa synthetic domain family (3 sequences,
    # unrelated to any planted sequence - the filter must PASS novA/novB)
    with open(ref / "pfam_seed.sto", "w") as fh:
        rr = random.Random(777)
        fh.write("# STOCKHOLM 1.0\n\n")
        cols = [(rr.choice("ACDEFGHIKLMNPQRSTVWY"),
                 rr.choice("ACDEFGHIKLMNPQRSTVWY"),
                 rr.choice("ACDEFGHIKLMNPQRSTVWY")) for _ in range(60)]
        for name, idx in (("seedseq1", 0), ("seedseq2", 1), ("seedseq3", 2)):
            fh.write(f"{name} {''.join(c[idx] for c in cols)}\n")
        fh.write("//\n")
    return genome_seq


def sam_records(genome_seq, rng):
    """Deterministic single-end 100-nt reads over the planted loci and the
    reference genes. Junction reads are anchored exactly at the splice
    sites (<l>M<intron>N<r>M) and carry the XS strand tag - StringTie drops
    spliced reads without XS on unstranded input (the tag every real
    aligner emits); exonic reads tile the exon bodies every ~20 nt (the
    assembly-coverage threshold found during fixture validation)."""
    recs = []

    def add(pos0, cigar, name, xs=""):
        seq, p = "", pos0
        num = ""
        for ch in cigar:
            if ch.isdigit():
                num += ch
                continue
            count, num = int(num), ""
            if ch == "M":
                seq += genome_seq[p:p + count]
                p += count
            elif ch == "N":
                p += count
        recs.append((name, pos0, cigar, seq, xs))

    def tile(ex_s, ex_e, prefix, stride=20):
        made, pos = 0, ex_s - 1
        while pos + READ_LEN - 1 <= ex_e:
            add(pos, f"{READ_LEN}M", f"{prefix}_t{made}")
            made += 1
            pos += stride

    for tid, locus in sorted(NOVEL_LOCI.items()):
        (s1, e1), (s2, e2) = locus["exons"]
        gap = s2 - e1 - 1
        for left in range(35, 66):                 # every anchor position
            add(e1 - left, f"{left}M{gap}N{READ_LEN - left}M",
                f"{tid}_j{left}", xs=locus["strand"])
        tile(s1, e1, f"{tid}e1")
        tile(s2, e2, f"{tid}e2")

    for k in range(5):
        start = 2000 + k * 60000
        gene = f"gene{k + 1:02d}"
        tile(start, start + 900, f"{gene}e1")
        tile(start + 3000, start + 3900, f"{gene}e2")
        for left in range(35, 66):                 # intron 3000-900-1
            add(start + 900 - left, f"{left}M2099N{READ_LEN - left}M",
                f"{gene}_j{left}", xs="+")
    return recs


def write_bams(bam_dir, genome_seq):
    samtools = os.environ.get("SAMTOOLS_BIN", "samtools")
    for sample in SAMPLES:
        rng = random.Random(SAMPLE_SEEDS[sample])
        recs = sam_records(genome_seq, rng)
        sam_lines = [f"@HD\tVN:1.6\tSO:coordinate",
                     f"@SQ\tSN:{CHROM}\tLN:{GENOME_LEN}"]
        for name, pos0, cigar, seq, xs in sorted(recs, key=lambda r: r[1]):
            tag = f"\tXS:A:{xs}" if xs else ""
            sam_lines.append(
                f"{sample}|{name}\t0\t{CHROM}\t{pos0 + 1}\t60\t{cigar}\t"
                f"*\t0\t0\t{seq}\t{'I' * len(seq)}{tag}")
        sam_path = bam_dir / f"{sample}.sam"
        bam_path = bam_dir / f"{sample}.markdup.sorted.bam"
        with open(sam_path, "w") as fh:
            fh.write("\n".join(sam_lines) + "\n")
        subprocess.run([samtools, "view", "-bS", "-o", str(bam_path),
                        str(sam_path)], check=True)
        subprocess.run([samtools, "index", str(bam_path)], check=True)
        sam_path.unlink()


def main():
    rng = random.Random(20261012)
    if FX.exists():
        import shutil
        print(f"removing existing fixture tree {FX}")
        shutil.rmtree(FX)

    ref = FX / "reference"
    ref.mkdir(parents=True)
    bam_dir = FX / "results" / "star_salmon"
    bam_dir.mkdir(parents=True)

    (FX / "samplesheet.csv").write_text(
        "\n".join(["sample_id"] + SAMPLES) + "\n", encoding="utf-8")

    genome_seq = write_reference(ref, rng)
    write_bams(bam_dir, genome_seq)

    samtools = os.environ.get("SAMTOOLS_BIN", "samtools")
    subprocess.run([samtools, "faidx", str(ref / "genome.fa")], check=True)

    # the designed content must actually discriminate (fail loudly if not)
    for tid, locus in sorted(NOVEL_LOCI.items()):
        spliced = splice(locus, genome_seq)
        orf = max_orf_len(spliced)
        if locus["planted_class"] == "coding":
            ok = orf >= 1200
        else:
            ok = orf <= 150
        print(f"{tid}: spliced {len(spliced)} nt, max ORF {orf} nt "
              f"({locus['planted_class']}) -> {'OK' if ok else 'DESIGN FAILURE'}")
        if not ok:
            return 1
    print(f"fixtures written under {FX}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
