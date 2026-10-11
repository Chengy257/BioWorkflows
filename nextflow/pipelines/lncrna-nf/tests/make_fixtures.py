#!/usr/bin/env python3
"""Create the synthetic fixture tree consumed by the `test` / `test_full`
profiles. Fixed seed; the layout mirrors the Stage-2 rnaseq handoff
({aligner}/ with markdup BAMs) plus the reference set the chain consumes.

Scientific gate (validated at the deferred real-container round): planted
NONCING transcripts must be recovered by the final set and planted CODING
transcripts must be excluded by the coding-potential + homology filters.

    fixtures/
    ├── samplesheet.csv                 (4 samples, single sample_id column)
    ├── reference/
    │   ├── genome.fa(.fai)             300 kb chr1
    │   ├── genes.gtf                   reference annotation (5 known genes)
    │   ├── genes_pcg.gtf               protein-coding-only CDS GTF
    │   ├── cpat_model.logit            placeholder CPAT logit model
    │   ├── hexamer.tsv                 placeholder hexamer table
    │   ├── pfam_db/Pfam-A.hmm          mini Pfam (empty-ish HMM library)
    │   └── nr.dmnd                     mini DIAMOND database (placeholder)
    └── results/star_salmon/
        └── {id}.markdup.sorted.bam(+.bai)   (placeholder; dry-run only)

Usage: python3 tests/make_fixtures.py
"""
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FX = HERE / "fixtures"
SAMPLES = ["w1", "w2", "m1", "m2"]
# three novel loci on chr1 (intergenic between the reference genes)
NOVEL_LOCI = {
    "novA": (30000, 33600, "+"),   # planted NONCODING: short ORF, no homology
    "novB": (60000, 65000, "-"),   # planted NONCODING on the minus strand
    "novC": (90000, 93600, "+"),   # planted CODING-like: long ORF (must be filtered)
}


def write_reference(ref):
    rng = random.Random(20261012)
    # genome: 300 kb chr1
    seq = "".join(rng.choice("ACGT") for _ in range(300000))
    with open(ref / "genome.fa", "w") as fa:
        fa.write(">chr1\n")
        for i in range(0, len(seq), 60):
            fa.write(seq[i:i + 60] + "\n")
    # .fai placeholder (correct shape; dry-run checks existence only)
    with open(ref / "genome.fa.fai", "w") as fai:
        fai.write("chr1\t300000\t6\t60\t61\n")

    # reference annotation: 5 genes with 2 exons each, every 60 kb
    with open(ref / "genes.gtf", "w") as gtf, open(ref / "genes_pcg.gtf", "w") as pcg:
        gtf.write("#!genome-build test\n")
        pcg.write("#!genome-build test (protein-coding CDS only)\n")
        for k in range(5):
            start = 2000 + k * 60000
            gene = f"gene{k + 1:02d}"
            tx = f"tx{k + 1:02d}"
            for feat, s, e in (("exon", start, start + 900), ("exon", start + 3000, start + 3900)):
                line = (f"chr1\ttest\t{feat}\t{s}\t{e}\t.\t+\t.\t"
                        f"gene_id \"{gene}\"; transcript_id \"{tx}\";\n")
                gtf.write(line)
                pcg.write(line)
            gtf.write(f"chr1\ttest\ttranscript\t{start}\t{start + 3900}\t.\t+\t.\t"
                      f"gene_id \"{gene}\"; transcript_id \"{tx}\";\n")
            pcg.write(f"chr1\ttest\ttranscript\t{start}\t{start + 3900}\t.\t+\t.\t"
                      f"gene_id \"{gene}\"; transcript_id \"{tx}\";\n")

    # CPAT placeholders (the real model/hexamer are species params; the
    # fixture run only needs existing files)
    (ref / "cpat_model.logit").write_text(
        "# placeholder CPAT logit model (fixture)\n", encoding="utf-8")
    (ref / "hexamer.tsv").write_text(
        "# placeholder hexamer table (fixture)\nAAAAA\t0.1\n", encoding="utf-8")

    # mini Pfam DB: a harmless small HMM header file (real run deferred)
    (ref / "pfam_db").mkdir(parents=True, exist_ok=True)
    (ref / "pfam_db" / "Pfam-A.hmm").write_text(
        "HMMER3/f [i1 placeholder]\nNAME fixture\nACC PF00000.1\nDESC fixture\n"
        "LENG 40\nCKSUM 0\nGA 25 25\nTC 25 25\nNC 22 22\n//\n", encoding="utf-8")

    # mini DIAMOND DB placeholder (real run deferred: built with diamond makedb)
    (ref / "nr.dmnd").write_bytes(b"")


def write_novel_gtfs(ref):
    """Write the two fixture-side novel-transcript artifacts the dry-run
    consumes downstream (the gffcompare/cpc2/cpat steps produce real ones at
    execution; for the dry-run only file names resolve)."""
    rows = []
    for tid, (start, end, strand) in sorted(NOVEL_LOCI.items()):
        mid = (start + end) // 2
        rows.append((start, mid))
        rows.append((mid + 1, end))
    with open(ref / "novel.expected.gtf", "w") as gtf:
        for tid, (start, end, strand) in sorted(NOVEL_LOCI.items()):
            mid = (start + end) // 2
            gtf.write(f"chr1\ttest\texon\t{start}\t{mid}\t.\t{strand}\t.\t"
                      f"gene_id \"{tid}\"; transcript_id \"{tid}\";\n")
            gtf.write(f"chr1\ttest\texon\t{mid + 1}\t{end}\t.\t{strand}\t.\t"
                      f"gene_id \"{tid}\"; transcript_id \"{tid}\";\n")
            gtf.write(f"chr1\ttest\ttranscript\t{start}\t{end}\t.\t{strand}\t.\t"
                      f"gene_id \"{tid}\"; transcript_id \"{tid}\";\n")


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

    for s in SAMPLES:
        (bam_dir / f"{s}.markdup.sorted.bam").touch()
        (bam_dir / f"{s}.markdup.sorted.bam.bai").touch()

    write_reference(ref)
    write_novel_gtfs(ref)
    print(f"fixtures written under {FX}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
