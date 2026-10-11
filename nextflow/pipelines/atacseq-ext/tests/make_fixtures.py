#!/usr/bin/env python3
"""Create the synthetic fixture tree consumed by the `test` / `test_full`
profiles (dry-run DAG baselines). Fixed seed; the layout mirrors the PINNED
nf-core/atacseq dev-HEAD publish contract (README "Consumed published
contract"):

    fixtures/
    ├── samplesheet.csv
    ├── reference/genome.fa(+.fai), motifs.jaspar
    └── results/bowtie2/merged_library/
        ├── {id}.mLb.clN.sorted.bam(+.bai)          (placeholder; dry-run)
        ├── samtools_stats/{id}.mLb.clN.sorted.bam.{stats,flagstat,idxstats}
        ├── picard_metrics/{id}.mLb.mkD.sorted.metrics.txt
        └── macs3/broad_peak/consensus/consensus_peaks.mLb.clN.bed

BAM/BAI files are empty placeholders: the dry-run baselines only resolve
names. Real-container fixture execution (deferred round) needs a real-BAM
helper mirroring chipseq-ext's make_real_bams.sh - folded into that round.
The idxstats text files carry organelle + spike contigs so the helper
semantics and the (deferred) container run have a documented example of the
expected format.

Usage: python3 tests/make_fixtures.py
"""
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FX = HERE / "fixtures"

SAMPLES = ["a1", "a2", "b1", "b2"]
GROUPS = {"a1": "ATAC_WT", "a2": "ATAC_WT", "b1": "ATAC_MUT", "b2": "ATAC_MUT"}
CONDITIONS = {"a1": "WT", "a2": "WT", "b1": "mut", "b2": "mut"}

# toy contigs: 200 kb chr1 + chloroplast + mitochondria + spike1
CONTIGS = [("chr1", 200000), ("chrc", 150000), ("chrm", 8000), ("spike1", 1000)]
PEAK_LOCI = [10000, 40000, 70000, 100000, 130000]  # on chr1


def main():
    rng = random.Random(20261012)
    if FX.exists():
        print(f"removing existing fixture tree {FX}")
        import shutil
        shutil.rmtree(FX)

    ml = FX / "results" / "bowtie2" / "merged_library"
    (ml / "samtools_stats").mkdir(parents=True)
    (ml / "picard_metrics").mkdir(parents=True)
    (ml / "macs3" / "broad_peak" / "consensus").mkdir(parents=True)
    ref = FX / "reference"
    ref.mkdir(parents=True)

    # sample sheet (no controls: two 2-replicate treat groups)
    lines = ["sample_id,group,condition,control,batch"]
    for s in SAMPLES:
        lines.append(f"{s},{GROUPS[s]},{CONDITIONS[s]},,")
    (FX / "samplesheet.csv").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # reference: toy FASTA + .fai (TOBIAS ATACorrect needs the indexed genome)
    fa_lines, fai_lines, offset = [], [], 0
    for name, length in CONTIGS:
        seq = "".join(rng.choice("ACGT") for _ in range(length))
        fa_lines.append(f">{name}")
        fa_lines += [seq[i:i + 60] for i in range(0, length, 60)]
        fai_lines.append(f"{name}\t{length}\t{offset + len(name) + 1}\t60\t61")
        offset += len(name) + len(seq) // 60 + (60 if length % 60 else 0) + (length % 60 > 0) * 0 + len(f">{name}\n") + length + 1
    (ref / "genome.fa").write_text("\n".join(fa_lines) + "\n", encoding="utf-8")
    # a real .fai is computed by `samtools faidx`; for the dry-run the exact
    # byte offsets do not matter (only existence is checked) - write a
    # correctly-shaped placeholder and note it.
    (ref / "genome.fa.fai").write_text("\n".join(fai_lines) + "\n", encoding="utf-8")

    # motifs: two small JASPAR PFMs (TOBIAS reads JASPAR format)
    (ref / "motifs.jaspar").write_text(
        ">toyA\tToy motif A\n"
        "A [  0.0  1.0  0.0  0.5 ]\n"
        "C [  0.0  0.0  1.0  0.5 ]\n"
        "G [  1.0  0.0  0.0  0.0 ]\n"
        "T [  0.0  0.0  0.0  0.0 ]\n"
        ">toyB\tToy motif B\n"
        "A [  0.5  0.0  0.0  1.0 ]\n"
        "C [  0.5  1.0  0.0  0.0 ]\n"
        "G [  0.0  0.0  1.0  0.0 ]\n"
        "T [  0.0  0.0  0.0  0.0 ]\n",
        encoding="utf-8",
    )

    # per-sample artifacts
    for s in SAMPLES:
        # BAM/BAI placeholders (empty; dry-run only resolves names)
        (ml / f"{s}.mLb.clN.sorted.bam").touch()
        (ml / f"{s}.mLb.clN.sorted.bam.bai").touch()

        # flagstat (parseable by qc_gates_summary.py): mapping rate per sample
        total, mapped = 10000, 9000 if s.startswith("a") else 8600
        (ml / "samtools_stats" / f"{s}.mLb.clN.sorted.bam.flagstat").write_text(
            f"{total} + 0 in total (QC-passed reads + QC-failed reads)\n"
            f"{mapped} + 0 mapped (90.00% : N/A)\n",
            encoding="utf-8",
        )
        (ml / "samtools_stats" / f"{s}.mLb.clN.sorted.bam.stats").write_text(
            "This file is not parsed by the ext stages.\n", encoding="utf-8")
        # idxstats: chr1 + organelle contigs + spike1 (spike reads a-b split
        # by group so a real container run would produce group-differentiated
        # scale factors)
        spike = 10 if s.startswith("a") else 40
        idx = (
            f"chr1\t200000\t{mapped - 100 - spike}\t0\n"
            f"chrc\t150000\t60\t0\n"
            f"chrm\t8000\t40\t0\n"
            f"spike1\t1000\t{spike}\t0\n"
            f"*\t0\t0\t{total - mapped}\n"
        )
        (ml / "samtools_stats" / f"{s}.mLb.clN.sorted.bam.idxstats").write_text(
            idx, encoding="utf-8")

        # picard MarkDuplicates metrics (PERCENT_DUPLICATION = column 9)
        (ml / "picard_metrics" / f"{s}.mLb.mkD.sorted.metrics.txt").write_text(
            "## METRICS CLASS\tpicard.sam.DuplicationMetrics\n"
            "LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\t"
            "SECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\t"
            "UNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\t"
            "READ_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\n"
            f"{s}\t100\t4000\t0\t1000\t20\t800\t0\t0.2000\n",
            encoding="utf-8",
        )

    # consensus peaks (broad_peak): one BED file published flat by MACS3_CONSENSUS
    bed = []
    for k, start in enumerate(PEAK_LOCI):
        bed.append(f"chr1\t{start}\t{start + 400}\tconsensus_{k + 1}\t100\t.")
    (ml / "macs3" / "broad_peak" / "consensus" / "consensus_peaks.mLb.clN.bed").write_text(
        "\n".join(bed) + "\n", encoding="utf-8")

    print(f"fixtures written under {FX}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
