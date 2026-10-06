#!/usr/bin/env bash
# Upgrade the fixture BAMs from empty placeholders to real tiny BAMs so the
# container-based extension stages can execute end-to-end on this machine.
#
# Usage: SAMTOOLS_BIN="<samtools command>" bash tests/make_real_bams.sh
#   SAMTOOLS_BIN may be a plain samtools on PATH or a container invocation,
#   e.g.: SAMTOOLS_BIN="apptainer exec /path/to/samtools.sif samtools"
#
# Each sample gets ~40 mapped reads on chr1, 10 on chrS (the spike-in contig)
# and 5 unmapped reads. Also writes tests/fixtures/genome/chr1.fa (the custom
# genome used by the HOMER stage, which cannot fetch prebuilt genomes
# offline).
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fx="${here}/fixtures"
samtools_bin="${SAMTOOLS_BIN:?usage: SAMTOOLS_BIN=<samtools command> bash make_real_bams.sh}"

mkdir -p "${fx}/genome" "${fx}/results/bwa/filtered_bam"

# Shared reference sequences: 10 kb chr1, 1 kb chrS (spike-in)
python3 - "$fx" <<'PYEOF'
import random, sys
fx = sys.argv[1]
rng = random.Random(42)
for name, length in (("chr1", 10000), ("chrS", 1000)):
    seq = ''.join(rng.choice('ACGT') for _ in range(length))
    with open(f"{fx}/genome/{name}.fa", "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")
print("reference sequences written")
PYEOF

seed=0
for s in s1 s2 s3 s4 ctl1; do
    seed=$((seed + 1))
    sam="${fx}/results/bwa/filtered_bam/${s}.mLb.clN.sorted.bam"
    samtmp="${fx}/results/bwa/filtered_bam/.${s}.sam"
    fx="$fx" seed="$seed" samfile="$samtmp" python3 - <<'PYEOF'
import os, random
seed, samfile = int(os.environ["seed"]), os.environ["samfile"]
rng = random.Random(seed * 7919)
read_seq = "ACGTTGCAAGGCTTACGATCGATCGATCGATCGATCGATCGATCGATCGA"
qual = "F" * 50
reads = []
if seed == 5:
    # ctl1 is the shared input control: uniform background plus three sharp
    # clusters so the SEACR empirical-FDR path has real peak/background
    # contrast (it degenerates on flat toy coverage)
    for i in range(40):
        pos = rng.randint(1, 9900)
        reads.append(f"bg{i}\t0\tchr1\t{pos}\t60\t50M\t*\t0\t0\t{read_seq}\t{qual}")
    for c, center in enumerate((1500, 4000, 7000)):
        for i in range(40):
            pos = center + rng.randint(-200, 200)
            reads.append(f"cl{c}_{i}\t0\tchr1\t{pos}\t60\t50M\t*\t0\t0\t{read_seq}\t{qual}")
    reads.sort(key=lambda r: (r.split("\t")[2], int(r.split("\t")[3])))
    with open(samfile, "w") as fh:
        fh.write("@HD\tVN:1.6\tSO:coordinate\n")
        fh.write("@SQ\tSN:chr1\tLN:10000\n")
        fh.write("@SQ\tSN:chrS\tLN:1000\n")
        fh.write("\n".join(reads) + "\n")
    raise SystemExit(0)
# reads biased into the fixture peak regions (chr1:1000-1200, 5000-5300) so
# DiffBind counts are non-zero; plus spike-in chrS and unmapped reads
for i in range(25):
    pos = rng.randint(1000, 1150)
    reads.append(f"r{i}\t99\tchr1\t{pos}\t60\t50M\t=\t{pos+120}\t{pos+170}\t{read_seq}\t{qual}")
    reads.append(f"r{i}\t147\tchr1\t{pos+120}\t60\t50M\t=\t{pos}\t-{pos+170}\t{read_seq}\t{qual}")
for i in range(10):
    pos = rng.randint(5000, 5250)
    reads.append(f"p{i}\t99\tchr1\t{pos}\t60\t50M\t=\t{pos+120}\t{pos+170}\t{read_seq}\t{qual}")
    reads.append(f"p{i}\t147\tchr1\t{pos+120}\t60\t50M\t=\t{pos}\t-{pos+170}\t{read_seq}\t{qual}")
for i in range(5):
    pos = rng.randint(2000, 4000)
    reads.append(f"b{i}\t0\tchr1\t{pos}\t60\t50M\t*\t0\t0\t{read_seq}\t{qual}")
for i in range(10):
    pos = rng.randint(100, 500)
    reads.append(f"sp{i}\t0\tchrS\t{pos}\t60\t50M\t*\t0\t0\t{read_seq}\t{qual}")
for i in range(5):
    reads.append(f"u{i}\t4\t*\t0\t0\t*\t*\t0\t0\t{read_seq}\t{qual}")
reads.sort(key=lambda r: (r.split("\t")[2], int(r.split("\t")[3])))
with open(samfile, "w") as fh:
    fh.write("@HD\tVN:1.6\tSO:coordinate\n")
    fh.write("@SQ\tSN:chr1\tLN:10000\n")
    fh.write("@SQ\tSN:chrS\tLN:1000\n")
    fh.write("\n".join(reads) + "\n")
PYEOF
    $samtools_bin view -b -o "${sam}.tmp.bam" "$samtmp"
    $samtools_bin sort -o "$sam" "${sam}.tmp.bam"
    rm -f "${sam}.tmp.bam" "$samtmp"
    $samtools_bin index "$sam"
    echo "  real BAM: $s"
done

echo "real fixture BAMs + genome FASTA ready under ${fx}"
