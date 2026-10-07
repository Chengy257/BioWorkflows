#!/usr/bin/env bash
# Upgrade the fixture BAMs from empty placeholders to real tiny BAMs so the
# container-based extension stages can execute end-to-end on this machine.
#
# Usage: SAMTOOLS_BIN="<samtools command>" bash tests/make_real_bams.sh
#   SAMTOOLS_BIN may be a plain samtools on PATH or a container invocation,
#   e.g.: SAMTOOLS_BIN="apptainer exec /path/to/samtools.sif samtools"
#
# Each sample gets mapped reads on chr1 spread over the 25 fixture peak
# loci, 10 on spike1 (the spike-in contig) and 5 unmapped reads. Also writes
# tests/fixtures/genome/chr1.fa (the custom genome used by the HOMER stage,
# which cannot fetch prebuilt genomes offline).
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fx="${here}/fixtures"
samtools_bin="${SAMTOOLS_BIN:?usage: SAMTOOLS_BIN=<samtools command> bash make_real_bams.sh}"

mkdir -p "${fx}/genome" "${fx}/results/bwa/filtered_bam"

# Shared reference sequences: 300 kb chr1 (MUST match the peak loci in
# tests/make_fixtures.sh), 1 kb spike1 (spike-in)
python3 - "$fx" <<'PYEOF'
import random, sys
fx = sys.argv[1]
rng = random.Random(42)
for name, length in (("chr1", 300000), ("spike1", 1000)):
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
import math, os, random
seed, samfile = int(os.environ["seed"]), os.environ["samfile"]
rng = random.Random(seed * 7919)
read_seq = "ACGTTGCAAGGCTTACGATCGATCGATCGATCGATCGATCGATCGATCGA"
qual = "F" * 50
# Peak loci MUST match tests/make_fixtures.sh: L(k) = 3000 + k * 11500.
loci = [3000 + k * 11500 for k in range(25)]
reads = []


def add_pair(name, pos):
    reads.append(f"{name}\t99\tchr1\t{pos}\t60\t50M\t=\t{pos+120}\t{pos+170}\t{read_seq}\t{qual}")
    reads.append(f"{name}\t147\tchr1\t{pos+120}\t60\t50M\t=\t{pos}\t-{pos+170}\t{read_seq}\t{qual}")


def add_single(name, pos, ref="chr1"):
    reads.append(f"{name}\t0\t{ref}\t{pos}\t60\t50M\t*\t0\t0\t{read_seq}\t{qual}")


if seed == 5:
    # ctl1 is the shared input control: uniform background over the whole
    # chromosome plus three sharp clusters at treated loci (k = 4, 14, 23).
    # The clusters give the SEACR normalized path a control feature set to
    # derive empirical FDR thresholds from, while most treated loci stay
    # control-sparse and pass.
    for i in range(300):
        add_single(f"bg{i}", rng.randint(1, 299000))
    for c, k in enumerate((4, 14, 23)):
        for i in range(40):
            add_single(f"cl{c}_{i}", loci[k] + rng.randint(-200, 200))
else:
    # s1/s2 = WT, s3/s4 = MUT. Locus baselines are drawn lognormally so
    # per-region counts span an order of magnitude with multiplicative
    # sample noise - Poisson-flat toy counts collapse DESeq2's dispersion
    # fit ("all gene-wise dispersion estimates within 2 orders of
    # magnitude"). Loci 21-23 are WT-only, 24-25 MUT-only, so the
    # differential-binding contrast is non-degenerate in both directions.
    def pairs(k):
        base = math.exp(rng.uniform(math.log(4), math.log(60)))
        if k < 20:
            n = base * math.exp(rng.gauss(0.0, 0.6))
        elif k < 23:
            n = base if seed <= 2 else rng.uniform(0.0, 1.5)
        else:
            n = rng.uniform(0.0, 1.5) if seed <= 2 else base
        return max(0, int(round(n)))

    for k, loc in enumerate(loci):
        for i in range(pairs(k)):
            add_pair(f"k{k}_r{i}", loc + rng.randint(-50, 50))
    for i in range(5):
        add_single(f"b{i}", rng.randint(30000, 270000))
for i in range(10):
    add_single(f"sp{i}", rng.randint(100, 500), ref="spike1")
for i in range(5):
    reads.append(f"u{i}\t4\t*\t0\t0\t*\t*\t0\t0\t{read_seq}\t{qual}")
reads.sort(key=lambda r: (r.split("\t")[2], int(r.split("\t")[3])))
with open(samfile, "w") as fh:
    fh.write("@HD\tVN:1.6\tSO:coordinate\n")
    fh.write("@SQ\tSN:chr1\tLN:300000\n")
    fh.write("@SQ\tSN:spike1\tLN:1000\n")
    fh.write("\n".join(reads) + "\n")
PYEOF
    $samtools_bin view -b -o "${sam}.tmp.bam" "$samtmp"
    $samtools_bin sort -o "$sam" "${sam}.tmp.bam"
    rm -f "${sam}.tmp.bam" "$samtmp"
    $samtools_bin index "$sam"
    echo "  real BAM: $s"
done

echo "real fixture BAMs + genome FASTA ready under ${fx}"
