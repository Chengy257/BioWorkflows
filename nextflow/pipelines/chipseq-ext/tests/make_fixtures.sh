#!/usr/bin/env bash
# Create the synthetic fixture tree consumed by the `test` / `test_full`
# profiles (dry-run DAG baselines; file contents are minimal but parseable,
# so the summary scripts can also run for real outside containers).
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fx="${here}/fixtures"

rm -rf "${fx}"
mkdir -p \
    "${fx}/results/bwa/filtered_bam" \
    "${fx}/results/bwa/merged_library/macs3/narrow_peak/consensus/H3K27ac" \
    "${fx}/results/bwa/bigwig" \
    "${fx}/results/bwa/flagstat" \
    "${fx}/results/bwa/duplicateMarked"

cat > "${fx}/samplesheet.csv" <<'EOF'
sample_id,group,condition,control,batch
s1,H3K27ac_WT,WT,ctl1,
s2,H3K27ac_WT,WT,ctl1,
s3,H3K27ac_MUT,mut,ctl1,
s4,H3K27ac_MUT,mut,ctl1,
ctl1,INPUT,input,,
EOF

# Minimal flagstat content (parseable by qc_gates_summary.py)
for s in s1 s2 s3 s4 ctl1; do
    cat > "${fx}/results/bwa/flagstat/${s}.mLb.clN.sorted.flagstat" <<EOF
1000 + 0 in total (QC-passed reads + QC-failed reads)
900 + 0 mapped (90.00% : N/A)
EOF
done

# Minimal picard MarkDuplicates metrics (PERCENT_DUPLICATION = column 9)
for s in s1 s2 s3 s4 ctl1; do
    cat > "${fx}/results/bwa/duplicateMarked/${s}.duplicate_metrics" <<EOF
## METRICS CLASS	picard.sam.DuplicationMetrics
LIBRARY	UNPAIRED_READS_EXAMINED	READ_PAIRS_EXAMINED	SECONDARY_OR_SUPPLEMENTARY_RDS	UNMAPPED_READS	UNPAIRED_READ_DUPLICATES	READ_PAIR_DUPLICATES	READ_PAIR_OPTICAL_DUPLICATES	PERCENT_DUPLICATION
${s}	100	400	0	100	20	80	0	0.2000
EOF
done

# BAM/BAI/bigWig placeholders (empty; dry-run only resolves names)
for s in s1 s2 s3 s4 ctl1; do
    : > "${fx}/results/bwa/filtered_bam/${s}.mLb.clN.sorted.bam"
    : > "${fx}/results/bwa/filtered_bam/${s}.mLb.clN.sorted.bam.bai"
    : > "${fx}/results/bwa/bigwig/${s}.mLb.clN.bigWig"
done

# Per-sample MACS3 narrowPeaks: 25 real 10-column rows per sample (the idr
# tool requires >= 20 peaks post-merge), coordinates jittered per sample
# around shared peak loci so the pairwise IDR stage sees overlap signal.
# Loci follow L(k) = 3000 + k * 11500 on a 300 kb chr1 (MUST match
# tests/make_real_bams.sh): the spacing keeps DiffBind summit recentering
# (summit_flank 250 -> 500 bp windows) from chain-merging adjacent windows,
# which on the old 10 kb toy chromosome collapsed every consensus region
# into one interval (a DiffBind degenerate case). Column 10 carries a
# realistic interior summit offset (200), not the peak width.
python3 - "$fx" <<'PYEOF'
import random, sys
fx = sys.argv[1]
for si, s in enumerate(("s1", "s2", "s3", "s4")):
    rng = random.Random(100 + si)
    rows = []
    for k in range(25):
        start = 3000 + k * 11500 + rng.randint(-40, 40)
        score = rng.randint(60, 600)
        rows.append(f"chr1\t{start}\t{start + 400}\t{s}_{k+1}\t{score}\t.\t5.5\t25\t25\t200")
    with open(f"{fx}/results/bwa/merged_library/macs3/narrow_peak/{s}_peaks.narrowPeak", "w") as fh:
        fh.write("\n".join(rows) + "\n")
PYEOF

# Consensus peak BED: 25 regions mirroring the per-sample peak loci
python3 - "$fx" <<'PYEOF'
import sys
fx = sys.argv[1]
rows = [f"chr1\t{3000 + k * 11500}\t{3000 + k * 11500 + 200}\tH3K27ac_{k+1}\t0\t+" for k in range(25)]
with open(f"{fx}/results/bwa/merged_library/macs3/narrow_peak/consensus/H3K27ac/H3K27ac.consensus_peaks.bed", "w") as fh:
    fh.write("\n".join(rows) + "\n")
PYEOF

echo "fixtures ready under ${fx}"
