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

# Per-sample MACS3 narrowPeaks: two real 10-column rows per sample
for s in s1 s2 s3 s4; do
    cat > "${fx}/results/bwa/merged_library/macs3/narrow_peak/${s}_peaks.narrowPeak" <<EOF
chr1	1000	1200	${s}_1	250	.	5.5	25	25	100
chr1	5000	5300	${s}_2	180	.	4.2	18	18	150
EOF
done

# Consensus peak BED for one antibody (native MACS3_CONSENSUS output shape)
cat > "${fx}/results/bwa/merged_library/macs3/narrow_peak/consensus/H3K27ac/H3K27ac.consensus_peaks.bed" <<'EOF'
chr1	1000	1200	H3K27ac_1	0	+
chr1	5000	5300	H3K27ac_2	0	+
EOF

echo "fixtures ready under ${fx}"
