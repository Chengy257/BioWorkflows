#!/bin/bash
#########################################################################
# Between-group DEG set comparison: generate pairwise combination jobs -> parallel GO enrichment comparison -> intersecting-gene annotation
# Usage: DEGgroupCompare.sh <results_dir> <control_group> <annotation_tsv> <sample_info.csv> <threads> [species]
#   species        osa | hsa (default osa), forwarded to run_deg_compare.R
# Parallel execution uses xargs -P instead of the former ParaFly (no personal software paths required)
#########################################################################

set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RSCRIPT="${RNASEQ_RSCRIPT:-Rscript}"
PYTHON="${RNASEQ_PYTHON:-python3}"

RD=$1
control=$2
anno=$3
samplefile=$4
cpu=$5
spe=${6:-osa}

## Temporary file for intersection annotation (mktemp + trap avoids leftovers in the run directory)
tmp_gene=$(mktemp)
trap 'rm -f "$tmp_gene"' EXIT

mkdir -p "$RD/6.DEGcompare/combination"
"$PYTHON" "$SCRIPT_DIR/getGroups.py" "$samplefile" "$RD/5.DEG/DEGs" "$control" "$RD/6.DEGcompare"

## Generate the job list (species / orgdb propagated to every job)
rm -f "$RD/6.DEGcompare/jobs"
ls "$RD"/6.DEGcompare/combination/* | while read line;
do
    printf '%q %q %q %q\n' "$RSCRIPT" "$SCRIPT_DIR/run_deg_compare.R" "$line" "$spe" >> "$RD/6.DEGcompare/jobs"
done

## Run in parallel (one command per line, -P controls concurrency; -r handles an empty job list when <2 treatment groups)
## A single failing job does not abort the whole run (same fault-tolerance semantics as enrich), but a warning is emitted
xargs -r -P "$cpu" -d '\n' -I CMD sh -c 'CMD' < "$RD/6.DEGcompare/jobs" \
    || echo "[DEGgroupCompare] [WARN] Some comparison jobs failed; results may be incomplete." >&2

## Functional annotation of intersecting genes
mkdir -p "$RD/6.DEGcompare/intersetionGene_anno"
ls "$RD"/6.DEGcompare/*_intersetion.xls | while read id;
do
    prefix=$(basename "$id")
    cat "$id" | fgrep ∩ | cut -f5 | tr " " "\n" | sed 's/,//g' > "$tmp_gene"
    cat "$anno" | fgrep -w -f "$tmp_gene" > "$RD/6.DEGcompare/intersetionGene_anno/${prefix%_*}_intersetionGene_anno.xls"
done

echo "$(date) : All Done!" > "$RD/6.DEGcompare/flag.log"
