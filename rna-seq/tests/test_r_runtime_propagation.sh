#!/bin/bash
# Verify nested enrichment shell calls inherit the configured Rscript and R library path.
set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TARGET="$REPO_DIR/workflow/scripts/enrich.sh"
T="$(mktemp -d)"
trap 'rm -rf "$T"' EXIT
mkdir -p "$T/deg/Diff_Expr_Analysis_Results" "$T/bin"
cat > "$T/deg/Diff_Expr_Analysis_Results/treat_vs_control_DESeq2.output.tsv" <<'TSV'
gene_id	baseMean	log2FoldChange	type
geneA	10	2	Up
geneB	10	-2	Down
TSV
printf 'geneA\tannoA\ngeneB\tannoB\n' > "$T/anno.tsv"
: > "$T/r.log"
cat > "$T/bin/custom-Rscript" <<EOS
#!/bin/bash
printf 'RLIB:%s\n' "\${R_LIBS_USER:-}" >> "$T/r.log"
printf 'SCRIPT:%s\n' "\${1:-}" >> "$T/r.log"
exit 0
EOS
chmod +x "$T/bin/custom-Rscript"

RNASEQ_RSCRIPT="$T/bin/custom-Rscript" R_LIBS_USER="$T/libA:$T/libB" \
    "$TARGET" "$T/deg" hsa "$T/anno.tsv" hsa >/dev/null

test "$(grep -c "^RLIB:$T/libA:$T/libB$" "$T/r.log")" = 4
test "$(grep -c '^SCRIPT:.*run_enrichment.R$' "$T/r.log")" = 3
test "$(grep -c '^SCRIPT:.*run_gsea.R$' "$T/r.log")" = 1
test -f "$T/deg/GO_KEGG_enrich/flag.log"
echo "[test] nested Rscript/R_LIBS propagation PASS"
