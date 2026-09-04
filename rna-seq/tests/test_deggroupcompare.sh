#!/bin/bash
#########################################################################
# Regression test for DEGgroupCompare.sh argument propagation (no R/snakemake needed)
#   Uses an argument-recording Rscript shim to verify that species is propagated to every
#   run_deg_compare.R job; also verifies no temporary files are left in the run directory.
# Usage: bash tests/test_deggroupcompare.sh
# Requires: bash, python3 (or python), mktemp; results are written to a temp directory and cleaned up automatically.
#########################################################################
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TARGET="$REPO_DIR/workflow/scripts/DEGgroupCompare.sh"

T="$(mktemp -d)"
trap 'rm -rf "$T"' EXIT
mkdir -p "$T/cwd" "$T/bin" "$T/results/5.DEG/DEGs" "$T/results/6.DEGcompare"
: > "$T/args.log"

## ---- Test fixtures: sample table (2 treatment groups + control), annotation table, two intersection files ----
## The genes of intersection 1 exist in the annotation table (fgrep success path); the genes of
## intersection 2 do not (fgrep failure path; the old implementation left a tmp.gene temp file behind)
cat > "$T/samples.csv" <<'EOF'
id,group
s1,treatA
s2,treatB
s3,control
EOF
printf 'geneX\tannotation of geneX\ngeneY\tannotation of geneY\n' > "$T/anno.tsv"
printf 'treatA\ttreatB\t2\ttreatA\342\210\251treatB\tgeneX geneY\n' \
    > "$T/results/6.DEGcompare/treatA_treatB_intersetion.xls"
printf 'treatA\ttreatC\t2\ttreatA\342\210\251treatC\tgeneZ\n' \
    > "$T/results/6.DEGcompare/treatA_treatC_intersetion.xls"

## ---- Shims: python works around local site-packages issues; Rscript records arguments line by line ----
PY=python3
command -v python3 >/dev/null 2>&1 || PY=python
printf '#!/bin/bash\nexec %s -S "$@"\n' "$PY" > "$T/bin/python"
printf '#!/bin/bash\nprintf "RLIB:%%s\\n" "${R_LIBS_USER:-}" >> %q\nprintf "ARG:%%s\\n" "$@" >> %q\n' "$T/args.log" "$T/args.log" > "$T/bin/Rscript"
chmod +x "$T/bin/python" "$T/bin/Rscript"

fails=0
check() {  ## check <description> <command...>
    local desc=$1; shift
    if "$@"; then echo "  PASS  $desc"; else echo "  FAIL  $desc"; fails=$((fails + 1)); fi
}

run_case() {  ## run_case <species>
    local spe=$1
    : > "$T/args.log"
    rm -f "$T/cwd/tmp.gene"
    (cd "$T/cwd" && PATH="$T/bin:$PATH" RNASEQ_RSCRIPT="$T/bin/Rscript" RNASEQ_PYTHON="$T/bin/python" R_LIBS_USER="$T/rlibA:$T/rlibB" "$TARGET" \
        "$T/results" control "$T/anno.tsv" "$T/samples.csv" 2 "$spe" \
        >/dev/null 2>&1)
}

echo "[test] case 1: species=hsa propagated to all jobs"
run_case hsa
## 2 treatment groups -> 1 pair x UP/DOWN/ALL 3x3 = 9 jobs
check "all 9 jobs received species=hsa" \
    test "$(grep -c '^ARG:hsa$' "$T/args.log")" = 9
check "all 9 R sub-jobs inherit R_LIBS_USER" \
    test "$(grep -c "^RLIB:$T/rlibA:$T/rlibB$" "$T/args.log")" = 9
check "flag.log produced" test -f "$T/results/6.DEGcompare/flag.log"
check "intersection gene annotation produced (fgrep success path)" \
    test -f "$T/results/6.DEGcompare/intersetionGene_anno/treatA_treatB_intersetionGene_anno.xls"
check "annotation content contains the intersection genes" \
    grep -q geneX "$T/results/6.DEGcompare/intersetionGene_anno/treatA_treatB_intersetionGene_anno.xls"
check "fgrep failure path leaves no tmp.gene behind" test ! -e "$T/cwd/tmp.gene"

echo "[test] case 2: species=osa argument propagation"
run_case osa
check "all 9 jobs received species=osa" \
    test "$(grep -c '^ARG:osa$' "$T/args.log")" = 9
check "flag.log produced" test -f "$T/results/6.DEGcompare/flag.log"

if [[ "$fails" -gt 0 ]]; then
    echo "[test] $fails assertions failed"
    exit 1
fi
echo "[test] all checks passed (DEGgroupCompare argument propagation + temp file hygiene)"
