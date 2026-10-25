#!/usr/bin/env bash
# Build the mini DIAMOND / Pfam databases and the CPAT model for the
# lncrna-nf real-container fixture round (design: the fixtures carry raw
# material - toy_proteins.fa / pfam_seed.sto / cpat_{coding,noncoding}.fa -
# and this script turns them into what the pipeline consumes).
#
# Everything runs inside the same containers the pipeline uses (images from
# the apptainer cache; pull routes recorded in the validation notes: depot
# direct, quay via the HTTP mixed-port proxy). Host needs: bash, apptainer,
# python3.
#
# Usage: bash tests/build_fixture_dbs.sh
set -euo pipefail

HERE=$(cd "$(dirname "${0}")" && pwd)
REF="$HERE/fixtures/reference"
CACHE="${APPTAINER_CACHE:-$HOME/soft/apptainer-cache}"

if [ -z "${APPTAINER_BIN:-}" ]; then
    for cand in apptainer "$HOME/soft/miniconda3/envs/apptainer/bin/apptainer" singularity; do
        if command -v "$cand" >/dev/null 2>&1; then APPTAINER_BIN="$cand"; break; fi
    done
fi
[ -n "${APPTAINER_BIN:-}" ] || { echo "apptainer/singularity not found" >&2; exit 1; }

img() { # img <cache-glob> - first matching cached image
    local found
    found=$(ls $CACHE/$1 2>/dev/null | head -1 || true)
    [ -n "$found" ] || { echo "image not in cache: $1" >&2; exit 1; }
    echo "$found"
}

DIAMOND_IMG="${DIAMOND_IMG:-$(img 'depot.galaxyproject.org-singularity-diamond-*')}"
PFAM_IMG="${PFAM_IMG:-$(img '*pfam_scan*')}"
CPAT_IMG="${CPAT_IMG:-$(img '*cpat*')}"

echo "== DIAMOND makedb ($DIAMOND_IMG)"
rm -f "$REF/nr.dmnd"
"$APPTAINER_BIN" exec --bind "$REF":/work "$DIAMOND_IMG" \
    bash -c 'cd /work && diamond makedb --in toy_proteins.fa -d nr && ls -la nr.dmnd'

echo "== Pfam hmmbuild + hmmpress ($PFAM_IMG)"
rm -rf "$REF/pfam_db" && mkdir -p "$REF/pfam_db"
"$APPTAINER_BIN" exec --bind "$REF":/work "$PFAM_IMG" \
    bash -c 'cd /work && hmmbuild --amino pfam_db/Pfam-A.hmm pfam_seed.sto \
        && hmmpress pfam_db/Pfam-A.hmm \
        && ls -la pfam_db/'
# pfam_scan.pl also requires non-empty Pfam-A.hmm.dat and active_site.dat
# (a minimal one-family entry matches the toy HMM)
printf 'PF99999.1\tPF99999\t1\tpfam_seed\tFixture synthetic domain\tPfam-A\t1\tfixture\t0\t0\t0\t0\t0\t0\t0\t0\n' \
    > "$REF/pfam_db/Pfam-A.hmm.dat"
printf 'PF99999\t1\n' > "$REF/pfam_db/active_site.dat"

echo "== CPAT hexamer + model ($CPAT_IMG)"
# CPAT 3.0.5 ships the builders as separate scripts: make_hexamer_tab.py
# writes the table to stdout, make_logitModel.py consumes it via -x
"$APPTAINER_BIN" exec --bind "$REF":/work "$CPAT_IMG" \
    bash -c 'cd /work &&
        make_hexamer_tab.py -c cpat_coding.fa -n cpat_noncoding.fa > hexamer.tsv &&
        make_logitModel.py -c cpat_coding.fa -n cpat_noncoding.fa -x hexamer.tsv -o cpat_model &&
        ls -la hexamer.tsv cpat_model*'
if [ ! -f "$REF/cpat_model.logit" ]; then
    found=$(ls "$REF"/cpat_model* 2>/dev/null | grep -i rdata | head -1 || true)
    [ -n "$found" ] && cp "$found" "$REF/cpat_model.logit"
fi
[ -f "$REF/cpat_model.logit" ] || { echo "CPAT model build failed (see output above)" >&2; exit 1; }
[ -s "$REF/hexamer.tsv" ] || { echo "hexamer table empty" >&2; exit 1; }

echo "fixture DBs ready under $REF (nr.dmnd, pfam_db/, cpat_model.logit, hexamer.tsv)"
