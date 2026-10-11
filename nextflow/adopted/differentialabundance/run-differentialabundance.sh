#!/usr/bin/env bash
# Example HPC launcher for the adopted nf-core/differentialabundance entry.
#
# Container-only on the HPC (design doc Section 15 amendment 5). The Nextflow
# executor and container cache directory come from machine-local, git-ignored
# configuration (~/.nextflow/config); this wrapper only pins the pipeline
# revision and applies the tracked shared config.
#
# Usage: bash run-differentialabundance.sh <project.params.yaml> [extra nextflow args]
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

revision="${DIFFAB_REVISION:-2.0.0}"
profile="${DIFFAB_PROFILE:-apptainer}"
outdir="${DIFFAB_OUTDIR:-results-da}"
cache="${DIFFAB_CACHE:-}"

params_file="${1:?usage: run-differentialabundance.sh <project.params.yaml> [extra args]}"
shift || true

extra=()
if [[ -n "${cache}" ]]; then
    extra+=(-w "${cache}")
fi

cd "$(dirname "${params_file}")"

# -ansi-log false keeps logs readable inside nohup/qsub output files.
exec nextflow run nf-core/differentialabundance \
    -revision "${revision}" \
    -profile "${profile}" \
    -c "${here}/../../conf/hpc.config" \
    -params-file "${params_file}" \
    --outdir "${outdir}" \
    -resume \
    -ansi-log false \
    "${extra[@]}" "$@"
