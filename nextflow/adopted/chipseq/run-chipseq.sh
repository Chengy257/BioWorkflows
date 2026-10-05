#!/usr/bin/env bash
# Example HPC launcher for the adopted nf-core/chipseq entry.
#
# Container-only on the HPC (design doc Section 15 amendment 5). The Nextflow
# executor and container cache directory come from machine-local, git-ignored
# configuration (~/.nextflow/config); this wrapper only pins the pipeline
# revision and applies the tracked shared config.
#
# Usage: bash run-chipseq.sh <project.params.yaml> [extra nextflow args]
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

revision="${CHIPSEQ_REVISION:-1b3078d71c4c8126a9b976d3b0248e9b9b525608}"
profile="${CHIPSEQ_PROFILE:-apptainer}"
outdir="${CHIPSEQ_OUTDIR:-results}"
cache="${CHIPSEQ_CACHE:-}"

params_file="${1:?usage: run-chipseq.sh <project.params.yaml> [extra args]}"
shift || true

extra=()
if [[ -n "${cache}" ]]; then
    extra+=(-w "${cache}")
fi

cd "$(dirname "${params_file}")"

# -ansi-log false keeps logs readable inside nohup/qsub output files.
exec nextflow run nf-core/chipseq \
    -revision "${revision}" \
    -profile "${profile}" \
    -c "${here}/../../conf/hpc.config" \
    -params-file "${params_file}" \
    --outdir "${outdir}" \
    -resume \
    -ansi-log false \
    "${extra[@]}" "$@"
