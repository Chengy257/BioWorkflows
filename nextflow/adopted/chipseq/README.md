# adopted/chipseq - nf-core/chipseq

Thin adoption entry for [nf-core/chipseq](https://nf-co.re/chipseq), pinned to
dev HEAD `1b3078d71c4c8126a9b976d3b0248e9b9b525608` (the 2.1.0 release is
incompatible with Nextflow 26; see `nextflow/versions.yml`). This is **not a
fork**: no pipeline code lives here. The entry carries only the convention
layer permitted by the design doc (Section 8): a params file, reference
configuration example, an HPC example command, and pointers to the scientific
decision and validation records.

## Contents

| File | Purpose |
|---|---|
| `params.example.yaml` | Example params file (rice shown as the example; every path is a placeholder - references are always user-configured, design doc D4) |
| `run-chipseq.sh` | Example HPC launcher wrapping `nextflow run nf-core/chipseq` with the pinned revision |

## Usage

```bash
# 1. Prepare a project params file from the example (fill in real reference
#    paths and the sample sheet); keep it in your project directory, not here.
cp params.example.yaml my-project.params.yaml

# 2. Launch (PBS via the machine-local executor setting; apptainer profile
#    comes from nextflow/conf/hpc.config):
bash run-chipseq.sh my-project.params.yaml
```

Environment variables honored by `run-chipseq.sh` (all optional, with
defaults): `CHIPSEQ_REVISION`, `CHIPSEQ_PROFILE`, `CHIPSEQ_OUTDIR`,
`CHIPSEQ_CACHE`, `CHIPSEQ_EXTRA_ARGS`.

## Local extensions

Capabilities that the upstream pipeline does not provide (DiffBind, spike-in,
IDR, SEACR caller, HOMER motif discovery, QC gate table, organelle QC) are
delivered by `nextflow/pipelines/chipseq-ext/` as default-off local stages
consuming this pipeline's published results. See that pipeline's README and
the decision record.

## Records

- Scientific decision record:
  `nextflow/docs/chipseq-scientific-decision-record.md` (ratified 2026-10-06).
- Validation notes (smoke run, container seeding, extension dry-run
  baselines): `nextflow/docs/chipseq-validation-notes.md`.
- Image acquisition recipe (per-source routing): implementation guidance
  Section 3.1 (`plans/2026-10-05-nextflow-migration-implementation-guidance.md`).
