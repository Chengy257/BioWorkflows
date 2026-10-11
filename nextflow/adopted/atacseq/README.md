# adopted/atacseq - nf-core/atacseq

Thin adoption entry for [nf-core/atacseq](https://nf-co.re/atacseq), pinned to
dev HEAD `6a9307709662bad3e7bcc10aa873f7e81fdd1646` (v2.2.0dev; the 2.1.2
release is incompatible with Nextflow 26 - pin gate 2026-10-12, see
`nextflow/docs/atacseq-scientific-decision-record.md`). This is **not a
fork**: no pipeline code lives here. The entry carries only the convention
layer permitted by the design doc (Section 8): a params file, reference
configuration example, an HPC example command, and pointers to the scientific
decision and validation records.

The pin is dev HEAD because the 2.1.2 release fails the Nextflow 26 config
parse gate (the same old-template failure class as chipseq 2.1.0). Re-run the
consumed-path contract audit when a Nextflow-26-compatible release ships.

## Contents

| File | Purpose |
|---|---|
| `params.example.yaml` | Example params file (rice shown as the example; every path is a placeholder - references are always user-configured, design doc D4) |
| `run-atacseq.sh` | Example HPC launcher wrapping `nextflow run nf-core/atacseq` with the pinned revision |

## Usage

```bash
# 1. Prepare a project params file from the example (fill in real reference
#    paths and the sample sheet); keep it in your project directory, not here.
cp params.example.yaml my-project.params.yaml

# 2. Launch (PBS via the machine-local executor setting; apptainer profile
#    comes from nextflow/conf/hpc.config):
bash run-atacseq.sh my-project.params.yaml
```

Environment variables honored by `run-atacseq.sh` (all optional, with
defaults): `ATACSEQ_REVISION`, `ATACSEQ_PROFILE`, `ATACSEQ_OUTDIR`,
`ATACSEQ_CACHE`.

## ATAC and FAIRE

FAIRE-seq runs through the SAME entry (decision record D4; the legacy
workflow also treated FAIRE identically to ATAC): identical params, with two
interpretive differences documented in `params.example.yaml` - the TSS
profile reads as a sonication/insertion-bias check for FAIRE, and the ext
footprint stage stays disabled. No separate pipeline exists or is needed.

## Local extensions

Capabilities that the upstream pipeline does not provide (organelle QC, QC
gate table, spike-in scaling, TOBIAS footprinting - the pinned pipeline ships
NO TOBIAS) are delivered by `nextflow/pipelines/atacseq-ext/` as default-off
local stages consuming this pipeline's published results. See that pipeline's
README, the decision record, and the design amendment 7 correction.

## Records

- Scientific decision record:
  `nextflow/docs/atacseq-scientific-decision-record.md` (ratified 2026-10-12).
- Validation notes (pin gate, publish contract, extension dry-run baselines;
  execution validation deferred): `nextflow/docs/atacseq-ext-validation-notes.md`.
- Image acquisition recipe (per-source routing): implementation guidance
  Section 3.1 (`plans/2026-10-05-nextflow-migration-implementation-guidance.md`).
