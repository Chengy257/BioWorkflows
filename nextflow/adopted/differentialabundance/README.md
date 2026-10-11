# adopted/differentialabundance - nf-core/differentialabundance

Thin adoption entry for
[nf-core/differentialabundance](https://nf-co.re/differentialabundance),
pinned to tag **2.0.0** (revision `30ed7741fc`; parse-validated under
Nextflow 26.04.6 on 2026-10-11). This is **not a fork**: no pipeline code
lives here. The entry carries only the convention layer permitted by the
design doc (Section 8): a params file, an HPC example command, and pointers
to the scientific decision and validation records.

## Why this pipeline (D2 gate, verified from the pinned source 2026-10-11)

The owner requires DESeq2. Differentialabundance 2.0.0 passes every gate
item, verified in `modules/nf-core/deseq2/differential/templates/deseq2_differential.R`
and `workflows/differentialabundance.nf`:

| Gate item | Finding |
|---|---|
| DESeq2 backend | `--differential_method deseq2` is the DEFAULT (also limma/dream/propd; edgeR absent - the legacy workflow is DESeq2-only, so nothing is lost) |
| Batch covariates | blocking variables (contrasts YAML `blocking_factors` or `--limma_block`-style columns) are validated against the samplesheet and appended to the design formula: `~ <variable> + <blocking...>` |
| Contrast semantics | contrasts YAML/CSV with explicit `comparison: [<variable>, <reference>, <target>]`; fold signs are target vs reference |
| Matrix contract | TSV, features as rows (optional gene metadata columns), observations as columns keyed by the samplesheet id column; the native rnaseq `salmon.merged.gene_counts.tsv` is DA's own test input (zero-conversion path) |

Recorded in the rnaseq scientific decision record (D2) and the validation
notes.

## Contents

| File | Purpose |
|---|---|
| `params.example.yaml` | Example params file (every path is a placeholder; pins `differential_method deseq2` and the legacy thresholds as overridable defaults, documents the contrast semantics and reference/target order) |
| `run-differentialabundance.sh` | Example HPC launcher wrapping `nextflow run nf-core/differentialabundance` with the pinned revision |

## Usage

```bash
# 1. Prepare project inputs: samplesheet.csv, gene matrix (see below),
#    contrasts.yaml (semantics documented in params.example.yaml).
cp params.example.yaml my-project.params.yaml

# 2. Launch (PBS via the machine-local executor setting; apptainer profile
#    comes from nextflow/conf/hpc.config):
bash run-differentialabundance.sh my-project.params.yaml
```

Environment variables honored by `run-differentialabundance.sh` (all
optional, with defaults): `DIFFAB_REVISION`, `DIFFAB_PROFILE`,
`DIFFAB_OUTDIR`, `DIFFAB_CACHE`, `DIFFAB_EXTRA_ARGS`.

## Input matrices (quantification-route switch, decision record D1)

- **route = salmon**: feed `{outdir}/salmon/salmon.merged.gene_counts.tsv`
  from the adopted rnaseq run directly - no conversion needed.
- **route = featurecounts**: run the rnaseq-ext featureCounts stage, then
  convert its matrix + the project samplesheet into DA inputs with
  `nextflow/pipelines/rnaseq-ext/bin/prepare_da_inputs.py` (matrix, DA
  samplesheet, contrasts skeleton).

Which route is authoritative is decided per project by the rnaseq-ext
`quantification_route` param; see the rnaseq-ext README.

## Records

- Scientific decision record (D2, D3):
  `nextflow/docs/rnaseq-scientific-decision-record.md` (ratified 2026-10-11).
- Validation notes: `nextflow/docs/rnaseq-ext-validation-notes.md`.
- OrgDb-based GO/KEGG enrichment is NOT provided by DA (`--functional_method`
  is GMT/GSEA-oriented) - it runs as the rnaseq-ext `enrichment` stage (D3).
