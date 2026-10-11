# adopted/rnaseq - nf-core/rnaseq

Thin adoption entry for [nf-core/rnaseq](https://nf-co.re/rnaseq), pinned to
tag **3.27.0** (revision `a1fcdddd3b826fe46eb46f0479f2ff8a7815af05`; parse-
validated under Nextflow 26.04.6 on 2026-10-11). This is **not a fork**: no
pipeline code lives here. The entry carries only the convention layer
permitted by the design doc (Section 8): a params file, reference
configuration example, an HPC example command, and pointers to the scientific
decision and validation records.

## Contents

| File | Purpose |
|---|---|
| `params.example.yaml` | Example params file (rice shown as the example; every path is a placeholder - references are always user-configured, design doc D4). Pins the ratified route: `star_salmon`, UMI opt-in default-off |
| `run-rnaseq.sh` | Example HPC launcher wrapping `nextflow run nf-core/rnaseq` with the pinned revision |

## Usage

```bash
# 1. Prepare a project params file from the example (fill in real reference
#    paths and the sample sheet); keep it in your project directory, not here.
cp params.example.yaml my-project.params.yaml

# 2. Launch (PBS via the machine-local executor setting; apptainer profile
#    comes from nextflow/conf/hpc.config):
bash run-rnaseq.sh my-project.params.yaml
```

Environment variables honored by `run-rnaseq.sh` (all optional, with
defaults): `RNASEQ_REVISION`, `RNASEQ_PROFILE`, `RNASEQ_OUTDIR`,
`RNASEQ_CACHE`, `RNASEQ_EXTRA_ARGS`.

## Published results consumed downstream

The adopted run publishes the artifacts the local extension pipeline and the
future lncRNA workflow consume (layout verified against the pinned
revision's `conf/modules/*.config`):

| Artifact | Location |
|---|---|
| Final (markdup) BAM + index | `{outdir}/star_salmon/{sample}.markdup.sorted.bam(+.bai)` |
| Samtools stats/flagstat/idxstats | `{outdir}/star_salmon/samtools_stats/{sample}.sorted.bam.{stats,flagstat,idxstats}` |
| Picard duplicate metrics | `{outdir}/star_salmon/picard_metrics/{sample}.markdup.sorted.metrics.txt` |
| Native Salmon gene matrix | `{outdir}/salmon/salmon.merged.gene_counts.tsv` (+ `_tpm.tsv`, ...) |
| Per-sample Salmon quant | `{outdir}/salmon/{sample}.salmon_quant/` |

Markdup BAMs publish unconditionally in 3.27.0 (verified against the pinned
source and a real smoke run's output tree), so the extension pipeline's
featureCounts stage and the lncRNA handoff work without extra flags;
`--save_align_intermeds` is only needed to ADD the pre-markdup
`{sample}.sorted.bam` (the legacy pre-dedup counting semantics -
`rnaseq-ext` selects the counting input via `ext_bam_suffix`). Layout
assumes the default `skip_quantification_merge = false`.

## Local extensions

- Capabilities the upstream pipeline does not provide (a featureCounts gene
  matrix with a per-project quantification-route switch, clusterProfiler
  GO/KEGG/GSEA enrichment) are delivered by `nextflow/pipelines/rnaseq-ext/`
  as default-off local stages consuming this pipeline's published results.
- Differential expression runs via `nextflow/adopted/differentialabundance/`
  (DESeq2 backend) on the switch-selected matrix.

## Records

- Scientific decision record:
  `nextflow/docs/rnaseq-scientific-decision-record.md` (ratified 2026-10-11).
- Validation notes: `nextflow/docs/rnaseq-ext-validation-notes.md`.
- Image acquisition recipe (per-source routing): implementation guidance
  Section 3.1 (`plans/2026-10-05-nextflow-migration-implementation-guidance.md`).
