# BioWorkflows

[![License: Apache-2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)

A personal modular workflow suite built on the [nf-core](https://nf-co.re/)
ecosystem for Nextflow, for reference-genome analysis on local Linux/WSL and
HPC (apptainer is the single supported container profile on the HPC). The
repository migrated off Snakemake: the legacy suite was retired on 2026-10-26
(see [Retired legacy suite](#retired-legacy-suite)).

## Assay routes

| Assay | Route | Entry points |
|---|---|---|
| ChIP-seq / CUT&Tag | adopted [nf-core/chipseq](https://nf-co.re/chipseq) + local extension stages | `nextflow/adopted/chipseq/` + `nextflow/pipelines/chipseq-ext/` |
| RNA-seq | adopted [nf-core/rnaseq](https://nf-co.re/rnaseq) + local extensions | `nextflow/adopted/rnaseq/` + `nextflow/pipelines/rnaseq-ext/` |
| Differential expression | adopted [nf-core/differentialabundance](https://nf-co.re/differentialabundance) | `nextflow/adopted/differentialabundance/` |
| ATAC-seq / FAIRE-seq | adopted [nf-core/atacseq](https://nf-co.re/atacseq) + local extension stages | `nextflow/adopted/atacseq/` + `nextflow/pipelines/atacseq-ext/` |
| lncRNA identification | local DSL2 pipeline (local-foundation route) | `nextflow/pipelines/lncrna-nf/` |
| small-RNA, BS-seq/methylation, seCLIP/eCLIP | NOT yet migrated - routes not started | - |

Extension stages are default-off and consume the adopted runs' published
results (see each pipeline's README for the input contract). Assay routing:
ChIP-seq and CUT&Tag run the chipseq route; ATAC-seq and FAIRE-seq run the
atacseq route (same entry, decision record D4); RNA-seq runs the rnaseq route;
de novo lncRNA identification runs lncrna-nf downstream of the rnaseq route's
published BAMs.

## Repository structure

```
BioWorkflows/
├── nextflow/                    # the workflow suite
│   ├── conf/                    #   shared execution config (container-only HPC)
│   ├── adopted/                 #   thin adoption entries for upstream nf-core pipelines
│   ├── pipelines/               #   local DSL2 pipelines (chipseq-ext, rnaseq-ext, atacseq-ext, lncrna-nf)
│   ├── docs/                    #   scientific decision + validation records
│   ├── tests/                   #   static-checks.sh (local static checks, see Testing)
│   └── versions.yml             #   tested revisions
└── plans/                       # frozen migration design documents (whitelisted)
```

## Records

- Scientific decision records + validation notes live under
  `nextflow/docs/` (one pair per route); the retired legacy suite's
  capabilities and recovery are recorded in
  `nextflow/docs/legacy-retirement-record.md`.
- Tested revisions: `nextflow/versions.yml` (rows are added only after a
  validated run - the ledger rule).
- The frozen migration mainline and implementation guidance:
  `plans/2026-10-05-nextflow-migration-design.md` and
  `plans/2026-10-05-nextflow-migration-implementation-guidance.md`.

## Requirements

- Nextflow 26.x (validated on 26.04.6; parse gates for the pinned upstream
  revisions are recorded in the validation notes).
- Apptainer on the HPC (the only supported execution profile there);
  container images follow the pinned upstream definitions, with per-process
  overrides for locally patched images.
- Machine-specific paths (executor, cache directories) live in machine-local,
  git-ignored Nextflow configuration - never in this repository.

## Testing

```bash
# One command for all Nextflow entries (former nextflow-static CI job, kept local):
bash nextflow/tests/static-checks.sh
# Unit suites only (no nextflow needed):
bash nextflow/tests/static-checks.sh --unit-only
```

There is no hosted CI: GitHub Actions was retired from the repository on
2026-10-12 (owner decision). Dry-run baselines for every pipeline are recorded
in the validation notes under `nextflow/docs/` and are re-derived with
`-stub-run` after wiring changes.

## Retired legacy suite

The five Snakemake workflows (`rna-seq/`, `chip_cuttag_atac_faire/`,
`seclip-seq/`, `srna-seq/`, `bs-seq/`) and their `shared/` layer were retired
on 2026-10-26 (owner decision; design doc D10). The final state is preserved
twice:

- annotated git tag `legacy-snakemake-suite` (recover with
  `git checkout legacy-snakemake-suite -- <dir>`);
- a local archive copy at
  `/home/chengyu/BioWorkflows_archive/2026-10-26_legacy_snakemake_suite/`
  (HPC machine-local, not in git).

Dropped capabilities and the replacement status of every workflow are
documented in `nextflow/docs/legacy-retirement-record.md`. Note that
small-RNA, BS-seq and seCLIP currently have NO Nextflow replacement - the
archived workflows remain the only implementation for those assays until
their migration routes are started.

## License

Released under the [Apache License 2.0](LICENSE). The workflow directories were previously distributed under MIT licenses (up to rna-seq v0.8.0 and chip_cuttag_atac_faire v0.4.0) and were relicensed to Apache-2.0 as part of the monorepo consolidation.
