# BioWorkflows

[![License: Apache-2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)

A collection of [Snakemake](https://snakemake.github.io/)-based bioinformatics workflows for reference-genome analysis on HPC clusters (SGE / SLURM / PBS) or local execution.

## Workflows

| Directory | Description | Version | Documentation |
|---|---|---|---|
| [`rna-seq/`](rna-seq/) | End-to-end bulk RNA-seq (reference-based): FastQC → Trim Galore → STAR → featureCounts, plus DESeq2 differential expression, GO/KEGG enrichment and GSEA, isoform assembly with StringTie, and de novo lncRNA identification (CPC2 / CNCI / Pfam / NR). Ships species presets (e.g. rice, human). | v0.8.0 | [README](rna-seq/README.md) · [docs](rna-seq/docs/) |
| [`chip_cuttag_atac_faire/`](chip_cuttag_atac_faire/) | Plant epigenomics pipeline with a single entry `workflow/Snakefile` supporting ChIP-seq, CUT&Tag, ATAC-seq and FAIRE-seq — including mixed-assay projects: bowtie2 alignment, assay-specific deduplication, MACS2 peak calling (narrow / broad / paired-end), ChIPseeker annotation, and FRiP / deepTools QC. | v0.4.0 | [README](chip_cuttag_atac_faire/README.md) · [docs](chip_cuttag_atac_faire/docs/) |

Each subdirectory is a self-contained project with its own configuration templates, user guide (`docs/`), example data (`example/`), regression tests (`tests/`), and changelog.

## Repository structure

```
BioWorkflows/
├── rna-seq/                     # Bulk RNA-seq workflow
│   ├── run.sh                   #   unified launcher / ops CLI
│   ├── workflow/                #   Snakefile, rules, scripts, scheduler profiles, env template
│   ├── config/                  #   analysis config, software env, species presets, sample sheet
│   ├── example/                 #   example project
│   ├── tests/                   #   synthetic-data regression tests + lint
│   └── docs/                    #   user guide, review reports, roadmap
├── chip_cuttag_atac_faire/      # ChIP-seq / CUT&Tag / ATAC-seq / FAIRE-seq workflow (same layout, plus legacy/)
└── LICENSE                      # Apache-2.0
```

## Requirements

- **Linux** — local execution or an HPC cluster (SGE / SLURM; a PBS profile is also provided by the epigenomics workflow).
- **Snakemake 7.x** — reference version 7.32.4. Snakemake 8.x passes parsing/lint but its new cluster-executor semantics are not yet validated for the cluster profiles.
- **Conda / Mamba (recommended)** — each workflow ships a pinned all-in-one environment template (`workflow/environment.yaml`). Existing conda environments or system PATH tools can be reused instead, configured via `config/software.yaml`.

## Quick start

```bash
# RNA-seq — run the "deg" pipeline on a project directory (10 threads)
bash rna-seq/run.sh deg /path/to/myproject myproject/config.yaml 10

# ChIP / CUT&Tag / ATAC / FAIRE — preflight-check the software environment first
bash chip_cuttag_atac_faire/run.sh -P /path/to/myproject --check-software
```

See each project's README for the full step-by-step tutorial (working directory → raw data → sample sheet → config → launch) and the `--dry-run` DAG preview.

## Design notes

- **Unified launcher (`run.sh`)** — scheduler profile auto-detection, preflight software and R-package checks, per-rule or global resource overrides, `--dry-run` preview, and `unlock` helpers.
- **Explicit environment model** — Snakemake never auto-creates per-rule conda environments. A single main environment is resolved at runtime from `config/software.yaml` (conda prefix/name or system PATH, including Rscript and R library paths), and the resolved tool versions are recorded into the results directory on every run.
- **Per-rule cluster resources** — threads / memory / walltime are declared per rule and can be overridden per project or globally (`run.sh --memory`, `--runtime`).
- **Data stays out of version control** — raw reads, references, and run outputs live in project working directories and are excluded via `.gitignore`.

## Testing

Both workflows ship deterministic synthetic-data regression tests, static lint checks, and dry-run DAG validation:

```bash
bash rna-seq/tests/run_test.sh

bash chip_cuttag_atac_faire/tests/run_test.sh    # dry-run regression
python chip_cuttag_atac_faire/tests/run_tests.py # 55 unit tests
bash chip_cuttag_atac_faire/tests/lint.sh
```

## License

This repository is released under the [Apache License 2.0](LICENSE).

The two workflow directories were previously distributed under MIT licenses (up to rna-seq v0.8.0 and chip_cuttag_atac_faire v0.4.0) and are relicensed to Apache-2.0 as part of the monorepo consolidation.
