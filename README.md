# BioWorkflows

[![License: Apache-2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)

A collection of [Snakemake](https://snakemake.github.io/)-based bioinformatics workflows for reference-genome analysis on HPC clusters (SGE / SLURM / PBS) or local execution. All workflows share one architecture: a unified launcher, an explicit software-runtime model, per-rule cluster resources in a dedicated config file, species presets, and a common `shared/` layer.

## Workflows

| Directory | Description | Version | Documentation |
|---|---|---|---|
| [`rna-seq/`](rna-seq/) | End-to-end bulk RNA-seq (reference-based): FastQC → Trim Galore → STAR → featureCounts, plus DESeq2 differential expression, GO/KEGG enrichment and GSEA, isoform assembly with StringTie, and de novo lncRNA identification (CPC2 / CNCI / Pfam / NR). Ships species presets (rice, human). | v0.8.0 | [README](rna-seq/README.md) · [User guide](rna-seq/docs/user-guide.md) |
| [`chip_cuttag_atac_faire/`](chip_cuttag_atac_faire/) | Plant epigenomics pipeline with a single entry `workflow/Snakefile` supporting ChIP-seq, CUT&Tag, ATAC-seq and FAIRE-seq — including mixed-assay projects: bowtie2 alignment, assay-specific deduplication, MACS2 peak calling (narrow / broad / paired-end), ChIPseeker annotation, and FRiP / deepTools QC. | v0.4.0 | [README](chip_cuttag_atac_faire/README.md) · [User guide](chip_cuttag_atac_faire/docs/user-guide.md) |
| [`seclip-seq/`](seclip-seq/) | Single-end enhanced CLIP (eCLIP-style) protein-RNA binding maps: UMI extraction, two-pass adapter trimming, an optional sncRNA/repeats pre-filter, unique STAR alignment, UMI-based deduplication, and PureCLIP peak calling (CLIPper optional), ending in a combined MultiQC report, plus optional replicate-consensus peak annotation (`--consensus`) and IP-vs-input-control (`--input-control`) scenarios. | v0.2.0 | [README](seclip-seq/README.md) · [User guide](seclip-seq/docs/user-guide.md) |
| [`srna-seq/`](srna-seq/) | Small-RNA sequencing with a bowtie1 cascade filter (rRNA → snoRNA → snRNA → tRNA → miRNA → mRNA → optional exogenous classes; order configurable, empty-fasta classes skipped), per-class count and RPM matrices, and optional differential expression (`--deg`) and novel miRNA discovery (`--novel-mirna`) stages. | v0.2.0 | [README](srna-seq/README.md) · [User guide](srna-seq/docs/user-guide.md) |
| [`bs-seq/`](bs-seq/) | Bisulfite sequencing (BS-seq) methylation calling: optional Trim Galore trimming, Bismark alignment, PCR-duplicate removal and methylation extraction with an optional merged CpG table, per-sample and run-level reports, and an optional group-design differential-methylation (`--dmr`) stage. | v0.2.0 | [README](bs-seq/README.md) · [User guide](bs-seq/docs/user-guide.md) |

## Repository structure

```
BioWorkflows/
├── rna-seq/                     # Bulk RNA-seq workflow
│   ├── run.sh                   #   unified launcher / ops CLI
│   ├── workflow/                #   Snakefile, rules, scripts, scheduler profiles, env template
│   ├── config/                  #   config + template, resources, software env, species presets, sample sheet
│   ├── example/  tests/  docs/  #   example project, regression tests, user guide
├── chip_cuttag_atac_faire/      # ChIP-seq / CUT&Tag / ATAC-seq / FAIRE-seq workflow (same layout)
├── seclip-seq/                  # Single-end eCLIP (protein-RNA binding) workflow (same layout)
├── srna-seq/                    # Small-RNA cascade-filter workflow (same layout)
├── bs-seq/                      # Bisulfite methylation-calling workflow (same layout)
├── nextflow/                    # Nextflow migration (in progress; see its README below)
│   ├── conf/                    #   shared execution config (container-only HPC)
│   ├── adopted/                 #   thin adoption entries for upstream nf-core pipelines
│   ├── pipelines/               #   local DSL2 pipelines (chipseq-ext, rnaseq-ext)
│   ├── docs/                    #   scientific decision + validation records
│   ├── tests/                   #   static-checks.sh (local static checks, see Testing)
│   └── versions.yml             #   tested revisions
├── shared/                      # Cross-project layer (consumed in place, never copied)
│   ├── python/                  #   software-runtime framework + version capture
│   └── lib/launcher.sh          #   shared launcher helpers for every project's run.sh
```

## Nextflow migration (in progress)

The repository is migrating to Nextflow/nf-core under a frozen scientific-first
mainline (design + guidance documents under `plans/`). Status:

- `nextflow/adopted/chipseq/` — adopted [nf-core/chipseq](https://nf-co.re/chipseq)
  entry (pinned dev HEAD), with an example params file and HPC launcher.
  The Snakemake `chip_cuttag_atac_faire/` legacy capabilities are delivered by
  `nextflow/pipelines/chipseq-ext/` as default-off extension stages
  (DiffBind, spike-in, IDR, SEACR, HOMER motifs, QC gates, organelle QC).
- `nextflow/adopted/rnaseq/` — adopted [nf-core/rnaseq](https://nf-co.re/rnaseq)
  entry (pinned 3.27.0), Stage-2 route: `star_salmon` with published BAMs,
  UMI opt-in. Local extensions (featureCounts quantification with a
  per-project salmon/featureCounts route switch, clusterProfiler GO/KEGG/GSEA
  enrichment) land in `nextflow/pipelines/rnaseq-ext/`; differential
  expression runs via `nextflow/adopted/differentialabundance/` (DESeq2).
- `nextflow/adopted/atacseq/` — adopted [nf-core/atacseq](https://nf-co.re/atacseq)
  entry (pinned dev HEAD `6a930770`; the 2.1.2 release fails the Nextflow 26
  parse gate), Stage-3 route for ATAC-seq AND FAIRE-seq (same entry, D4).
  Local extensions (organelle QC, QC gates, spike-in scaling, TOBIAS
  footprinting - the pinned pipeline ships no TOBIAS) land in
  `nextflow/pipelines/atacseq-ext/`.
- Assay routing: ChIP-seq and CUT&Tag run the chipseq route
  (`adopted/chipseq` + `pipelines/chipseq-ext`); ATAC-seq and FAIRE-seq run
  the atacseq route (`adopted/atacseq` + `pipelines/atacseq-ext`); RNA-seq
  runs the rnaseq route above.
- Records: `nextflow/docs/chipseq-scientific-decision-record.md` (ratified
  2026-10-06), `nextflow/docs/chipseq-ext-validation-notes.md`,
  `nextflow/docs/rnaseq-scientific-decision-record.md` (ratified 2026-10-11),
  `nextflow/docs/rnaseq-ext-validation-notes.md`, and
  `nextflow/versions.yml` (tested revisions).
- The Snakemake suites below remain the validated daily drivers; retirement
  happens only when the new routes cover real analysis needs (design doc
  Stage 4).

## Shared architecture

All workflows follow the same model (see each project's README for specifics):

- **Output layout** — raw inputs (`1.rawdata/`) stay at the project working-directory root; every derived artifact goes under a configurable `results_dir` (default `results/`) using numbered stage directories (`2.cleandata/`, `3.align/`, `4.peak/` or `4.expression/`, `5.QC/` or `5.DEG/`, ...).
- **Config layering** — repo defaults → project config → `config.local.yaml` in the working directory (auto-detected; later files win). Species presets (`config/species.yaml`) fill unset reference keys.
- **Resources** — per-rule `threads / mem_mb / runtime_min` live in a dedicated `config/resources.yaml` (copy into a project to override); a legacy top-level `threads` still caps everything globally.
- **Unified launcher (`run.sh`)** — scheduler auto-detection, preflight software and R-package checks (`--check-software` / `--check-r`), parse-time config validation with aggregated errors, per-rule or global resource overrides, `--dry-run` DAG preview, `--unlock`.
- **Explicit software runtime** — Snakemake never auto-creates environments. A single main environment is resolved from `config/software.yaml` (Conda prefix/name or system PATH, including Rscript and R library paths); the resolved versions are recorded into the results directory on every run (`software_versions.yaml`).

## Requirements

- **Linux** — local execution or an HPC cluster (SGE / SLURM; a PBS profile is also provided by the epigenomics workflow).
- **Snakemake 7.x** — reference version 7.32.4. Snakemake 8.x passes parsing/lint but its cluster-executor semantics are not yet validated for the cluster profiles.
- **Conda / Mamba (recommended)** — each workflow ships a pinned all-in-one environment template (`workflow/environment.yaml`). Existing conda environments or system PATH tools can be reused instead via `config/software.yaml`.

## Quick start

```bash
# RNA-seq — run the "deg" pipeline on a project directory (10 threads)
bash rna-seq/run.sh deg /path/to/myproject myproject/config.yaml 10

# ChIP / CUT&Tag / ATAC / FAIRE — preflight-check the software environment first
bash chip_cuttag_atac_faire/run.sh -P /path/to/myproject --check-software
```

## Testing

```bash
# rna-seq: runtime resolver tests + lint; full regression needs the analysis environment
cd rna-seq && make check && make lint
bash tests/run_test.sh --pipeline deg        # synthetic data, end-to-end

# chip: 62 unit tests + lint + synthetic-data dry-run (no analysis tools needed)
cd chip_cuttag_atac_faire && make test
bash tests/run_test.sh                       # dry-run DAG regression

# seclip-seq / srna-seq / bs-seq: check + lint + synthetic-data regression (no analysis tools needed)
cd seclip-seq && make test && bash tests/run_test.sh --reads 2000
cd ../srna-seq && make test && bash tests/run_test.sh --reads 2000
cd ../bs-seq && make test && bash tests/run_test.sh --reads 2000
```

There is no hosted CI: GitHub Actions was retired from the repository on
2026-10-12 (owner decision — it was unnecessary during development). All checks
run locally through the commands above, and the Nextflow migration entries are
checked with a single command:

```bash
bash nextflow/tests/static-checks.sh   # former nextflow-static CI job, kept local
```

## License

Released under the [Apache License 2.0](LICENSE). The workflow directories were previously distributed under MIT licenses (up to rna-seq v0.8.0 and chip_cuttag_atac_faire v0.4.0) and were relicensed to Apache-2.0 as part of the monorepo consolidation.
