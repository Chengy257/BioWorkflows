# RNA-seq Transcriptome Analysis Workflow

A **Snakemake** end-to-end analysis pipeline for bulk RNA-seq (reference-genome based), covering the complete path from raw fastq to differential expression and functional enrichment, with alternative splicing analysis and de novo lncRNA identification.

> **Status**: v0.8.0 — the runtime follows a "single main environment + `software.yaml` overrides" model: existing server Conda environments or the system PATH can be reused directly, and Rscript, R version, R library paths, external tools, and databases are configured in one place; Snakemake no longer creates a dedicated Conda environment per rule. The v0.7 per-rule CPU/memory/walltime scheduler model is retained, now in the dedicated `config/resources.yaml`; output directories are flattened under the configurable `results_dir` (default `results/`).

## Feature overview

| pipeline (selected via the first `run.sh` argument) | Content |
|---|---|
| `upstream` | FastQC -> Trim Galore (incl. FastQC) -> strandedness inference -> STAR alignment -> featureCounts quantification -> TPM/count matrices |
| `deg` | upstream + DESeq2 differential expression -> GO/KEGG enrichment (clusterProfiler; optional aPEAR network plots) -> GSEA -> between-group comparison |
| `as` | upstream (STAR with assembly-optimized parameters) + StringTie transcript assembly/merge/gffcompare -> isoform quantification |
| `lncrna` | upstream + StringTie assembly -> coding-potential filtering (CPC2/CNCI/length/Pfam/NR) -> lncRNA expression matrix |

Runtime environment: Linux + local / SGE (qsub) / SLURM (sbatch) + Snakemake >= 7. It is recommended to maintain a unified `rna-seq` Conda environment on the server, referenced from `config/software.yaml` by prefix/name; the system PATH also works. `run.sh` resolves the runtime, verifies the software/R packages, and selects the scheduler profile in the order SGE -> SLURM -> local.

Scheduler resources are declared per rule in the dedicated `config/resources.yaml` (threads / mem_mb / runtime_min). For example, STAR index defaults to 12 threads / 48 GB / 480 min, STAR alignment to 12 threads / 32 GB / 360 min, and DESeq2 to 4 threads / 16 GB / 240 min. A project can override them by copying `config/resources.yaml` into the project directory as `resources.yaml` (auto-detected by `run.sh`; `RNASEQ_RESOURCES_CONFIG` also works), while `run.sh --memory` / `--runtime` force global cluster-request overrides. The legacy top-level `threads` in `config.yaml` still caps every rule's threads.

## Directory layout

```
rna-seq/
├── run.sh                    # unified launcher (SGE/SLURM/local profile selected automatically or manually)
├── Makefile                  # make check / lint / test development checks
├── workflow/
│   ├── Snakefile             # single entry point (mode selected by the pipeline parameter)
│   ├── rules/                # common / meta / align / quant / deg / as / lncrna
│   ├── scripts/              # R/Python/Shell scripts (stable names, no date suffixes)
│   ├── environment.yaml      # optional all-in-one Conda environment template (not auto-created)
│   ├── profile/              # scheduler profiles (sge / slurm / default)
│   └── multiqc_config.yaml
├── config/
│   ├── config.yaml           # analysis parameters and reference resources
│   ├── config.template.yaml  # project config template (or copy as config.local.yaml, layered last)
│   ├── software.yaml         # main environment, R runtime, tool/database overrides
│   ├── resources.yaml        # per-rule scheduler resources (threads/mem_mb/runtime_min)
│   ├── species.yaml          # osa/hsa species resource presets
│   └── samples.csv           # sample table template
├── example/                  # example project (real sample table reference + one-command start)
├── tests/                    # regression tests: make_testdata.py (synthetic data) / run_test.sh / lint.sh
├── docs/                     # user guide
├── CONTRIBUTING.md           # contribution guide (doc-sync checklist / CHANGELOG requirements)
└── results/                  # run outputs (gitignored; root configurable via results_dir)
```

## Quick start

Full steps: [docs/user-guide.md](docs/user-guide.md). Summary:

1. **Create the project**: put `1.rawdata/` (paired-end `{sample}_1/_2(.fastq|.fq).gz` or `{sample}_R1/_R2(.fastq|.fq).gz`, single-end `{sample}.fastq.gz`/`{sample}.fq.gz`; the first complete pair wins, see `docs/user-guide.md` §4.1) and a sample table (template `config/samples.csv`, must include the control group) in the project directory.
2. **Write the configuration**: copy `config/config.yaml` (or the minimal `config/config.template.yaml`) and `config/software.yaml` into the project directory. The first configures analysis/reference resources, the second the server software environment; with an existing unified Conda environment, usually only `environment.conda_prefix` needs to be filled in. Configuration layers in order: repo defaults -> project config (`-c` or positional) -> `config.local.yaml` in the project directory (auto-detected) or `-l/--extra-config FILE`.
3. **Check the environment**: `bash run.sh -p deg -P /path/to/myproject --software /path/to/myproject/software.yaml --check-software`.
4. **Run**: `bash run.sh deg /path/to/myproject myproject/config.yaml 10`; a `bash run.sh -p deg -P /path/to/myproject -c /path/to/myproject/config.yaml --dry-run` preview of the DAG is recommended first.

> Snakemake no longer creates Conda environments automatically. Without a ready environment, create the recommended all-in-one one with `mamba env create -f workflow/environment.yaml`; with an existing environment, reuse it directly. The environment, Rscript/R library, and tool resolution actually used by each run are recorded in `results/software_versions.yaml`.

## Documentation index

| Document | Content |
|---|---|
| [docs/user-guide.md](docs/user-guide.md) | Data preparation, sample table constraints, configuration, run monitoring, results layout, FAQ |
| [CHANGELOG.md](CHANGELOG.md) | Version change history |

## Known issues (v0.8.0)

Stage 1 fixed all 8 P0 issues; stage 2 completed structure normalization, unified trim/FastQC (the empty MultiQC report is fixed), the volcano plot fix, a parameterized `results_dir` outputs root, species preset mapping, and version recording; stage 3 delivered one-shot regression tests, static checks, CI configuration, and the slurm profile; v0.8.0 further separated the software/R environment from analysis parameters into a unified `software.yaml` runtime.

Remaining major items:

1. KEGG enrichment requires internet access; offline nodes skip it automatically (GO results are unaffected).
2. The `lncrna` pipeline is not covered by the full end-to-end automated regression (it depends on CPC2/CNCI/pfam_scan.pl and large Pfam/NR databases); these dependencies are configured via `software.yaml`, and CNCI's Python2 interpreter can be specified separately.
3. CI lives at the repository root `.github/workflows/ci.yml` (rna-seq job: lint + runtime resolver tests + deg end-to-end regression); the per-subproject `.github/` was removed. End-to-end regression and lint can also be run locally (`make check` / `make lint` / `make test` or `tests/`).

## Testing

```bash
make check                                # runtime resolver + R runtime propagation tests + shell syntax
make lint                                 # static checks: bash/shellcheck/Python/R/snakemake --lint
make test                                 # everything CI runs (= check + lint)
bash tests/test_deggroupcompare.sh        # DEGgroupCompare argument-propagation regression (no R/snakemake needed)
bash tests/run_test.sh --pipeline deg     # end-to-end regression: synthetic data -> dry-run -> run -> 35 assertions
```

Test data is generated deterministically by `tests/make_testdata.py` (2 x 100 kb chromosomes, 60 simulated genes, 2 groups x 2 samples of simulated reads with adapters/sequencing errors/expected differential-gene truth) and is not committed. Assertions cover key file existence, mapping_stat consistency and numeric checks, count matrix shape, and DEG direction checks (against the truth up/down genes).

## Contributing and changes

- Record all changes in `CHANGELOG.md` (for large changes, update the docs before the code).
- Line endings are uniformly LF (see `.gitattributes`); run outputs are never committed (see `.gitignore`).

## License

[Apache-2.0](LICENSE)
