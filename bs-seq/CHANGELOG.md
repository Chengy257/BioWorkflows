# Changelog

All notable changes to this project are documented in this file. Format based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.1.0] - 2026-09-05

Initial scaffold of the bs-seq subproject (Bismark bisulfite-seq calling + QC, Snakemake 7): config layer, scheduler profiles, boilerplate, and the synthetic test-data contract. Workflow rules, launcher, and tests land in subsequent tasks.

### Added
- Config layer: `config/config.yaml` (repository defaults), `config/config.template.yaml` (annotated project template), `config/resources.yaml` (per-rule scheduler resources), `config/software.yaml` (unified runtime; bismark / bowtie2 / trim_galore / multiqc resolve from PATH), `config/species.yaml` (rice Oryza sativa IRGSP-1.0 and human Homo sapiens GRCh38.p14 genome presets), `config/samples.csv` (s1 / s2). `trim.enabled: false` restores the legacy raw-read alignment; `methylation_extractor.buffer_frac` derives the extractor `--buffer_size` from the rule's `mem_mb`.
- `workflow/multiqc_config.yaml` and the four scheduler profiles `workflow/profile/{default,pbs,sge,slurm}` (Snakemake 7 classic `--cluster` interface) plus `workflow/profile/README.md`.
- `tests/make_testdata.py`: deterministic synthetic test-data generator (pure standard library, fixed seed 42, gzip mtime pinned to 0 for byte-reproducibility): 2 x 20 kb genome and paired-end 100 bp reads simulated post-bisulfite (read 1 = C->T-converted forward substrate, read 2 = G->A-converted reverse-complement substrate, 150-300 bp fragments, 1% sequencing errors, 30% of pairs carry a 3' adapter tail beyond the 100 bp genomic body so trimming recovers an alignable read) plus a ready-to-run miniature `config.yaml` (relative `ref/` path, full repository key set).
- Project boilerplate: `Makefile` (`make check` / `lint` / `test`), `LICENSE` (Apache-2.0), `CONTRIBUTING.md`, `.gitignore`, `docs/TODO.md` (Bismark derived output names to verify at Phase D real-run validation; v0.2 backlog: differential methylation / DMR, unit-test suite).

### Added (2026-09-05: launcher, runtime, tests, docs)

- `run.sh` unified launcher (`BSSEQ_` env prefix family: `BSSEQ_CONFIG` / `BSSEQ_RESOURCES_CONFIG` / `BSSEQ_EXTRA_CONFIG` / `BSSEQ_SOFTWARE_CONFIG` / `BSSEQ_PYTHON`): four scheduler profiles with auto detection, software preflight (`--check-software`), per-rule resource overrides (`--memory`/`--runtime`), retries/latency/throttling options, `--unlock`, `--` passthrough to Snakemake, and the positional form `bash run.sh <project_dir> [config.yaml] [jobs]`.
- `workflow/scripts/runtime_config.py` (software.yaml resolver/preflight over the shared runtime framework; the full bs-seq tool set from bismark to multiqc) and `workflow/scripts/collect_versions.py` (writes `results/5.QC/software_versions.yaml`).
- `workflow/environment.yaml`: pinned all-in-one environment template `bs-seq` (snakemake-minimal 7.32.4 / bismark 0.24.0 / bowtie2 2.5.2 / samtools 1.17 / trim-galore 0.6.10 / fastqc 0.11.9 / multiqc 1.21; created explicitly by the user).
- `tests/lint.sh` (bash -n + shellcheck + py_compile + `snakemake --lint` with the known-baseline filter) and `tests/run_test.sh` (one-command regression: synthetic data -> assembled working directory -> dry-run with DAG assertions / `--real-run` with output assertions -> cleanup, optional DAG SVG).
- Documentation: `README.md` (workflow overview, Bismark output naming table, differences from the legacy pipeline, testing), `docs/user-guide.md` (installation, quick start, config reference, cluster submission, output interpretation, FAQ), and `example/` (two-sample project with a filled config and start guide). `docs/TODO.md` gains an environments-and-CI section.

### Regression baseline (2026-09-05)

- Baseline dry-run DAG (regression reference): **20 jobs** — synthetic dataset
  (2 samples x PE reads, `species: "none"`, miniature `ref/` genome):
  rule all (1), trim_pe (2), bismark_genome_prep (1), bismark_align (2),
  deduplicate (2), bam2nuc_genome (1), bam2nuc_sample (2),
  methylation_extractor (2), coverage2cytosine (2), bismark2report (2),
  bismark2summary (1), multiqc (1), software_versions (1). Compare against
  this count after refactors.
