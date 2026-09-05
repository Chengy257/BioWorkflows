# Changelog

All notable changes to this project are documented in this file. Format based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.1.0] - 2026-09-05

Initial scaffold of the bs-seq subproject (Bismark bisulfite-seq calling + QC, Snakemake 7): config layer, scheduler profiles, boilerplate, and the synthetic test-data contract. Workflow rules, launcher, and tests land in subsequent tasks.

### Added
- Config layer: `config/config.yaml` (repository defaults), `config/config.template.yaml` (annotated project template), `config/resources.yaml` (per-rule scheduler resources), `config/software.yaml` (unified runtime; bismark / bowtie2 / trim_galore / multiqc resolve from PATH), `config/species.yaml` (rice Oryza sativa IRGSP-1.0 and human Homo sapiens GRCh38.p14 genome presets), `config/samples.csv` (s1 / s2). `trim.enabled: false` restores the legacy raw-read alignment; `methylation_extractor.buffer_frac` derives the extractor `--buffer_size` from the rule's `mem_mb`.
- `workflow/multiqc_config.yaml` and the four scheduler profiles `workflow/profile/{default,pbs,sge,slurm}` (Snakemake 7 classic `--cluster` interface) plus `workflow/profile/README.md`.
- `tests/make_testdata.py`: deterministic synthetic test-data generator (pure standard library, fixed seed 42, gzip mtime pinned to 0 for byte-reproducibility): 2 x 20 kb genome and paired-end 100 bp reads simulated post-bisulfite (read 1 = C->T-converted forward substrate, read 2 = G->A-converted reverse-complement substrate, 150-300 bp fragments, 1% sequencing errors, 30% of pairs carry a 3' adapter tail beyond the 100 bp genomic body so trimming recovers an alignable read) plus a ready-to-run miniature `config.yaml` (relative `ref/` path, full repository key set).
- Project boilerplate: `Makefile` (`make check` / `lint` / `test`), `LICENSE` (Apache-2.0), `CONTRIBUTING.md`, `.gitignore`, `docs/TODO.md` (Bismark derived output names to verify at Phase D real-run validation; v0.2 backlog: differential methylation / DMR, unit-test suite).
