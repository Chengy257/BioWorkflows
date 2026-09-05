# Changelog

All notable changes to this project are documented in this file. Format based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.1.0] - 2026-09-05

### Changed (2026-09-05/06, environment solve validation in WSL)
- environment.yaml: python 3.10 -> 3.9 and umi-tools via pip (`umi-tools==1.1.5`) — the classic
  `umi-tools` conda package was removed from bioconda (only the Rust rewrite `umi-tools-rs` remains),
  and the umi-tools sdist cannot build under python >= 3.10 (it bootstraps setuptools 10.0 from 2015).
  `umi_tools` 1.1.5 validated working in the recreated env (STAR 2.7.10b, snakemake 7.32.4).
- multiqc pinned to `--filename multiqc_report.html` (config-file `title:` renames the report otherwise;
  multiqc >=1.21). The seclip real-run smoke was NOT executed (user-adjudicated early scope closure
  2026-09-06); lint + dry-run regression remain the v0.1 reference bar.

Initial scaffold of the seclip-seq subproject (single-end enhanced CLIP, Snakemake 7): config layer, scheduler profiles, boilerplate, and the synthetic test-data contract. Workflow rules, launcher, and tests land in subsequent tasks.

### Added
- Config layer: `config/config.yaml` (repository defaults), `config/config.template.yaml` (annotated project template), `config/resources.yaml` (per-rule scheduler resources), `config/software.yaml` (unified runtime; CLIPper resolved as an external absolute path via `paths.clipper`, empty = auto-skip), `config/species.yaml` (hsa GRCh38.p14 / GENCODE v44 preset), `config/samples.csv` (FC_rep1 / FC_rep2).
- `workflow/multiqc_config.yaml` and the four scheduler profiles `workflow/profile/{default,pbs,sge,slurm}` (Snakemake 7 classic `--cluster` interface) plus `workflow/profile/README.md`.
- `tests/make_testdata.py`: deterministic synthetic test-data generator (pure standard library, fixed seed 42, gzip mtime pinned to 0 for byte-reproducibility): 2 x 20 kb genome, 20 genes (gene + exon GTF), 3 snRNA-like 300 bp repeats, and SE UMI reads (10 N UMI + 30 bp genomic/repeat/random body + shifted 3' adapter tail, 1% substitution errors, 60/25/15% source mix) plus a ready-to-run miniature `config.yaml` (tiny STAR indices, CLIPper disabled).
- Project boilerplate: `Makefile` (`make check` / `lint` / `test`), `LICENSE` (Apache-2.0), `CONTRIBUTING.md`, `.gitignore`, `docs/TODO.md` (v0.2 backlog: cross-sample reproducible peaks, peak annotation, IP vs input, ea-utils fallback).

### Scaffold progress (2026-09-05)

- v0.1.0 scaffold, rules, launcher, runtime resolver, and regression tests
  landed on branch `feature/seclip-srna-bs-seq-v0.1`.
- Baseline dry-run DAG (regression reference): **23 jobs** — synthetic
  dataset (2 samples x 2000 SE reads, filter_repeats=true, PureCLIP on,
  CLIPper absent/auto-skipped). Compare against this count after refactors.
