# Changelog

All notable changes to this project are documented in this file. Format based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.1.0] - 2026-09-05

Initial scaffold of the srna-seq subproject (small-RNA cascade filter + quantification, Snakemake 7): config layer, scheduler profiles, boilerplate, and the synthetic test-data contract. Workflow rules, launcher, and tests land in subsequent tasks.

### Added
- Config layer: `config/config.yaml` (repository defaults), `config/config.template.yaml` (annotated project template), `config/resources.yaml` (per-rule scheduler resources), `config/software.yaml` (unified runtime; bowtie / trim_galore / multiqc resolve from PATH), `config/species.yaml` (rice Oryza sativa IRGSP-1.0 preset with per-class sncRNA fastas), `config/samples.csv` (s1 / s2). The ordered `cascade:` list (rRNA / snoRNA / snRNA / tRNA / miRNA / mRNA / rhizo) defines the filter order; an entry whose fasta is empty after the species-preset fallback is skipped entirely.
- `workflow/multiqc_config.yaml` and the four scheduler profiles `workflow/profile/{default,pbs,sge,slurm}` (Snakemake 7 classic `--cluster` interface) plus `workflow/profile/README.md`.
- `tests/make_testdata.py`: deterministic synthetic test-data generator (pure standard library, fixed seed 42, gzip mtime pinned to 0 for byte-reproducibility): 2 x 20 kb genome, 2 x 200 bp rRNA, 3 x 75 bp tRNA, 4 x 21 nt miRNA references, and 21-26 nt SE reads (55/15/10/20% miRNA/rRNA/tRNA/intergenic mix, 0-1 mismatches in miRNA-derived reads, 35% shifted 3' adapter tails) plus a ready-to-run miniature `config.yaml` (rRNA/tRNA/miRNA cascade filled; snoRNA/snRNA/mRNA/rhizo empty -- exercises the empty-fasta skip path).
- Project boilerplate: `Makefile` (`make check` / `lint` / `test`), `LICENSE` (Apache-2.0), `CONTRIBUTING.md`, `.gitignore`, `docs/TODO.md` (v0.2 backlog: DESeq2 differential expression, novel miRNA discovery with miRDeep-P2, unit-test suite).

### Scaffold progress (2026-09-05)

- v0.1.0 scaffold, rules, launcher, runtime resolver, and regression tests
  landed on branch `feature/seclip-srna-bs-seq-v0.1`.
- Baseline dry-run DAG (regression reference): **27 jobs** — synthetic
  dataset (2 samples x 2000 SE reads, species none: rRNA/tRNA/miRNA cascade
  + genome, four classes skipped via empty fasta). Compare against this
  count after refactors.
