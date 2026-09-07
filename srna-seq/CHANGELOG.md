# Changelog

All notable changes to this project are documented in this file. Format based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.1.0] - 2026-09-05

Initial scaffold of the srna-seq subproject (small-RNA cascade filter + quantification, Snakemake 7): config layer, scheduler profiles, boilerplate, and the synthetic test-data contract. Workflow rules, launcher, and tests land in subsequent tasks.

### Changed (2026-09-07, D2 window: multiqc filename pin, trim artifact guard renames, real-run verification)
- multiqc rule pins the report filename on the command line with `--filename multiqc_report.html` (`workflow/rules/meta.smk`): multiqc >= 1.21 otherwise renames the report from the config `title:`, which would miss the declared output `results/5.QC/multiqc/multiqc_report.html`.
- trim rule guard-renames all three trim_galore artifacts to the contract paths `{sample}_trimmed.fq.gz`, `{sample}_trimming_report.txt`, and `fastqc/{sample}_trimmed_fastqc.zip` (`workflow/rules/upstream.smk`): trim_galore 0.6.x names every output after the input basename, so `_R1`-style raw names (e.g. `s1_R1.fq.gz` -> `s1_R1_trimmed.fq.gz`) left the declared outputs missing and failed the job. The R1 fix batch (e3ff082) added the trimming-report rename; the R2 follow-up (03178d6) extended it to the trimmed FASTQ and the FastQC zip. Canonical-name inputs reduce to same-file moves; the `mv ... 2>/dev/null || true` guards keep them non-fatal.
- Real-run verification (2026-09-07, WSL, srna-seq conda env): end-to-end regression `bash tests/run_test.sh --real-run --reads 2000` passed with all six EXPECTED output assertions PASS (`results/4.expression/miRNA/miRNA_counts.tsv`, `miRNA_RPM.tsv`, `all_classes_counts.tsv`, `results/5.QC/cascade_summary.tsv`, `results/5.QC/multiqc/multiqc_report.html`, `results/5.QC/software_versions.yaml`); the dry-run baseline DAG remains 27 jobs.

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
