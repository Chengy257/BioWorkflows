# Changelog

All notable changes to this project are documented in this file. Format based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.2.0] - 2026-09-08

### Added

- Differential-methylation stage (default off): optional `dmr:` config section
  (`enabled` / `control_group` / `qvalue` / `min_diff` / `tile_len` /
  `tile_step` / `min_cpg` / `batch_correction` / `min_cov` / `max_cov`) and a
  parse-time-guarded `workflow/rules/dmr.smk` with a single wildcard-free
  `dmr_methylkit` job (resources: 4 threads / 16000 MB / 240 min, overridable
  via `resources.dmr`). It runs `workflow/scripts/run_dmr.R` (methylKit) over
  every sample's merged CpG table
  (`5.methylation/{sample}/{sample}.CpG_merged.CpG_report.merged_CpG_evidence.cov.gz`,
  read via `pipeline="bismarkCoverage"`): per contrast `{treat}_vs_{control}`
  one pairwise methylKit analysis produces
  `6.DMR/{treat}_vs_{control}_DMC_{all,hyper,hypo}.tsv` and
  `_DMR_tiles.tsv` (columns: chr, start, end, strand, meth.diff in percent
  (treatment minus control), pvalue, SLIM qvalue) plus `6.DMR/DMR_summary.tsv`
  and `sessionInfo.txt`; outputs stay inside `6.DMR/`.
- Sample-table design columns: headers `sample_id`, `sample_id,group`, or
  `sample_id,group,batch` (same name rules as `sample_id`, non-empty values);
  `load_sample_table` still returns the ordered id list and fills
  `SAMPLE_GROUPS` / `SAMPLE_BATCH` in the same parse. Design validation for
  `dmr.enabled: true` (aggregated `WorkflowError`): group column required,
  `control_group` present, every group >= 2 replicates (hard error — methylKit
  cannot fit fewer), `batch_correction: "T"` requires the batch column, and
  `methylation_extractor.merge_cpg` must be true.
- R plumbing: `workflow/environment.yaml` gains `r-base=4.3`,
  `bioconductor-methylkit`, `r-getopt`; `config/software.yaml` gains the `r:`
  section (resolved `Rscript` exported as `BSSEQ_RSCRIPT`); `rules/common.smk`
  gains `RSCRIPT = os.environ.get("BSSEQ_RSCRIPT", "Rscript")`.
- Tests: `tests/make_testdata.py --dmr` emits the 2x2 differential-methylation
  scenario (4 samples, s1/s2 control + s3/s4 treat, group column, dmr-enabled
  miniature config; default output byte-identical to before);
  `bash tests/run_test.sh --dmr` dry-runs it and asserts `dmr_methylkit` in
  the DAG (35 jobs; default scenario stays 20 jobs and asserts the rule is
  absent). `tests/test_common.py` grows to 93 tests: group/batch loader
  coverage, dmr section validation, and design-validation pins.
- Smoke-tested `run_dmr.R` end-to-end against methylKit 1.20.0 in WSL
  (2x2 synthetic tables with low-/high-coverage artifact sites): 299 hypo
  DMCs recovered as simulated, `filterByCoverage` bounds and the batch
  covariate path (`-b T`) verified.

### Changed

- `workflow/environment.yaml` aligned to the standard all-in-one template
  style (standard header comment + block-form `channels:`), keeping every
  existing pin unchanged (per the repository guide this file was due its
  style alignment on its next touch).
- `config/config.yaml`, `config/config.template.yaml`, and
  `example/config.yaml` carry the `dmr:` section (enabled: false);
  `example/samples.csv` documents the optional `group` column (s1/s2
  control; the two-sample example cannot enable dmr as shipped).

### Fixed

- Final-review fix: `rules/common.smk` resolves the sample table once into the
  module-level `_SAMPLE_TABLE` (absolute > project-relative > repository-relative)
  and `rules/dmr.smk` consumes that resolved path instead of the raw
  `config["SampleListFile"]` value — a table declared relative to the repository
  root would previously have failed at job runtime (matches the srna-seq deg
  pattern).
- Real-run round: `--real-run --dmr` (2x2, methylKit inside the workflow) passes
  all assertions including `results/6.DMR/DMR_summary.tsv` after the pinned
  environment gained `r-base` + `bioconductor-methylkit`.

## [0.1.0] - 2026-09-05

### Changed (2026-09-05/06, WSL real-run reconciliation against bismark 0.24.0)
- Bismark tool flag surface reconciled to 0.24: `--genome_folder`/`--output_dir` (not `--genome`/`--od`),
  `-1/-2` mate flags, no `--pe` on bismark (paired-ness inferred), `-p` on deduplicate,
  no per-job `--parallel` (mutually exclusive with `--basename` in 0.24; parallelism across samples via Snakemake).
- Derived-name chain redeclared to the observed names (see README "Bismark output naming" and docs/TODO.md §1):
  `{sample}.deduplicated.bam` + `{sample}.deduplication_report.txt`; extractor outputs
  `{sample}.deduplicated.bismark.cov.gz` / `.bedGraph.gz` / `_splitting_report.txt` / `.M-bias.txt`;
  coverage2cytosine `{sample}.CpG_merged.CpG_report.merged_CpG_evidence.cov.gz`;
  bam2nuc `genomic_nucleotide_frequencies.txt` + `{sample}.deduplicated.nucleotide_stats.txt`;
  bismark2report `{sample}.html`; bismark2summary `bismark2summary.html` (fed `_pe`/`_se` BAM symlinks).
- Genome sentinel `Bisulfite_Genome/GA_conversion/BS_GA.1.bt2` (0.24 bowtie2 conversion naming).
- multiqc pinned to `--filename multiqc_report.html` on the command line (config-file `title:` renames the
  report otherwise; multiqc >=1.21 behavior).
- Real-run status at v0.1 close (user-adjudicated early scope closure 2026-09-06): 13 reconciliation rounds;
  15/20 jobs green in the last recorded run incl. trim/align/dedup/extraction/merge/summary/multiqc; the
  `bam2nuc_sample` genome-folder path fix and the bismark2report naming were applied after that run and are
  DAG-verified (dry-run) but not re-run end-to-end (recorded honestly; see docs/TODO.md §1).

### Changed (2026-09-07, first full end-to-end real-run green)

- `bam2nuc_sample` no longer passes `--genomic_composition <totals>`: bam2nuc has no such option, and Perl
  Getopt abbreviation silently resolves it to `--genomic_composition_only` (genome composition, exit 0, BAM
  never processed — the declared `{sample}.deduplicated.nucleotide_stats.txt` was never written). The
  composition file is auto-detected from the genome folder; the totals input stays as the DAG ordering edge.
- `bismark2report` passes `--output {sample}.html`: `--output` uses the file name VERBATIM (passing
  `{sample}` writes a file literally named `{sample}`, missing the declared `.html` output).
- Real-run regression verified end-to-end in WSL on 2026-09-07 (`bash tests/run_test.sh --real-run
  --reads 3000`, bs-seq-pinned env): all 20 jobs green, all four EXPECTED output assertions PASS; the two
  fixes above close the previously recorded 15/20 gap, and also verify the `bam2nuc_sample` genome-folder
  path fix and the `bismark2report` output name at runtime.
- Documentation synced to the reconciled naming chain (final-review TODO §5 backlog): README results
  quick-reference / QC notes / limitations, user-guide §4/§7/§8 + FAQ Q5, rule header comments, and the
  example config `merge_cpg` comment. The lightweight CI job (lint + `--reads 2000` dry-run regression) is
  wired at the repository root `.github/workflows/ci.yml` (see docs/TODO.md §4).

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
