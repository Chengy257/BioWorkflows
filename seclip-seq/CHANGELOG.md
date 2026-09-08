# Changelog

All notable changes to this project are documented in this file. Format based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.2.0] - 2026-09-08

### Added

- Optional `reproducible_peaks` stage (default off; closes TODO item 1): per-condition
  cross-sample consensus of the ip-role PureCLIP beds via `bedtools multiinter`
  (`rule consensus_peaks` in `workflow/rules/consensus.smk`); a site is kept when
  present in >= `min_replicates` beds and the support count is written to column 4
  of `results/6.reproducible_peaks/{condition}.consensus.bed`.
- Optional `annotate_peaks` stage (default off; closes TODO item 2): GTF-based peak
  annotation for every peak set (per-sample PureCLIP beds plus each consensus BED).
  `rule gtf_gene_regions` derives gene/exon BEDs and a gene_id/gene_name/gene_biotype
  table from the already-required GTF (new stdlib script
  `workflow/scripts/gtf_to_gene_regions.py`: quoted/unquoted attribute parsing,
  gene spans merged from gene + exon records, 1-based GTF -> 0-based half-open BED);
  the `annotate_sample_peaks` / `annotate_consensus_peaks` rules classify each peak
  with `bedtools intersect -u` (exon / gene overlap) and `bedtools closest -d -t first`
  (nearest gene + signed distance) and merge everything into
  `results/6.annotation/{set}.annotation.tsv` (chrom, start, end, score, nearest_gene,
  nearest_gene_id, distance, feature_class exon|gene|intergenic, gene_biotype) via the
  new stdlib script `workflow/scripts/annotate_peaks.py`.
- Optional `condition`/`role` sample-table columns driving the consensus: the sample
  table now accepts either exactly `[sample_id]` (unchanged; single-column tables keep
  working) or exactly `[sample_id, condition, role]`. `role` is restricted to
  `ip|input`, `condition` follows the sample-id name rules (values become file names);
  `load_sample_table` keeps returning the ordered sample-id list, with the grouping
  exposed as the module-level `SAMPLE_CONDITIONS` / `SAMPLE_ROLES` dicts
  (`workflow/rules/common.smk`).
- IP vs input control for the consensus stage (default off; closes TODO item 3): the
  new `reproducible_peaks.input_control` switch peak-calls the `role: input` samples
  too (without it, once the consensus stage is enabled, only ip samples join the
  PureCLIP target set) and unions their PureCLIP beds per condition (`rule
  input_background`: `bedtools multiinter` with support >= 1, i.e. the plain union;
  column 4 = number of inputs covering the feature) into
  `results/6.reproducible_peaks/{condition}.input_background.bed`. For conditions with
  inputs, `rule flag_input_background` produces the final
  `{condition}.consensus.bed` as BED5 — the W7 columns 1-4 plus the binary
  `in_input_background` flag (1 = the site overlaps the condition's input background);
  the raw pre-flagging consensus is kept as the `{condition}.consensus.raw.bed`
  intermediate (`rule consensus_peaks_raw`). Conditions without input samples keep the
  byte-identical W7 BED4 consensus pipeline (`rule consensus_peaks`, scoped to them by
  per-rule wildcard constraints). `reproducible_peaks.filter_by_input` (requires
  `input_control`) adds `rule filter_input_background`, writing
  `results/6.reproducible_peaks/{condition}.consensus.filtered.bed` (flagged sites
  dropped, BED4 again); the consensus annotation consumes the filtered BED under
  `filter_by_input`, otherwise the flagged BED5 (`annotate_peaks.py` already accepts
  arbitrary peak-column counts, so the flag column rides through the annotation
  machinery and is dropped per its TSV column contract).
- Parse-time validation for the input-control switches (aggregated in
  `validate_config`): both keys must be booleans; `input_control` requires the
  condition/role sample-table columns (and warns while `reproducible_peaks.enabled`
  is false); `filter_by_input` requires `input_control`; a condition with input
  samples but no ip samples is a hard error (no consensus to build or flag); an ip
  condition without inputs warns (the background is simply absent and its consensus
  stays unflagged).
- New scheduler-resource defaults (`RESOURCE_DEFAULTS` in common.smk): `input_background`,
  `flag_input_background`, and `filter_input_background` at 1 thread / 2048 MB /
  30 min each (overridable per rule via a project resources.yaml).
- `tests/make_testdata.py --with-inputs`: opt-in input-control test samples (FC_in1 /
  FC_in2, appended to the RNG stream after the ip samples, so the default output stays
  byte-identical); `tests/test_gtf_regions.py` gains a flagged-consensus (BED5) case
  and a 7-column PureCLIP-layout case (28 tests total).
- Parse-time validation for the new switches (aggregated in `validate_config`): both
  sections are optional (absent = disabled, so v0.1 project configs keep working);
  `reproducible_peaks.{enabled,min_replicates}` and `annotate_peaks.enabled` shapes;
  enabled consensus requires the condition/role columns, `callpeak.pureclip=true`,
  and at least `min_replicates` ip samples per condition (hard error, consensus
  support could never be reached); enabled annotation requires `callpeak.pureclip=true`.
- `bedtools=2.31.0` pinned in `workflow/environment.yaml` (seclip had no bedtools
  before; matches the chip sibling pin; invoked bare from PATH like the other tools).
- New scheduler-resource defaults (`RESOURCE_DEFAULTS` in common.smk):
  `consensus_peaks` 1 thread / 2048 MB / 30 min, `gtf_gene_regions` and
  `annotate_peaks` 1 / 4096 / 60 (overridable per rule via a project resources.yaml).
- `tests/test_gtf_regions.py`: pytest unit tests for both scripts (attribute
  parsing, interval extraction, output contracts, merge/feature-class logic, error
  paths; 28 cases after the real-run round additions); `make unit` target
  (`python -m pytest tests -q`) and `make test` now runs check + lint + unit.
- `tests/run_test.sh --consensus` scenario: rewrites the test sample table to the
  condition/role form (both samples one ip condition), enables both stages in the
  test config, dry-runs, and asserts `consensus_peaks`, `gtf_gene_regions`,
  `annotate_sample_peaks`, and `annotate_consensus_peaks` join the DAG.

### Changed

- README results tables now document the PureCLIP output contract (the workflow's
  first written record of it): `results/5.callpeak/{sample}.pureclip.bed` carries 7
  columns — BED6 (chromosome, start, end, site name, crosslink-site score, strand)
  plus a trailing score-attributes field (`[score_CL=...;...]`), corrected from the
  first draft's BED6 assumption by the 40k-read real-run round (see Fixed below).
- Final-review fix: the four bedtools multiinter rules sort their input beds with
  `LC_COLLATE=C sort -k1,1 -k2,2n` (numeric start, stable collation) — the same form
  chip's `callpeak.smk` codified from a real 2026-09-05 deployment failure; the
  lexicographic first draft mis-sorted starts past a 9999/10000 boundary and would
  have corrupted the consensus support counts on real genomes.
- Dry-run job-count baseline unchanged at **23 jobs** with both stages off (the new
  sections are optional in validation; the default test scenario keeps the generated
  single-column sample table). The `--consensus` scenario dry-run adds 5 jobs
  (28 total: consensus 1 + GTF prep 1 + per-sample annotation 2 + consensus
  annotation 1). The new `--input-control` scenario (4 samples: 2 ip + 2 input,
  `input_control` + `filter_by_input` on, both annotation stages off) dry-runs 45
  jobs (18 per-sample jobs for the two extra samples + background/flag/filter/raw
  consensus). `config/resources.yaml` (repository mirror) intentionally not
  extended — unlisted rules keep the built-in defaults per the documented override
  semantics; a project copy can override the new rules by name.
- `tests/run_test.sh`: new `--input-control` scenario (4-sample condition/role
  table, reproducible_peaks enabled with input_control + filter_by_input; dry-run
  asserts `consensus_peaks_raw`, `input_background`, `flag_input_background`, and
  `filter_input_background` join the DAG; `--real-run` additionally asserts the
  background/consensus/filtered outputs); `--consensus` and `--input-control` are
  mutually exclusive; `--consensus --real-run` / `--input-control --real-run`
  additionally need bedtools (documented, not rejected, matching the existing
  convention.

### Fixed

- Real-run round (40k reads, both v0.2 scenarios green end-to-end): PureCLIP 1.3.1
  emits 7 columns (BED6 + a trailing score-attributes field), not BED6 as first
  documented — `annotate_sample_peaks` passes `--peak-cols 7` (the merge script's
  column-count check caught the 14-vs-13 closest mismatch); README, user-guide,
  and the script docs carry the corrected contract and a regression test pins the
  7-column closest layout.
- Real-run round: `tests/make_testdata.py` plants deterministic crosslink hotspots
  (shared ones in every sample, ip-only ones in ip samples) — with uniformly random
  30 bp reads, crosslink sites essentially never coincide across replicates, so the
  consensus filtered to an empty file and the input control flagged every site;
  the scenarios now produce non-empty consensus, background, flagged, and filtered
  outputs under real runs.

## [0.1.0] - 2026-09-05

### Changed (2026-09-05/06, environment solve validation in WSL)
- environment.yaml: python 3.10 -> 3.9 and umi-tools via pip (`umi-tools==1.1.5`) — the classic
  `umi-tools` conda package was removed from bioconda (only the Rust rewrite `umi-tools-rs` remains),
  and the umi-tools sdist cannot build under python >= 3.10 (it bootstraps setuptools 10.0 from 2015).
  `umi_tools` 1.1.5 validated working in the recreated env (STAR 2.7.10b, snakemake 7.32.4).
- multiqc pinned to `--filename multiqc_report.html` (config-file `title:` renames the report otherwise;
  multiqc >=1.21). The seclip real-run smoke was initially not executed (user-adjudicated early scope
  closure 2026-09-06); it was executed and passed on 2026-09-07 (see below).

### Changed (2026-09-07, first end-to-end real-run validation)

- `rule fastq_sort` switched from ea-utils `fastq-sort --id` to `seqkit sort -n` (`seqkit=2.13.0` pinned
  in `workflow/environment.yaml`, ea-utils entry removed): the current bioconda ea-utils build (1.1.2.779)
  ships no `fastq-sort` binary. The runtime resolver and version-collector tool tables follow the rename
  (`fastq_sort`/`fastq-sort` -> `seqkit`).
- STAR rules moved to isolated `--outTmpDir {resources.tmpdir}/STAR_*` directories (cleaned before each
  run — STAR refuses a pre-existing `--outTmpDir`): concurrent index/align jobs previously raced on the
  shared `./_STARtmp` in the working directory, and STAR's FIFO temp files fail outright on non-FIFO
  filesystems (e.g. WSL `/mnt/*` NTFS mounts).
- `star_align` reads its input with `--readFilesCommand cat` when `filter_repeats: true`: the repeats
  filter emits a PLAIN (uncompressed) unmapped Fastx, and `zcat` rejects non-gzip input, so the genome
  alignment silently saw zero reads.
- `umi_dedup` sorts + indexes the STAR BAM before `umi_tools dedup` (pysam fetch requires an indexed,
  coordinate-sorted BAM; STAR writes unsorted BAMs). Dedup preserves input order, so only the final
  index is needed afterwards.
- `fastqc` renames the tool's derived outputs (`{sample}_clean.fqTrTr.sorted_fastqc.*`) to the declared
  `{sample}_fastqc.html/.zip` contract paths (fastqc names outputs after the input file).
- `tests/run_test.sh` assembles its working directory under `/tmp` (overridable via
  `SECLIP_TEST_WORK_DIR`): STAR finishes writing some outputs asynchronously, and on WSL DrvFs mounts a
  reader opening the file within seconds of the writer exiting can observe an empty view until the host
  flushes.
- Real-run regression verified end-to-end in WSL on 2026-09-07 (`bash tests/run_test.sh --real-run`,
  default `--reads 50000`, seclip-seq env): all 23 jobs green, all four EXPECTED output assertions PASS
  (CLIPper absent/auto-skipped path exercised). Note: PureCLIP's HMM parameter learning needs on the
  order of tens of thousands of reads — with `--reads 2000` (the CI dry-run size) it crashes inside
  boost's binomial distribution (`success fraction -nan`); that is a data-size limitation, not a
  workflow defect, so real-run validation uses the script default.

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
