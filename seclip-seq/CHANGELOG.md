# Changelog

All notable changes to this project are documented in this file. Format based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

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
