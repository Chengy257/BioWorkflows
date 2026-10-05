# chipseq-ext validation notes

Running record of validation for the chipseq-ext local extension pipeline.
Append-only; every entry carries an absolute date.

## 2026-10-06 - skeleton validation (HPC head node)

Environment: Nextflow 26.04.6 (conda env `nf`), apptainer 1.3.2, no proxy.

- `nextflow config -flat -profile test`: OK.
- Dry-run `-profile test` (all stages off): exit 0, **0 tasks** - exercises
  the runtime sample-sheet parse (5 samples, 1 control row), the strict BAM
  resolution (every sheet sample has exactly one BAM; no stray BAMs) and the
  per-directory contract globs.
- Dry-run `-profile test_full` (every stage on): exit 0. DAG resolves the 17
  extension processes across the seven stages:
  `SAMTOOLS_IDXSTATS`, `ORGANELLE_SUMMARY`, `QC_GATES_SUMMARY`,
  `SAMTOOLS_IDXSTATS` (spike-in), `SPIKEIN_SUMMARY`, `SPIKEIN_RESCALE`,
  `IDR`, `IDR_UNION`, `SAMTOOLS_MERGE`, `SAMTOOLS_MERGE_CTL`,
  `DEEPTOOLS_BDGC_RAW`, `DEEPTOOLS_BDGC_RAW_CTL`, `SEACR_CALLPEAK`,
  `SEACR_CONVERT`, `HOMER_FINDMOTIFS`, `DIFFBIND_SHEET`, `RUN_DIFFBIND`.
- Script-level smoke (system python, no containers, fixture data):
  - `organelle_summary.py`: 5-sample table, organelle fraction 0.0000, OK.
  - `spikein_summary.py` + `spikein_factor.py`: spike fraction 0.1111, scale
    factor 10000.0 (1e6/100 spike reads), OK.
  - `qc_gates_summary.py`: mapping rate 0.9000, dup rate 0.2000, organelle
    0.0000, all PASS; sample join on pipeline-declared suffixes, OK.
    (Two real bugs found and fixed by this smoke: prefix-matching of sample
    ids was unreliable across suffix conventions - replaced by explicit
    `--dup-suffix` / `--flagstat-suffix` joining.)
  - `seacr_to_narrowpeak.py`: 6-column SEACR bed to 10-column narrowPeak,
    score clamped to 1000, OK.
  - `idr_union.py`: pairwise merge + distinct-replicate support (support 4
    for a two-pair interval, 2 for a one-pair interval at min-replicates 2),
    OK. (One real bug found and fixed: interval-merge tuple unpacking.)
  - `diffbind_sheet.py`: 4 treat samples across `WT_vs_mut`, control BAM
    attached per group, replicates numbered per condition, OK.

## Known limitations at this point

- No container-level execution yet: `bioconductor-diffbind`,
  `homer:5.1` and `python:3.12.12` images are declared but not pulled
  (acquisition follows the Section 3.1 recipe at the next pull session).
- The end-to-end chain has not run against real adopted-run outputs; the
  published-path contract (README) comes from the dev HEAD source and is
  confirmed only against fixtures.
- DiffBind R code is the legacy WSL-validated logic ported verbatim; its
  first real-data run is the scientific acceptance gate (decision record,
  Section 7).
