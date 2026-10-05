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

## 2026-10-06 - container seeding + real-fixture preparation (PAUSED HERE)

- Six new images seeded into the shared apptainer cache (blob layer cache,
  offline rebuild at runtime): idr 2.0.4.2, seacr mulled, deeptools mulled
  (all three direct from depot), bioconductor-diffbind 3.20.0, python 3.12.12,
  samtools 1.24 (all three via the proxy quay route). Pull logs:
  ~/soft/build-cache/chipseq-ext-{depot,quay}-pull.log; lists:
  ~/soft/nf-pipelines/chipseq-ext-{depot,quay}.txt (0 failures).
- Machine-local artifacts: `~/soft/biowf-chipseq-ext-overrides.config`
  (rewrites the installed samtools modules' wave blob URL to the quay
  samtools:1.24 image), `~/soft/biowf-ext-smoke.config` (local executor +
  head-node throttle for execution smokes).
- Real-BAM fixtures: `tests/make_real_bams.sh` upgrades the placeholder BAMs
  to real tiny BAMs (reads biased into the fixture peak regions so DiffBind
  counts are non-zero; 10 spike reads on chrS per sample) and writes the
  custom genome FASTA the HOMER stage uses offline. Verified: idxstats on
  the generated BAMs returns chr1:75 / chrS:10 / unmapped:5 for s1.
- First real execution attempt (`-profile test_full` + overrides + local
  executor) STOPPED at runtime sample-sheet validation:
  "control 'ctl1' is declared by more than one group". The fixture sheet
  legitimately shares one control (ctl1) across the WT and MUT groups;
  the validation in `workflows/ext.nf` (groupControl count check) is
  stricter than intended and must allow one control sample shared by
  several treat groups (legacy semantics). **Resume point: relax that
  single check, re-run the smoke.**
- Task paused by owner request 2026-10-06; no processes left running; all
  work committed at the checkpoint below.
