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

## 2026-10-06/07 - real container execution debugging (runs 2-21), PAUSED

Owner scope decision: HPC-only ("不管 WSL 侧，全部以本机为主"); owner then
paused ("暂停任务，报告目前的进展") at run 21.

Images: 7/7 seeded (6 + homer 5.1, verified to carry findMotifsGenome.pl).
Runtime offline pulls from the seeded blob cache confirmed in the logs
(no network at run time).

Nextflow 26.04.6 runtime traps found and fixed (all recorded for
re-derivation avoidance):

1. `collect()` flattens tuples by default - value-list aggregation after a
   tuple stream must use `.map { ... }.toList()` (organelle/spikein summary
   inputs).
2. Workflow `take` accepts channels only; scalars are read from `params`
   inside subworkflows instead of being passed as val arguments.
3. A single-element `tuple(map)` makes downstream `meta.field` a cross-tuple
   spread list - emit bare maps for single-value payloads (DiffBind contrast).
4. Shell globs in script templates must stay quoted when the helper expands
   them itself (`--dup-metrics 'dup/*'`); unquoted, the shell expands them
   into many argv words and argparse rejects the positionals.
5. `join(..., remainder: true)` pads a missing right-hand side with null
   (4-tuple with null) rather than shortening the tuple - handle null
   explicitly before path inputs.
6. `join` pairs duplicate keys ONE-TO-ONE: a shared control sample joined by
   sample id reaches only the first declaring group. Shared-control BAMs are
   resolved synchronously in the workflow body and passed per group instead.
7. Module-level `bin/` auto-staging did not fire in this layout - helper
   scripts moved to the pipeline-level `bin/` (auto-staged + PATH'd).
8. A lost `meta.control` field silently routed the no-control group into the
   normalized SEACR branch - carry the field through every channel hop.

Scientific routing decision implemented: SEACR `norm` normalization requires
a control bedGraph; groups without a control route to a `non`-mode alias
(SEACR_CALLPEAK_TH) with the numeric FDR threshold.

Real execution state at the pause (all in real containers, local executor,
throttled head-node config; 39 tasks incl. 5 cached):

- PASS with published outputs: organelle QC, QC gates, spike-in (summary +
  5 rescaled bigWigs), IDR (both groups, 15 reproducible peaks for WT),
  HOMER motifs (knownResults + homerMotifs), SEACR for groups WITH control
  (H3K27ac_WT narrowPeak; MUT convert was killed mid-run at the pause).
- Remaining for a full PASS (resume points, in order):
  1. RUN_DIFFBIND has never completed - submitted then killed at each
     earlier abort; first resume lets it finish (bioconductor-diffbind
     container already seeded).
  2. SEACR_CONVERT (H3K27ac_MUT) - trivial follow-on after CALLPEAK (done).
  3. SEACR_CALLPEAK_TH (INPUT, the no-control group): SEACR's empirical-FDR
     mode degenerates on toy data - it derives a threshold above the max
     feature AUC (threshold.txt = max AUC), the thresholded set is empty and
     its internal awk divides by zero. Tool minimum-data edge, NOT a pipeline
     bug. Options on resume: (a) give the fixture INPUT real peak/background
     separation; (b) skip SEACR for groups whose only treat IS the shared
     control (scientifically defensible - peaks are not called on the input
     alone) - owner to choose.
- Task paused by owner request; no processes left running.
