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

## 2026-10-08 - runs 22-26: DiffBind completion + PBS end-to-end, ALL GREEN

Owner resumed the task; validation scope still HPC-only.

### RUN_DIFFBIND root causes (two independent blockers, both fixed)

1. **bioconda bioconductor-diffbind:3.20.0 image is missing GenomeInfoDb**
   (and UCSC.utils + GenomeInfoDbData). GenomicRanges/GreyListChIP are
   installed but their hard dependency is not, so DiffBind's blacklist step
   died ("there is no package called 'GenomeInfoDb'"). Fixed by sandbox
   surgery: `apptainer build --sandbox` from the SIF, `R CMD INSTALL` of
   UCSC.utils 1.6.1 -> GenomeInfoDbData 1.2.15 -> GenomeInfoDb 1.46.2
   (Bioconductor 3.22 source tarballs, matching the image's R 4.5.2) into
   /usr/local/lib/R/library, repack to SIF, swap in place at
   ~/soft/apptainer-cache/quay.io-biocontainers-bioconductor-diffbind-3.20.0--r45ha27e39d_0.img
   (original kept as *.img.broken until real-data validation). The repo
   pipeline config is unchanged - the patched image carries the fix.
2. **`dba.contrast(block = NULL)` is an error** ("attribute must be a DBA_
   attribute, a logical vector, or a list of logical vectors"): `block` has
   no default and NULL fails attribute validation; omitting the argument is
   fine. run_diffbind.R now builds the contrast argument list and attaches
   `block` only when the Batch column really carries >= 2 levels.

### Fixture surgery (tests/make_fixtures.sh + tests/make_real_bams.sh)

The 10 kb toy chromosome was itself a third blocker, on two fronts:

- Peak loci at 300 bp spacing chain-merged under DiffBind summit
  recentering (summit_flank 250 -> 500 bp windows): pv.Recenter collapsed
  all 25 regions into ONE consensus interval and pv$called lost matrix
  shape (trace: `pv.merge exit: dim(merged)= 1x3`). Peaks now sit at
  L(k) = 3000 + k * 11500 on a 300 kb chr1 (both generator scripts share
  the formula; narrowPeak column 10 is a realistic interior summit
  offset 200, not the width).
- Poisson-flat toy counts broke DESeq2's dispersion fit ("all gene-wise
  dispersion estimates are within 2 orders of magnitude"). Counts are now
  lognormal per locus (4-60 read-pair baseline, multiplicative sample
  noise), loci 21-23 WT-only / 24-25 MUT-only, so the contrast is
  non-degenerate in both directions.
- Spike contig renamed chrS -> **spike1**: the spike-in summary matches
  contig patterns by substring ("spike", 5 chars) and chrS never matched -
  this was the real cause of the all-NA spikein_summary.tsv (NOT stale
  cache; run 23 recomputed it still-NA before the rename).

Standalone DiffBind chain on the new fixture (patched container): count 25
consensus regions, DESeq2 report OK, 4 significant at FDR<=0.05, Fold in
[-4.4, 4.9] with the expected per-locus direction.

### Run outcomes

- run 25 (local executor, -resume): exit 0. Published under
  tests/fixtures/results-ext/: diffbind/WT_vs_mut/ (results+significant+3
  plots+sessionInfo), seacr/ (WT 17 peaks, MUT 10 peaks via norm path),
  idr/ (15 reproducible peaks per group), spike_in/ (real numbers:
  spike_mapped=10, scale_factor=100000), organelle_qc/, qc_gates/
  (5/5 PASS), homer/.
- run 26 (PBS executor, fresh work dir ~/soft/build-cache/nf-pbs-work):
  exit 0 with 41/41 tasks through qsub (queue workq, 4 cpus / 16 GB per
  task, queueSize 8) - the PBS path is validated end-to-end. Machine-local
  config: ~/soft/biowf-ext-pbs.config.
- **The no-control SEACR edge self-resolved**: with the enriched ctl1
  (uniform background + 3 sharp 40-read clusters at loci 4/14/23), the
  INPUT-as-treat `non`+numeric-FDR path exits 0 and calls a valid minimal
  peak set (1 peak, score 1000) - option (a) "enrich the fixture" is
  effectively realized; SEACR_CALLPEAK_TH no longer needs the
  errorStrategy=ignore override (kept only as belt-and-braces in
  ~/soft/biowf-chipseq-ext-seacr-ignore.config).

Machine-local DiffBind image patch provenance: UCSC.utils 1.6.1 +
GenomeInfoDbData 1.2.15 + GenomeInfoDb 1.46.2 from bioconductor.org/3.22
(3-layer proxy), tarballs + build logs under
~/soft/build-cache/diffbind-fix/.

Fixture dry-run DAG baselines are unchanged (content-only edits).
