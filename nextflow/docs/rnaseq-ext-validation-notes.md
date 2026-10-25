# rnaseq-ext validation notes

Validation record for the Stage-2 RNA-seq route (adopted nf-core/rnaseq +
`nextflow/pipelines/rnaseq-ext/` + adopted nf-core/differentialabundance).
Companion to `rnaseq-scientific-decision-record.md` (ratified 2026-10-11).
Dates are absolute; every entry states what ran and what was observed.

## 2026-10-11 - pin gate + publish-contract extraction

- nf-core/rnaseq tag **3.27.0** (revision `a1fcdddd3b826fe46eb46f0479f2ff8a7815af05`)
  parses under Nextflow 26.04.6: `nextflow config -flat -profile test` exit 0.
  No dev-HEAD fallback needed (unlike the chipseq 2.1.0 precedent).
- nf-core/differentialabundance tag **2.0.0** (revision `30ed7741fc`) parses
  under Nextflow 26.04.6 the same way.
- Publish contract extracted from the PINNED source (`conf/modules/*.config`,
  not docs):
  - markdup BAMs `{outdir}/{aligner}/{id}.markdup.sorted.bam(+.bai)` -
    publish UNCONDITIONALLY in 3.27.0 (no `--save_align_intermeds` needed for
    them; that flag only adds the pre-markdup `{id}.sorted.bam`).
  - samtools stats `{outdir}/{aligner}/samtools_stats/{id}.sorted.bam.{stats,flagstat,idxstats}`.
  - Salmon merged matrix `{outdir}/salmon/salmon.merged.gene_counts.tsv` (+ variants).
  - RSeQC strandedness is NOT published by default (publishDir disabled) -
    the ext stage runs its own inference.
- D2 gate verified from the pinned DA source (not docs):
  - `--differential_method deseq2` is the DEFAULT backend (edgeR absent;
    legacy is DESeq2-only).
  - `deseq2_differential.R` validates blocking variables against the
    samplesheet and appends them to the design formula
    (`~ <variable> + <blocking...>`) - batch correction supported.
  - contrasts YAML `comparison: [<variable>, <reference>, <target>]`;
    fold signs are target vs reference (adapter writes control as reference).
  - matrix contract: the native rnaseq Salmon matrix is DA's own test input
    (zero conversion); the featureCounts matrix rides the adapter.

## 2026-10-11 - HPC smoke run (adopted rnaseq 3.27.0), ALL GREEN

- `nextflow run nf-core/rnaseq -revision 3.27.0 -profile test,apptainer`,
  head node, throttled local executor (6 cpu / 24 GB machine-local smoke
  config), work dir under `~/soft/nf-work/rnaseq-smoke/` (never the repo -
  apptainer hardlink lesson).
- Image seeding: `nf-core pipelines download` CANNOT render the 3.27.0
  config (`nextflow config -o json` fails with "No such variable: task" -
  the pipeline-level default `publishDir` closure in `nextflow.config`
  references `task`/`meta`; 10 amend retries all failed at the same step),
  and the wave blob URLs are not pullable via `apptainer pull` directly.
  Route taken instead: LAZY PULLS AT SMOKE TIME - Nextflow fetched every
  executed image into the shared cache itself (`apptainer.enabled` +
  proxy env in the launcher shell). This mirrors the Stage-1 finding that
  runtime pulls from the wave endpoints work; the per-image pre-pull
  script route is retired for rnaseq. 50 unique images were resolved by
  `nextflow inspect -format json -profile test,apptainer` (78 processes);
  the clara-parabricks GPU image never executes in the test profile and
  was not pulled.
- Outcome: pipeline completed successfully; published outputs verified:
  `multiqc/multiqc_report.html`, `star_salmon/{id}.markdup.sorted.bam(+.bai)`,
  `star_salmon/picard_metrics/`, `salmon/salmon.merged.gene_counts.tsv` (+
  `gene_tpm.tsv` - note the underscore variant), `salmon/deseq2_qc/`.
  The consumed contract above matches the real output tree.

## 2026-10-25 - review fix batch (rnaseq-ext P0 6-8, plans/2026-10-14 code review)

- **P0-6 (silent strand corruption) fixed via the S2 convention**: the
  duplicated `-s` injection is gone - the workflow layer now sets
  `meta.strandedness` (forward/reverse/unstranded, mapped from the RSeQC
  legacy codes 1/2/0 or the explicit param) and the vendored
  subread/featurecounts module derives `-s` from it as the single source of
  truth (identical to the lncrna-nf copy; both consumers now share the
  convention). Previously getopt last-wins made every project count
  unstranded while `*.strand.txt` claimed otherwise.
- **P0-7 (>= 2 contrasts broke enrichment) fixed**: `ENRICHMENT_EXT` now
  consumes ONE channel of `(meta, de_table)` tuples (PREPARE_ENRICHMENT
  _INPUTS input reshaped to the nf-core tuple form) - previously a raw Java
  list staged ALL tables in every task and argparse failed on the extras.
  The fixture grew a second contrast (drugB group, 6 samples) to pin the
  pairing permanently.
- **P0-8 / S1 fixed**: `merge_featurecounts.py` accepts the BAM-basename
  count-column labels featureCounts really emits (`<id>.markdup.sorted.bam`)
  while still rejecting other samples' columns; unit fixtures emit the real
  headers (plus a bare-id acceptance case). Previously every real
  featurecounts_quant run aborted on the header check.
- Baselines re-derived by REAL local stub execution (the safe method from
  the lncrna-nf 2026-10-25 entry): `test` (all-off) exit 0, **0 tasks**;
  `test_full` exit 0, **29 tasks over 7 processes** (RSEQC_INFEREXPERIMENT 6,
  PARSE_STRANDEDNESS 6, SUBREAD_FEATURECOUNTS 6, MERGE_FEATURECOUNTS 1,
  PREPARE_ENRICHMENT_INPUTS 2, RUN_ENRICHMENT 6, RUN_GSEA 2). The previous
  18-task baseline encoded the 1-contrast/4-sample fixture; the change is
  fixture growth plus the per-contrast pairing fix, not new stages.
- pytest 19 -> 20 (BAM-basename header accept/reject cases).
- Same-day P2 sweep: the merged-matrix column order now rides the adopted
  sheet order via a value channel (previously task-completion order
  contradicted the helper's documented contract); a paired sheet with
  featurecounts_quant.paired=false warns (fragment double-counting footgun);
  the README commands drop the nonexistent `-profile apptainer` and use
  `-stub-run` instead of `--preview` (on the HPC --preview submits REAL
  cluster jobs via the machine-local executor config).
- Baselines unchanged by the sweep (test 0 tasks; test_full 29 tasks over
  7 processes, re-run after the changes).

## 2026-10-11 - rnaseq-ext dry-run baselines (Nextflow 26.04.6, fixtures)

- Fixtures: `tests/make_fixtures.py` (fixed seed) mirrors the pinned
  publish layout; 4 samples (WT_1/WT_2 control, drugA_1/drugA_2), 40 genes
  with a planted DE signal, one flat DA table
  `da_out/tables/differential/drugA_vs_control.deseq2.results.tsv` (DA
  publishes flat under `tables/differential/`).
- `-profile test` (all stages off): exit 0, **0 tasks**.
- `-profile test_full` (quantification_route=featurecounts + enrichment):
  exit 0, **18 tasks over 7 processes** (RSEQC_INFEREXPERIMENT 4,
  PARSE_STRANDEDNESS 4, SUBREAD_FEATURECOUNTS 4, MERGE_FEATURECOUNTS 1,
  PREPARE_ENRICHMENT_INPUTS 1, RUN_ENRICHMENT 3, RUN_GSEA 1). Per the
  Stage-1 lesson, `-preview` prints progressive snapshots and cannot verify
  skips; real-container fixture execution is the follow-up check.

## 2026-10-11 - unit tests + helper fixes

- `tests/unit/` pytest suite (19 tests, run by the local static-checks
  script `nextflow/tests/static-checks.sh`)
  pins: the verbatim legacy strand mapping (margin 0.4, secondstrand -> 1 /
  firststrand -> 2 / unstranded -> 0), raw featureCounts merging (gene order
  anchored to the sample list, NA fill, status aggregation), the DA adapter
  (metadata-ordered subsetting, contrast derivation with control as
  reference, blocking factors, loud failure modes), and the enrichment DEG
  split (legacy thresholds).
- Two helper bugs the tests caught immediately: the strand parser first
  assumed values on their own lines (real infer_experiment.py puts them on
  the label lines, which the legacy `awk $NF` semantics also imply), and
  the merge header parser had to take the sample column as the LAST column
  of the raw featureCounts layout.

## Known limitations at this point (deferred, not blockers)

- Real-container execution of the ext stages on fixtures (quant chain needs
  real BAMs via a `make_real_bams.sh`-style helper; the R enrichment stage
  needs the clusterProfiler container + an OrgDb) - folded into the WP7
  real-data round together with the D1-route dual-matrix validation.
- differentialabundance images are not yet seeded (no validated DA run
  yet; `versions.yml` gets its row only after one, per the ledger rule).
- The `run_enrichment.R`/`run_gsea.R` ports have no local clusterProfiler
  to smoke against on this machine (the legacy R runtime lives elsewhere);
  their real validation is the WP7 enrichment round. KEGG REST offline
  behavior is documented in the module (auto-skip, GO unaffected).
- Adopted-pipeline re-pins: rnaseq 3.27.x point releases and DA 2.x each
  re-run the consumed-path contract audit (Stage-1 watch-item rule).
