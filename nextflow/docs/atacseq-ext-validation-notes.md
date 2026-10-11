# atacseq-ext validation notes

Validation record for the Stage-3 ATAC/FAIRE route (adopted nf-core/atacseq +
`nextflow/pipelines/atacseq-ext/`). Companion to
`atacseq-scientific-decision-record.md` (ratified 2026-10-12). Dates are
absolute; every entry states what ran and what was observed.

## 2026-10-12 - pin gate + publish-contract extraction

- nf-core/atacseq tag **2.1.2** FAILS the Nextflow 26.04.6 parse gate:
  `nextflow config -flat -profile test` dies with "Unable to parse config
  file" (the same old-template failure class as chipseq 2.1.0).
- nf-core/atacseq **dev HEAD `6a9307709662bad3e7bcc10aa873f7e81fdd1646`**
  (v2.2.0dev, 2026-09-30, `nextflowVersion = '!>=25.10.4'`) parses: exit 0.
  Pinned per the decision record D1; the revision source was fetched to
  /tmp for the contract extraction (not a repo fork).
- Publish contract extracted from the PINNED source (`conf/modules.config`):
  filtered BAMs `{outdir}/{aligner}/merged_library/{id}.mLb.clN.sorted.bam(+.bai)`
  publish UNCONDITIONALLY (no `--save_align_intermeds` caveat on this
  route); stats/flagstat/idxstats under `merged_library/samtools_stats/`
  (both clN and mkD prefixes); picard metrics under `merged_library/picard_metrics/`;
  MACS3 peaks under `merged_library/macs3/{narrow_peak|broad_peak}/`;
  consensus BED `.../{type}/consensus/consensus_peaks.mLb.clN.bed`;
  bigWigs under `merged_library/bigwig/`; unmapped FASTQs only with
  `--save_unaligned` (default false). Defaults: `aligner = 'bwa'`,
  `narrow_peak = false`, `min_reps_consensus = 1`.
- Source-verification corrections recorded in the decision record:
  NO TOBIAS anywhere in atacseq (2.1.2 or dev; design amendment 7 premise
  wrong - resolved by the D2 override ext port), and UMI extraction is
  hardcoded OFF in the main workflow (no opt-in flag; D5 correction).

## 2026-10-12 - atacseq-ext local validation

- Config syntax: `nextflow config -flat -profile test` and `-profile
  test_full` exit 0 for `nextflow/pipelines/atacseq-ext/` under Nextflow
  26.04.6.
- Fixtures: `tests/make_fixtures.py` (fixed seed 20261012) mirrors the
  pinned publish layout - 4 samples (ATAC_WT a1/a2, ATAC_MUT b1/b2, no
  controls), 200 kb chr1 + chrc/chrm organelle contigs + spike1, idxstats
  with organelle/spike rows, parseable flagstat/picard metrics, a 5-region
  consensus BED under `macs3/broad_peak/consensus/`, toy genome FASTA(+.fai)
  and a two-motif JASPAR file. BAM/BAI files are empty placeholders
  (dry-run only resolves names).
- Dry-run baselines: `-profile test` (all stages off) - real run exits 0
  with 0 tasks. `-profile test_full` (all four stages on,
  rescale_bigwigs=true) - `-preview` resolves 11 process slots over 9
  processes (SAMTOOLS_IDXSTATS x2 subworkflow slots, ORGANELLE_SUMMARY,
  QC_GATES_SUMMARY, SPIKEIN_SUMMARY, SPIKEIN_RESCALE, TOBIAS_MERGE_GROUPS,
  TOBIAS_ATACORRECT, TOBIAS_SCOREBIGWIG, TOBIAS_BINDETECT). Per the
  Stage-1 lesson `-preview` prints progressive snapshots and cannot verify
  skips or task multiplicity; the expected task counts at execution are
  24 (8 idxstats, 1 organelle summary, 1 gates, 1 spikein summary, 4
  rescale, 2 merge, 2 ataccorrect, 2 scorebigwig, 1 bindetect) - to be
  confirmed by the deferred real-container fixture run.
- pytest: `tests/unit/` 13 tests green (ms env) pinning the ported helper
  semantics (organelle/spike contig pattern convention incl. the
  spike~spike1 substring rule, QC-gates longest-prefix join + thresholds).
  Full local sweep: `nextflow/tests/static-checks.sh` green under the `nf`
  env (config x3 pipelines, YAML 27 files, doc links, py_compile) and
  `--unit-only` under `ms` (59 tests: 27 chipseq + 19 rnaseq + 13 atacseq).

## Known limitations at this point (deferred, not blockers)

- DEFERRED EXECUTION VALIDATION (owner decision 2026-10-12, design
  amendment 8): the adopted test-profile smoke on the HPC (public test
  data + ~30 images via lazy pulls), real-container fixture execution of
  the ext stages (needs a make_real_bams.py-style helper and seeded
  samtools/tobias images - note BINDetect's sklearn-NaN trap on degenerate
  toy signal), and the exact all-on task-count confirmation. versions.yml
  rows follow the ledger rule - added only after those runs.
- The vendored tobias modules carry upstream container
  `tobias:0.17.5--py310h3479294_0` (depot + quay); not yet pulled.
- The pin is a dev HEAD; re-run the consumed-path contract audit when a
  Nextflow-26-compatible atacseq release ships (Stage-1 watch-item rule).
- Development-box note: during implementation this CentOS 7 machine
  intermittently served stale file contents to freshly spawned processes
  (md5-verified writes invisible to next readers) and an orphaned
  directory entry had to be renamed away before deletion worked; the
  tobias modules are therefore stored under `modules/local/` and every
  file was verified by md5 after writing. Re-verify the tree on the HPC
  at the deferred execution round.
