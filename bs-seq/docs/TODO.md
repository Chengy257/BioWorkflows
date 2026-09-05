# bs-seq backlog (TODO)

> Items identified during v0.1 and deferred to a later release; none are scheduled yet.

## 1. Bismark output naming (reconciled against pinned 0.24.0, 2026-09-05/06)

The v0.1 rules declare Bismark's *derived* output names. All items below were
reconciled in 13 WSL real-run rounds against bismark 0.24.0 (bs-seq-pinned
env, `tests/run_test.sh --real-run`); the rules, `TARGETS`, and the README
naming table now declare the observed names:

- [x] Deduplication report: `deduplicate_bismark` writes
  `{sample}.deduplicated.bam` (input `.bam` extension replaced) and
  `{sample}.deduplication_report.txt` into `--output_dir`.
- [x] Nucleotide stats: `bam2nuc` writes `{sample}.deduplicated.nucleotide_stats.txt`
  (input extension replaced by `nucleotide_stats.txt`) into `--dir`; the
  genome-wide composition is `0.index/bismark_genome/genomic_nucleotide_frequencies.txt`.
  NOTE: the `bam2nuc_sample` genome-folder path fix (totals sit one dirname up,
  not three) was applied after the last recorded real-run and is verified at
  DAG level only.
- [x] Per-sample report HTML: `bismark2report --output {sample}` writes
  `{sample}.html`.
- [x] `coverage2cytosine`: `-o {prefix}` + `--merge_CpG` writes
  `{prefix}.CpG_report.merged_CpG_evidence.cov` (plain); the rule gzips it to
  the declared `.cov.gz`.
- [x] `bismark --basename {sample}` (0.24) writes `{sample}_pe.bam`/`_se.bam`
  and `{sample}_PE_report.txt`/`_SE_report.txt`; the rule renames the BAM and
  copies the report to the layout-neutral contract paths, keeping the native
  copies (plus a `{sample}_pe.bam` symlink) for `bismark2summary`.
- [x] `bismark2summary` takes alignment BAMs whose basenames end in `_pe`/`_se`
  and reads the native report beside each; the rule passes absolute
  `_pe`/`_se` symlink paths and `-o bismark2summary` after `cd` into `5.QC/`.
  Dedup/splitting stats are skipped unless reports sit beside the BAM under
  the tool's own names (v0.1 limitation).
- [ ] Per-sample `--parallel` alignment: Bismark 0.24 rejects `--basename`
  together with `--multicore`; v0.1 aligns single-threaded per sample and
  parallelizes across samples via Snakemake (the legacy script used ParaFly
  the same way). Revisit when Bismark lifts the restriction.

## 2. Differential methylation / DMR

- v0.1 stops at per-sample cytosine reports and merged CpG tables
  (`results/5.methylation/`). Differential methylation between conditions
  (methylKit/DSS-style DMR calling; needs a group/batch design in the sample
  table) is out of scope for v0.1 and is the first v0.2 candidate.

## 3. Unit-test suite

- The v0.1 regression is a synthetic-data dry-run (`tests/run_test.sh`). Add a
  unit-test suite (e.g. pytest) for the pure-Python pieces — the workflow
  scripts and the `common.smk` validation/species-merge logic — so refactors
  do not depend on the full snakemake dry-run.

## 4. Environments and CI

- [x] Pinned all-in-one environment solves (validated 2026-09-05 in WSL as
  `bs-seq-pinned`; a pre-existing legacy user env named `bs-seq` was left
  untouched). Record: `mamba env create -f workflow/environment.yaml`.
- [ ] Wire the lightweight CI job (plan Task D3, repository-root
  `.github/workflows/ci.yml`): the `tests/lint.sh` suite plus the
  synthetic-data dry-run regression `bash tests/run_test.sh --reads 2000`.
  CANCELLED for v0.1 by user adjudication (2026-09-06, early scope closure);
  re-open with the next release. Until it lands, `make test` +
  `bash tests/run_test.sh --reads 2000` locally are the reference bar.
- [x] After the first pinned-environment real run, revisit the §1 derived
  names above and the `snakemake --lint` baseline (conda-env advice +
  helper-style warnings are filtered in `tests/lint.sh`); the §1 names are
  reconciled (see above); the lint baseline was re-checked green under
  snakemake 7.32.4.

## 5. Documentation-sync backlog (final-review findings, archived by user adjudication 2026-09-06)

The independent final review of the carrier branch (fix-first verdict,
receipt 2026-09-05) verified the executable name chain (rules / `TARGETS` /
test assertions / the README naming table above / TODO §1) as fully
reconciled, and found scope discipline, red lines, and verification evidence
all clean. Its remaining findings are user-visible documentation that still
carries the PRE-reconciliation name chain — deferred to TODO by explicit user
adjudication (interrupt all execution; archive, do not fix). Sync these on
the next doc pass:

- [ ] `README.md` (~L134-203, outside the naming table): the results
  quick-reference, QC notes, and known-limitations sections still use the old
  names (`{sample}.bam.deduplicated.bam`, `{sample}.dedup_report.txt`,
  `{sample}.nucleotide_stats.txt`, `{sample}.CpG_merged.tsv.gz`,
  `{sample}_seq_context.html`, `genomic_nucleotide_totals.txt`) and say the
  names are "to be confirmed at Phase D" — replace with the reconciled chain
  of the naming table.
- [ ] `docs/user-guide.md` (~§4/§7/§8/FAQ): the directory trees, results
  quick-reference, and FAQ carry the same stale chain; the FAQ also explains
  the obsolete ".bam.deduplicated infix" and links to the README table that
  has since been reconciled.
- [ ] `workflow/rules/align.smk` / `methylation.smk` / `common.smk` header
  comments: still narrate the pre-reconciliation chain and "Phase D must
  confirm" caveats in places.
- [ ] `example/config.yaml` `merge_cpg` comment: references the old
  `{sample}.CpG_merged.tsv.gz` name.
- [ ] `srna-seq/CHANGELOG.md`: record the D2-window behavior changes
  (multiqc `--filename` pin in `meta.smk`; R1/R2 trim artifact guard renames
  in `upstream.smk`) per its own CONTRIBUTING rules.
- [ ] `seclip-seq/workflow/environment.yaml` comment: references
  `docs/env-validation.md`, which was never created (D1 deliverable
  superseded by CHANGELOG records) — point the comment at the CHANGELOG.
