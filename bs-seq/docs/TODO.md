# bs-seq backlog (TODO)

> Items identified during v0.1 and deferred to a later release; none are scheduled yet.

## 1. Bismark output naming (verify at Phase D real-run validation)

The v0.1 rules declare Bismark's *derived* output names, which depend on the
input file name (`{sample}.bam.deduplicated.bam`) and on tool-version
behavior. Confirm each of the following against the pinned environment
(`tests/run_test.sh --real-run --keep`) and reconcile the declared outputs
and final targets before release:

- [ ] Deduplication report name: `deduplicate_bismark` is expected to write
  `{sample}.bam.dedup_report.txt` beside the deduplicated BAM (the human
  readable report also goes to the log); confirm the exact on-disk name.
- [ ] Nucleotide-stats name: `bam2nuc` is expected to derive
  `{sample}.bam.deduplicated.nucleotide_stats.txt` from the input BAM name;
  confirm the derived name.
- [ ] Per-sample report HTML name: `bismark2report --output {sample}` may
  produce `{sample}.html` rather than `{sample}_seq_context.html`; align the
  declared output and the final targets with whatever the tool writes.
- [ ] `coverage2cytosine` `.gz` suffix behavior: whether a `-o {prefix}.gz`
  output argument makes the tool gzip by itself or the extra `gzip -f` step
  is required; keep the fallback until confirmed on the pinned version.
- [ ] `bismark --basename {sample}` support (needs Bismark >= 0.23): confirm
  it yields `{sample}.bam` + `{sample}_report.txt` inside `--od`; documented
  fallback is a post-rename of `{sample}_pe.bam`.
- [ ] `bismark2summary` writes its HTML into the working directory; the
  workflow `cd`s into `5.QC/` to keep `results/` clean (a `--basename`
  alternative exists).

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
