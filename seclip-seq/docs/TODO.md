# seclip-seq backlog (TODO)

> Items identified during v0.1 and deferred to a later release; none are scheduled yet.

## 1. Cross-sample reproducible peaks

- v0.1 calls peaks per sample (PureCLIP, plus CLIPper when configured). A reproducible peak set across replicates (IDR-style ranking or overlap-based merging of `results/5.callpeak/*.bed`) is out of scope for v0.1 and is the first v0.2 candidate.

## 2. Peak annotation

- Annotate called peaks against gene models / repeat features (GTF-overlap based, or a HOMER/ChIPseeker-style step) so every `results/5.callpeak/*.bed` ships with nearest-gene and biotype tables.

## 3. IP vs input control

- v0.1 has no input/background channel concept. Supporting an paired input control (PureCLIP background estimation from a control BAM, or CLIPper contrast modes) needs a sample-table design decision first (how to declare IP/input pairs).

## 4. fastq_sort tool switch: ea-utils -> seqkit (resolved 2026-09-07)

- [x] `rule fastq_sort` now uses `seqkit sort -n` (seqkit 2.13.0 pinned in
  `workflow/environment.yaml`, ea-utils entry removed): the current bioconda
  ea-utils build (1.1.2.779) ships no `fastq-sort` binary (only fastq-clipper /
  fastq-join / fastq-mcf / fastq-multx / fastq-stats), so the pre-specified
  seqkit fallback was applied during the 2026-09-07 real-run validation. The
  runtime resolver and version-collector tool tables follow the rename
  (`fastq_sort`/`fastq-sort` -> `seqkit`); see CHANGELOG.md.
