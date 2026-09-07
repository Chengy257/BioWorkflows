# seclip-seq backlog (TODO)

> Items identified during v0.1 and deferred to a later release.

## 1. Cross-sample reproducible peaks (resolved 2026-09-08, v0.2.0)

- [x] Implemented as the optional, default-off `reproducible_peaks` stage: an
  overlap-based consensus per condition via `bedtools multiinter` over the
  ip-role PureCLIP beds (support >= `min_replicates`, written to column 4 of
  `results/6.reproducible_peaks/{condition}.consensus.bed`), driven by the
  optional `condition`/`role` sample-table columns. IDR-style ranking of
  replicates remains a possible future refinement.

## 2. Peak annotation (resolved 2026-09-08, v0.2.0)

- [x] Implemented as the optional, default-off `annotate_peaks` stage: the
  already-required GTF is parsed into gene/exon BEDs plus a gene attribute
  table (stdlib script `workflow/scripts/gtf_to_gene_regions.py`), and every
  peak set (per-sample PureCLIP beds plus each consensus BED) is classified
  with `bedtools intersect` / `bedtools closest` and merged into
  `results/6.annotation/{set}.annotation.tsv` (chrom/start/end/score,
  nearest_gene/nearest_gene_id/distance, feature_class exon|gene|intergenic,
  gene_biotype). Repeat-feature annotation (beyond the GTF) is future work.

## 3. IP vs input control

- v0.1 has no input/background channel concept. Supporting an paired input control (PureCLIP background estimation from a control BAM, or CLIPper contrast modes) needs a sample-table design decision first (how to declare IP/input pairs).
- The v0.2 sample-table extension already reserves the `role` column
  (`ip`|`input`) that this feature will build on; the input-specific
  background logic itself is still unscheduled (role values other than
  grouping are validated but no input-specific rules exist yet).

## 4. fastq_sort tool switch: ea-utils -> seqkit (resolved 2026-09-07)

- [x] `rule fastq_sort` now uses `seqkit sort -n` (seqkit 2.13.0 pinned in
  `workflow/environment.yaml`, ea-utils entry removed): the current bioconda
  ea-utils build (1.1.2.779) ships no `fastq-sort` binary (only fastq-clipper /
  fastq-join / fastq-mcf / fastq-multx / fastq-stats), so the pre-specified
  seqkit fallback was applied during the 2026-09-07 real-run validation. The
  runtime resolver and version-collector tool tables follow the rename
  (`fastq_sort`/`fastq-sort` -> `seqkit`); see CHANGELOG.md.
