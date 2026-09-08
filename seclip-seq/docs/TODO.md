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

## 3. IP vs input control (resolved 2026-09-08, v0.2.0)

- [x] Implemented as `reproducible_peaks.input_control` (+ `filter_by_input`),
  both default off, building on the `role` column reserved in W7: role=input
  samples join the PureCLIP target set, one `input_background` job per
  condition unions their PureCLIP beds (`bedtools multiinter`, support >= 1 =
  plain union) into `results/6.reproducible_peaks/{condition}.input_background.bed`,
  and the ip consensus of conditions with inputs gains the binary
  `in_input_background` column (0/1) as column 5 of
  `{condition}.consensus.bed`; `filter_by_input: true` additionally writes
  `{condition}.consensus.filtered.bed` without the flagged sites (BED4) and
  routes the consensus annotation to it.
- Design rationale: the background is **condition-scoped** (an input control
  only pairs meaningfully with the ip replicates declared under the same
  `condition`, not across conditions); the union threshold is the lowest
  possible (support >= 1 — every interval ever seen in any input of the
  condition counts as background); and flagging is the default over filtering
  so both views remain available — the flagged BED5 records background
  membership per reproducible site, while the filtered BED4 is the opt-in
  cleaned product. Conditions without input samples keep the byte-identical
  W7 BED4 consensus (warning only), and a condition with inputs but no ip
  samples is a hard parse-time error (nothing to build or flag).
- Not pursued: PureCLIP-level background estimation from a control BAM and
  CLIPper contrast modes stay external options; IDR-style contrast ranking
  remains possible future work (see item 1).

## 4. fastq_sort tool switch: ea-utils -> seqkit (resolved 2026-09-07)

- [x] `rule fastq_sort` now uses `seqkit sort -n` (seqkit 2.13.0 pinned in
  `workflow/environment.yaml`, ea-utils entry removed): the current bioconda
  ea-utils build (1.1.2.779) ships no `fastq-sort` binary (only fastq-clipper /
  fastq-join / fastq-mcf / fastq-multx / fastq-stats), so the pre-specified
  seqkit fallback was applied during the 2026-09-07 real-run validation. The
  runtime resolver and version-collector tool tables follow the rename
  (`fastq_sort`/`fastq-sort` -> `seqkit`); see CHANGELOG.md.
