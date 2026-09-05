# rna-seq follow-up TODO

> Created during the enhancer_lncRNA_2026 real-data preparation stage (2026-09-05).

## 1. UMI library support — RESOLVED (2026-09-05)

- Trim Galore arguments are no longer hardcoded: the `trim:` config section
  (`quality` / `stringency` / `error_rate` / `extra`) mirrors the chip workflow
  and injects into both trim rules (`workflow/rules/align.smk`).
- UMI clipping is a first-class `umi:` config section (`enabled` / `read1_len` /
  `read2_len`) that appends per-sample `--clip5pNbases` to the STAR command
  (PE renders two comma-separated values, required by STAR 2.7.10b and
  data-verified on the server deployment; previously done via raw
  `star_extra_args`).
- `umi_tools extract` (fastq-level UMI extraction into read names) and
  `umi_tools dedup` (alignment-based deduplication): evaluated and deferred —
  umi_tools is absent from the environments/CI and adding it creates a solve
  risk, while clip-without-dedup is accepted by the project for now.
- Data basis: mRNA-ZH11-0H-1 first 500k reads per-position base composition —
  R1/R2 positions 1-8 random, position 9 = 99.5% T, genomic composition from
  position 10 (insert starts there); the caRNA control shows no such pattern.

## 2. FASTQ naming convention compatibility — RESOLVED (2026-09-05)

- `_raw_reads` (`workflow/rules/common.smk`) now recognizes the
  sequencer-delivery names `{id}_R1/_R2.fastq.gz` and `{id}_R1/_R2.fq.gz`
  in addition to the original `{id}_1/_2` patterns; the original patterns keep
  priority so existing projects resolve identically. The manual `{id}_1/_2`
  symlink workaround is no longer needed. All accepted names are documented in
  `docs/user-guide.md`.
- `validate_samples.py` validates the sample table only (it never inspects
  FASTQ file names), so no validator change was required.

## 3. PBS cluster support — resolved (2026-09-05)

- Status (resolved 2026-09-05): chip's `profile/pbs` was ported and `run.sh`'s
  profile whitelist and auto-detection updated (qsub without SGE_ROOT -> pbs);
  verified working on the server; to be committed together with the repository.
