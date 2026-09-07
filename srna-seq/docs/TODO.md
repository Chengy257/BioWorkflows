# srna-seq backlog (TODO)

> Items identified during v0.1 and deferred to a later release.

## 1. Differential expression — DONE in v0.2.0 (2026-09-08)

- Landed in v0.2.0 as the optional, default-off `deg` stage: DESeq2 over the per-class count matrices (`results/4.expression/{class}/{class}_counts.tsv` -> `results/6.DEG/{class}/`), driven by optional group/batch columns in the sample table (headers `sample_id[, group[, batch]]`) and the `deg:` config section (`workflow/rules/deg.smk` + `workflow/scripts/run_deseq2.R`, ported from the rna-seq reference). Enable with `deg.enabled: true`; see the user guide §6. Multi-factor designs beyond group/batch (e.g. interaction terms, time courses) remain future work.

## 2. Novel miRNA discovery

- miRDeep-P2-style novel miRNA prediction from the genome-aligned reads (hairpin excision from the reference genome, mapping-signature scoring); v0.1 only quantifies the mature sequences listed in the configured `miRNA` cascade fasta.

## 3. Unit-test suite

- The v0.1 regression is a synthetic-data dry-run (`tests/run_test.sh`). Add a unit-test suite (e.g. pytest) for the pure-Python pieces — `count_features.py`, `merge_counts.py`, `cascade_summary.py`, and the `common.smk` validation/species-merge logic — so refactors do not depend on the full snakemake dry-run. (Partially landed alongside v0.2.0: `tests/test_common.py` covers the sample-table loader, config validation, deg design checks, and the species-merge parity; the counting scripts are still uncovered.)
