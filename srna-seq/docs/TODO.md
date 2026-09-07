# srna-seq backlog (TODO)

> Items identified during v0.1 and deferred to a later release; none are scheduled yet.

## 1. Differential expression

- v0.1 produces per-class count/RPM matrices (`results/4.expression/`). Differential miRNA expression across conditions (DESeq2 on the combined count matrix; needs a group/batch design in the sample table) is out of scope for v0.1 and is the first v0.2 candidate.

## 2. Novel miRNA discovery

- miRDeep-P2-style novel miRNA prediction from the genome-aligned reads (hairpin excision from the reference genome, mapping-signature scoring); v0.1 only quantifies the mature sequences listed in the configured `miRNA` cascade fasta.

## 3. Unit-test suite

- The v0.1 regression is a synthetic-data dry-run (`tests/run_test.sh`). Add a unit-test suite (e.g. pytest) for the pure-Python pieces — `count_features.py`, `merge_counts.py`, `cascade_summary.py`, and the `common.smk` validation/species-merge logic — so refactors do not depend on the full snakemake dry-run.
