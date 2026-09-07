# srna-seq backlog (TODO)

> Items identified during v0.1 and deferred to a later release.

## 1. Differential expression — DONE in v0.2.0 (2026-09-08)

- Landed in v0.2.0 as the optional, default-off `deg` stage: DESeq2 over the per-class count matrices (`results/4.expression/{class}/{class}_counts.tsv` -> `results/6.DEG/{class}/`), driven by optional group/batch columns in the sample table (headers `sample_id[, group[, batch]]`) and the `deg:` config section (`workflow/rules/deg.smk` + `workflow/scripts/run_deseq2.R`, ported from the rna-seq reference). Enable with `deg.enabled: true`; see the user guide §6. Multi-factor designs beyond group/batch (e.g. interaction terms, time courses) remain future work.

## 2. Novel miRNA discovery — DONE in v0.2.0 (2026-09-08)

- Landed in v0.2.0 as the optional, default-off `novel_mirna` stage: miRDeep-P2 (bioconda `mirdeep-p2=1.1.4`, shipped as `miRDP2-v1.1.4_pipeline.bash`) run per sample over the trimmed reads (`workflow/rules/novel_mirna.smk`). Candidate hairpins are excised from the reference genome and scored by mapping signature; the tool-native result tree lands under `results/6.novel_mirna/{sample}/` (main table `{sample}_filter_P_prediction`). Enable with `novel_mirna.enabled: true` + `novel_mirna.mature_fasta` (a configured `genome.fasta` is mandatory); see the user guide §9. The rule wraps the tool's verified 1.1.4 surface — including its bundled mature-miRNA index (no user-facing mature flag), its required pre-built bowtie genome index (reused from `0.index/genome`), and the collapsed-fasta input format — with the full evidence recorded in the rule header comment. Known packaging limitation (documented, not worked around): the package does not ship the rfam ncRNA index the script references, so that filter silently no-ops.

## 3. Unit-test suite

- The v0.1 regression is a synthetic-data dry-run (`tests/run_test.sh`). Add a unit-test suite (e.g. pytest) for the pure-Python pieces — `count_features.py`, `merge_counts.py`, `cascade_summary.py`, and the `common.smk` validation/species-merge logic — so refactors do not depend on the full snakemake dry-run. (Partially landed alongside v0.2.0: `tests/test_common.py` covers the sample-table loader, config validation, deg design checks, and the species-merge parity; the counting scripts are still uncovered.)
