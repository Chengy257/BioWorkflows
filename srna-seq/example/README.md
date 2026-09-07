# Example project: legacy rice small-RNA

This directory contains everything needed to start an srna-seq project modeled on the legacy rice sRNA dataset (root and leaf, two replicates each).

## Contents

| File | Description |
|---|---|
| `samples.csv` | sample table for the example project (`root_rep1`, `root_rep2`, `leaf_rep1`, `leaf_rep2`) with `sample_id,group` columns (the group feeds the optional DE stage; a single `sample_id` column remains valid) |
| `config.yaml` | filled copy of `config/config.template.yaml` for that project: `species: "osa"` (rice IRGSP-1.0 preset), `threads: 24` (legacy PBS cap), the legacy Trim Galore parameters, the legacy cascade class order (rRNA → snoRNA → snRNA → tRNA → miRNA → mRNA → rhizo), and the genome alignment of the trimmed reads; the class FASTAs and the genome still hold `/path/to/` placeholders, with the legacy project's `0.index/` references from the old `srna-seq.sh` script as comments |

For the miniature synthetic dataset (no real data needed) see `tests/make_testdata.py` (`bash tests/run_test.sh` assembles the working directory and dry-runs automatically); a real analysis project must supply its own FASTQ files.

## Adapt it to your project

1. Replace the sample ids in `samples.csv` with your own (naming rules in `docs/user-guide.md` §3) and point `SampleListFile` at it;
2. Replace the `/path/to/` reference placeholders in `config.yaml` with real server paths (the commented lines show where the legacy project kept its pre-built bowtie indices; this workflow builds the indices itself from the FASTAs);
3. Adjust `threads`, the `trim` block, and the `cascade` list if your design differs (drop a class by removing its entry, add one anywhere in the order) — every key is documented in `docs/user-guide.md` §4.

## Start the example project in three commands

```bash
# 1) Create the working directory with the example table and config
mkdir -p ~/work/rice
cp example/samples.csv example/config.yaml ~/work/rice/

# 2) Place the raw data (SE fastq: <sample_id>.fastq.gz or <sample_id>_R1.fq.gz)
mkdir -p ~/work/rice/1.rawdata
#   cp /data/raw/root_rep1.fq.gz ~/work/rice/1.rawdata/
#   ... (root_rep2, leaf_rep1, leaf_rep2)

# 3) Preflight -> dry-run -> run (execute from the repository's srna-seq
#    directory; replace the reference placeholders in the config first)
bash run.sh -P ~/work/rice --check-software
bash run.sh -P ~/work/rice -n
bash run.sh -P ~/work/rice -j 10
```

Notes:

- the config and sample table are auto-detected when they sit in the working directory as `config.yaml` / `samples.csv`; otherwise pass them explicitly (`-c FILE`, positional config);
- before launching, fill the reference paths — until then parse-time warnings point at the placeholders;
- cluster submission example: `bash run.sh -P ~/work/rice --profile pbs --queue workq --memory 16G --runtime 600`; the full option list is in `run.sh --help`.
