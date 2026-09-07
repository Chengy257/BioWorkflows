# Example project: two-sample bisulfite sequencing

This directory contains everything needed to start a bs-seq project (Bismark methylation calling) with two samples.

## Contents

| File | Description |
|---|---|
| `samples.csv` | sample table for the example project (`s1`, `s2`), single `sample_id` column — replace them with your own ids |
| `config.yaml` | filled copy of `config/config.template.yaml` with the full bs-seq schema: `species: "osa"` (rice IRGSP-1.0 preset), the Trim Galore parameters, the Bismark alignment extras, and the methylation-extraction options (`cx_report` / `merge_cpg` / `buffer_frac`); the genome still holds a `/path/to/` placeholder |

For the miniature synthetic dataset (no real data needed) see `tests/make_testdata.py` (`bash tests/run_test.sh` assembles the working directory and dry-runs automatically); a real analysis project must supply its own FASTQ files.

## Adapt it to your project

1. Replace the sample ids in `samples.csv` with your own (naming rules in `docs/user-guide.md` §3) and point `SampleListFile` at it;
2. Replace the `/path/to/` genome placeholder in `config.yaml` with a real server path (an explicit non-empty `genome` wins over the species preset; with `species: "none"` the key is mandatory); the FASTQ phred encoding must match `bismark.align_extra`;
3. Adjust `threads`, the `trim` block, and the `methylation_extractor` block if your design differs (`trim.enabled: false` restores the legacy raw-read alignment; every key is documented in `docs/user-guide.md` §4). Paired-end and single-end samples can coexist — the layout is auto-detected per sample from the `1.rawdata/` file names.

## Start the example project in three commands

```bash
# 1) Create the working directory with the example table and config
mkdir -p ~/work/bsseq
cp example/samples.csv example/config.yaml ~/work/bsseq/

# 2) Place the raw data (PE: <sample_id>_1.fastq.gz + <sample_id>_2.fastq.gz;
#    SE: <sample_id>.fastq.gz)
mkdir -p ~/work/bsseq/1.rawdata
#   cp /data/raw/s1_1.fastq.gz ~/work/bsseq/1.rawdata/
#   ... (s1_2, s2_1, s2_2)

# 3) Preflight -> dry-run -> run (execute from the repository's bs-seq
#    directory; replace the genome placeholder in the config first)
bash run.sh -P ~/work/bsseq --check-software
bash run.sh -P ~/work/bsseq -n
bash run.sh -P ~/work/bsseq -j 10
```

Notes:

- the config and sample table are auto-detected when they sit in the working directory as `config.yaml` / `samples.csv`; otherwise pass them explicitly (`-c FILE`, positional config);
- before launching, fill the genome path — until then parse-time warnings point at the placeholder;
- cluster submission example: `bash run.sh -P ~/work/bsseq --profile pbs --queue workq --memory 32G --runtime 720`; the full option list is in `run.sh --help`.
