# Example project: legacy human FBL-CLIP

This directory contains everything needed to start a seclip-seq project modeled on the legacy human FBL eCLIP dataset (five `FC` replicates).

## Contents

| File | Description |
|---|---|
| `samples.csv` | sample table for the legacy project's five replicates (`FC_rep1` .. `FC_rep5`), single `sample_id` column |
| `config.yaml` | filled copy of `config/config.template.yaml` for that project: `species: "hsa"`, `threads: 24`, legacy STAR/cutadapt/UMI defaults, `filter_repeats: true`, PureCLIP + CLIPper (`GRCh38_v40`); the three reference files still hold `/path/to/` placeholders with the legacy GRCh38.p14 / GENCODE v44 / sncRNA paths as comments |

For the miniature synthetic dataset (no real data needed) see `tests/make_testdata.py` (`bash tests/run_test.sh` assembles the working directory and dry-runs automatically); a real analysis project must supply its own FASTQ files.

## Adapt it to your project

1. Replace the sample ids in `samples.csv` with your own (naming rules in `docs/user-guide.md` §3) and point `SampleListFile` at it;
2. Replace the three `/path/to/` reference placeholders in `config.yaml` with real server paths (the commented lines show the legacy file set as an example);
3. Adjust `threads`, the `cutadapt`/`star` blocks, and the `callpeak` switches if your library design differs (UMI length, read length, species) — every key is documented in `docs/user-guide.md` §4.

## Start the example project in three commands

```bash
# 1) Create the working directory with the example table and config
mkdir -p ~/work/fbl
cp example/samples.csv example/config.yaml ~/work/fbl/

# 2) Place the raw data (SE fastq: <sample_id>.fastq.gz or <sample_id>_R1.fq.gz)
mkdir -p ~/work/fbl/1.rawdata
#   cp /data/raw/FC_rep1.fq.gz ~/work/fbl/1.rawdata/
#   ... (FC_rep2 .. FC_rep5)

# 3) Preflight -> dry-run -> run (execute from the repository's seclip-seq
#    directory; replace the reference placeholders in the config first)
bash run.sh -P ~/work/fbl --check-software
bash run.sh -P ~/work/fbl -n
bash run.sh -P ~/work/fbl -j 10
```

Notes:

- the config and sample table are auto-detected when they sit in the working directory as `config.yaml` / `samples.csv`; otherwise pass them explicitly (`-c FILE`, positional config);
- before launching, fill the reference paths — until then parse-time warnings point at the placeholders;
- cluster submission example: `bash run.sh -P ~/work/fbl --profile pbs --queue workq --memory 16G --runtime 600`; the full option list is in `run.sh --help`.
