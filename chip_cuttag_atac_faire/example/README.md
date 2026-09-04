# Example project

This directory provides all templates needed to start a new analysis project.

## Contents

| File | Description |
|---|---|
| `samples.csv` | Sample table from a real project (myc / IgG, two samples, chip narrow, group `myc_vs_IgG`); use it as a format reference |
| `config.yaml` | Project-level override example, passable directly via `-c`; keys not listed keep the repository `config/config.yaml` defaults |

For the miniature test dataset (synthetic fastq) see `tests/make_testdata.py` (`bash tests/run_test.sh` assembles the working directory and dry-runs automatically); real analysis projects must supply their own fastq files.

## Start a new project in three commands

```bash
# 1) Create the working directory and copy the sample table (you can point -c at
#    example/config.yaml in the repository, or copy it and edit it)
mkdir -p ~/work/demo
cp example/samples.csv ~/work/demo/samples.csv

# 2) Place the raw data (PE paired-end: <sample_id>_1.fastq.gz / <sample_id>_2.fastq.gz;
#    common R1/R2 suffixes can be renamed in batch via bash run.sh -r)
mkdir -p ~/work/demo/1.rawdata
#   cp /data/raw/myc_R1.fastq.gz ~/work/demo/1.rawdata/myc_1.fastq.gz
#   cp /data/raw/myc_R2.fastq.gz ~/work/demo/1.rawdata/myc_2.fastq.gz

# 3) Dry-run preview -> real run (execute from the repository root; replace the
#    reference genome paths in the config first)
bash run.sh -P ~/work/demo -n
bash run.sh -P ~/work/demo -c example/config.yaml --profile auto -j 10
```

Notes:

- `-c example/config.yaml` names the project config explicitly; if you copy `example/config.yaml` into the working directory as `config.yaml`, `-c` can be omitted (auto-detected). `config.grouplist` resolves relative to the **working directory**, so the sample table must live in the working directory (or use an absolute path).
- Before launching, run `bash run.sh -P ~/work/demo --check-software` to preflight the tools and the R environment.
- Cluster submission example: `--profile pbs --queue workq --memory 16G --runtime 600`; for the full option list see `run.sh --help`.

## Sample table format

For the 6-column schema template and row-level validation rules see [config/samples.csv](../config/samples.csv):

```csv
sample_id,role,group,seqtype,layout,peak_type
myc,treat,myc_vs_IgG,chip,PE,narrow
IgG,control,myc_vs_IgG,chip,PE,narrow
```

- `sample_id`: sample name; `-`, spaces, `/`, and double underscores are forbidden, and samples must not be prefixes of each other;
- `role`: `treat` / `control` (groups without a control may list treat rows only);
- `group`: group name; `seqtype`, `peak_type` must be consistent within each group;
- `seqtype`: `chip` / `cuttag` / `atac` / `faire`; mixed-assay projects are routed row by row automatically;
- `layout`: only `PE` (paired-end) is supported in this version;
- `peak_type`: required for chip only (`narrow` / `broad`); atac/faire use `none`.
