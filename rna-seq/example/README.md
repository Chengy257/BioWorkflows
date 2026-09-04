# Example project

This directory provides all the templates needed to start a new analysis project.

## Directory contents

| File | Description |
|---|---|
| `config.yaml` | Example project configuration (copy into the project directory, then edit paths) |
| `samples_231107XTL.csv` | Real project sample table example (project 231107XTL, 27 samples / 8 treatment groups + control), usable as a sample table format reference |

A miniature test dataset (2 groups x 2 samples, downsampled reads) can be generated deterministically by `tests/make_testdata.py` (used by `make test`); example projects still need to supply their own real fastq files.

## Start a new project with one command

```bash
# 1. Create the project directory, copy the configuration and sample table
mkdir -p myproject && cd myproject
cp ../example/config.yaml .
cp ../example/samples_231107XTL.csv samples.csv   # or replace with your own sample table

# 2. Place the raw data (PE paired-end example)
mkdir -p 1.rawdata
#   cp /data/raw/{sample}_1.fastq.gz 1.rawdata/
#   cp /data/raw/{sample}_2.fastq.gz 1.rawdata/

# 3. Edit config.yaml: replace /path/to/... with real cluster paths, verify the sample table file name

# 4. Validate the sample table (optional but recommended)
python3 ../workflow/scripts/validate_samples.py samples.csv control

# 5. Dry-run preview
RNASEQ_PIPELINE=deg RNASEQ_CONFIG=$PWD/config.yaml \
  snakemake -s ../workflow/Snakefile --profile ../workflow/profile/default -n --quiet

# 6. Real run (SGE / local selected automatically)
bash ../run.sh deg . config.yaml 10
```

## Sample table format

```csv
id,group,layout
WT_1,control,auto
drugA_1,drugA,auto
```

- `id`: sample name; `-`, spaces, and `/` are forbidden; avoid names that are prefixes of each other;
- `group`: group name; the control group name must match `control_group` in the config;
- `layout`: `PE`/`SE`/`auto` (the current version detects from the actual files, so `auto` is fine).
