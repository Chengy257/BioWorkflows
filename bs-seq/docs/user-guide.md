# bs-seq workflow user guide

> Updated: 2026-09-05 (v0.1.0)
> Intended audience: analysts running this workflow on their cluster/server
> The workflow never creates Conda environments on its own: the runtime environment is created explicitly from `workflow/environment.yaml`, or reuses the server's existing environment via `config/software.yaml`. Cluster scheduling resources are declared per rule in `config/resources.yaml` and support project-level overrides. Raw inputs live in `1.rawdata/` at the working-directory root; every derived output is consolidated under the project's `results/` directory (`results_dir` in config).

---

## 1. Installation and environment setup

### 1.1 Prerequisites

| Dependency | Notes |
|---|---|
| Linux (PBS/SGE/SLURM cluster or a single machine; WSL works) | the scheduler is auto-detected by `run.sh` or set via `--profile` |
| Snakemake | reference version **7.32.4** (pinned in `workflow/environment.yaml`); see the version matrix in §1.4 |
| Python 3 + PyYAML | the bootstrap dependency `run.sh` needs to resolve `software.yaml` and run preflight checks; must be on the main PATH |
| Analysis tools | Bismark 0.24.x suite (bismark / bismark_genome_preparation / deduplicate_bismark / bismark_methylation_extractor / coverage2cytosine / bam2nuc / bismark2report / bismark2summary) / bowtie2 / samtools / Trim Galore (brings Cutadapt + FastQC) / MultiQC |

Launcher bootstrap order: `run.sh` first uses the `python3` on the system PATH (with PyYAML) to resolve `config/software.yaml`, then injects the main environment/tools into the current process — so even on the conda_prefix reuse route, the login-node PATH must have `python3` (snakemake is supplied by the main environment after resolution).

### 1.2 Three ways to set up the environment

**Path 1: fresh server — create the environment from the all-in-one template**

```bash
mamba env create -f workflow/environment.yaml   # environment name bs-seq
conda activate bs-seq
```

This command is executed explicitly by the user; Snakemake never creates or modifies software environments on its own. Pinned versions in the template: snakemake-minimal 7.32.4, bismark 0.24.0, bowtie2 2.5.2, samtools 1.17, trim-galore 0.6.10, fastqc 0.11.9, multiqc 1.21.

**Path 2: reuse an existing conda environment on the server**

Copy `config/software.yaml` into the project directory (or edit the repository default) and set:

```yaml
environment:
  type: conda
  conda_prefix: "/share/conda/envs/bs-seq"   # recommended on HPC
  # conda_name: "bs-seq"                     # or resolve by environment name
  strict: true
```

No `conda activate` needed: after resolution, `run.sh` injects that prefix's `bin` into PATH. An environment created via path 1 can also be connected this way (`environment.yaml` then only serves to "create the environment"; the runtime always goes through `software.yaml`).

**Path 3: system mode**

When all tools already come from the system PATH (the repository default):

```yaml
environment:
  type: system
  strict: true
```

`strict: true` makes the preflight fail loudly when a required executable is missing instead of continuing silently. All Bismark suite executables (bismark / bismark_genome_preparation / deduplicate_bismark / bismark_methylation_extractor / coverage2cytosine / bam2nuc / bismark2report / bismark2summary), bowtie2, samtools, trim_galore, fastqc, multiqc, and python3 resolve from PATH automatically; `software.yaml` needs no `tools:` overrides unless a binary lives outside PATH. Like the other v0.1 workflows there is no R section.

### 1.3 Preflight checks

```bash
# Check that all required executables resolve (the Bismark suite, bowtie2,
# samtools, trim_galore, fastqc, multiqc, python3); then exit
bash run.sh -P /path/to/workdir --check-software

# Parse the workflow + all configs (sample-table/config validation runs at parse time),
# list the rules, then exit without starting jobs
bash run.sh -P /path/to/workdir -c /path/to/workdir/config.yaml --validate-only
```

Notes:

- `--check-software` requires `-P` to name the project directory (created automatically when missing);
- before a real run, `run.sh` executes a full preflight **automatically**; `--skip-software-check` skips it, and `-n` (dry-run) skips it automatically;
- `software.yaml` lookup order: `--software FILE` > `<project dir>/software.yaml` > repository `config/software.yaml`.

### 1.4 Snakemake version matrix

| snakemake version | support | notes |
|---|---|---|
| **7.32.4** | reference version | pinned in `workflow/environment.yaml`; cluster profiles (pbs/sge/slurm) target the 7.x classic `--cluster` interface |
| **8.x** | cluster semantics unverified | parsing/dry-run work; 8.x replaced the classic `--cluster` interface with the executor plugin system — test on a small scale before real cluster runs |
| <7 or >=9 | unverified | test before use |

---

## 2. Quick start

Six steps to a minimal project (the two-sample example from [example/](../example/README.md)):

**Step 1: create the working directory**

```bash
mkdir -p ~/work/bsseq/1.rawdata
```

The working directory (data) and the analysis code directory (this repository) are separate; `run.sh` generates all results inside the working directory.

**Step 2: add the raw data**

FASTQ files named `{sample}_1.fastq.gz` + `{sample}_2.fastq.gz` (paired-end), `{sample}_1.fq.gz` + `{sample}_2.fq.gz` (paired-end), or `{sample}.fastq.gz` / `{sample}.fq.gz` (single-end) — exact names, checked when the workflow resolves the sample:

```bash
cp s1_1.fastq.gz s1_2.fastq.gz s2_1.fastq.gz s2_2.fastq.gz ~/work/bsseq/1.rawdata/
```

**Step 3: write the sample table**

```bash
cp example/samples.csv ~/work/bsseq/samples.csv
```

Edit the file with your own samples (a single `sample_id` column; validation rules in §3).

**Step 4: write the project config**

```bash
cp example/config.yaml ~/work/bsseq/config.yaml
```

Fill in the genome FASTA path (or keep `species: "osa"` / `"hsa"` and let the `config/species.yaml` preset supply it — an unset `genome` inherits the preset; see §4.2/§4.3). `run.sh` auto-detects `<project dir>/config.yaml` when `-c` is omitted, so placing it in the working directory is enough. Full key reference in §4.

**Step 5: preflight and dry-run**

```bash
cd /path/to/repo/bs-seq
bash run.sh -P ~/work/bsseq --check-software   # environment preflight
bash run.sh -P ~/work/bsseq -n                 # dry-run: builds the DAG and prints planned jobs without executing
```

**Step 6: launch**

```bash
bash run.sh -P ~/work/bsseq -j 10              # local / auto-detected scheduler
# positional form:
bash run.sh ~/work/bsseq ~/work/bsseq/config.yaml 10
# cluster example:
bash run.sh -P ~/work/bsseq --profile pbs --queue workq --memory 32G --runtime 720
```

After a successful launch the main log is `~/work/bsseq/snakemake.logs.txt` (change with `--log FILE`), and per-rule logs live under `results/logs/`. See §7 for the results layout.

---

## 3. Sample table in detail

The sample table is a **single-column CSV** (template `config/samples.csv`):

```csv
sample_id
s1
s2
```

Every row is one sample; v0.1 has no other columns (no group/condition/batch design yet — see the [TODO backlog](TODO.md)). Each `sample_id` must match the prefix of FASTQ file(s) in `1.rawdata/`.

### 3.1 Validation rules (at parse time, errors carry line numbers)

Validation runs centrally at Snakemake parse time (`load_sample_table` in `workflow/rules/common.smk`); any violation aborts the workflow before any job runs:

1. The header must be **exactly one column named `sample_id`** — otherwise: `Sample table <path> must have exactly one column with header 'sample_id' (got [...]); see config/samples.csv`;
2. `sample_id` must not be empty — `Sample table line <n>: sample_id must not be empty`;
3. `sample_id` allows only **alphanumerics plus `. _ -`**, must start alphanumeric, must not start with `-`, and must not contain consecutive underscores `__` (ids become file names and Bismark command-line arguments) — `Sample table line <n>: sample_id='...' contains illegal characters; only alphanumerics and . _ - are allowed (no leading '-', no '__')`;
4. `sample_id` must be unique — `Sample table line <n>: duplicate sample_id '...'`;
5. The table must have at least one data row — `Sample table <path> has no data rows`.

### 3.2 Raw FASTQ resolution and library layout

For each declared sample the workflow detects the library layout by looking for exact file names under `1.rawdata/`, **paired-end first**: `{sample}_1.fastq.gz` + `{sample}_2.fastq.gz`, then `{sample}_1.fq.gz` + `{sample}_2.fq.gz`; if no PE pair exists, single-end `{sample}.fastq.gz` / `{sample}.fq.gz`. If nothing matches, the workflow fails with the supported-name list (`sample <s>: no raw FASTQ files were found under 1.rawdata/ ...`).

The detected layout drives the whole per-sample chain: the PE or SE trim definition, the `--pe`/`-p` flags of `bismark`/`deduplicate_bismark`/`bismark_methylation_extractor`, and the `--1/--2` vs positional FASTQ arguments. PE and SE samples can coexist in one run.

---

## 4. Configuration in detail

### 4.1 Config chain and layering

The `workflow/Snakefile` assembles the configuration in this order; **later files override earlier keys**:

1. Repository `config/species.yaml` (species preset map; loaded first);
2. Repository `config/resources.yaml` (per-rule scheduler resources; the project copy is injected by run.sh instead when present);
3. Project config: `-c FILE` explicitly; when omitted, `<project dir>/config.yaml` is auto-detected; when neither exists, the repository `config/config.yaml` defaults apply;
4. Overlay config: `-l FILE` explicitly; when omitted, `<project dir>/config.local.yaml` is auto-layered last (an INFO line is printed when found).

Typical usage: the project config holds project-specific values (genome path, samples); machine-specific differences (e.g. a lower thread cap) go into `config.local.yaml` (kept out of version control). `config/config.template.yaml` is the annotated template — copy it into the working directory and change only what differs (unlisted keys keep the repository defaults).

Scheduler resources are layered separately: see §4.4 for `config/resources.yaml` and the project-level `resources.yaml`.

### 4.2 All config.yaml keys

| Key | Type | Default | Effect |
|---|---|---|---|
| `results_dir` | string | `"results"` | root directory for all derived outputs, relative to the working directory (raw inputs `1.rawdata/` stay at the root) |
| `SampleListFile` | string | `"config/samples.csv"` | sample-table path; resolves absolute > working-directory-relative > repository-relative (§3) |
| `threads` | int >= 1 | `12` | legacy global thread cap: every rule's thread request (§4.4) is capped at this value |
| `species` | string | `"osa"` | selects a reference preset from `config/species.yaml`; supported values are `osa` (rice IRGSP-1.0), `hsa` (human GRCh38.p14), and `none` (no preset, no fallback — anything else errors at parse time) |
| `genome` | non-empty string | from the species preset | plain genome FASTA for `bismark_genome_preparation`; **unset = inherit the species preset** (`/path/to/` placeholders only warn); under `species: "none"` you must set it explicitly or validation fails |
| `trim.enabled` | bool | `true` | `false` = align the **raw** reads directly (legacy behavior); the trimmed outputs are then never requested and no trim job runs |
| `trim.quality` | int >= 0 | `20` | 3' quality-trimming cutoff passed to Trim Galore (`-q`) |
| `trim.min_len` | int >= 1 | `20` | discard trimmed reads shorter than this (`--length`) |
| `trim.adapter` | string of `ACGTN` | `"AGATCGGAAGAGC"` | adapter sequence (`-a`, and `-A` in paired-end mode) |
| `trim.stringency` | int >= 1 | `3` | adapter overlap required to trigger trimming (`--stringency`) |
| `trim.error_rate` | float in (0, 1] | `0.1` | fraction of mismatches allowed in an adapter match (`-e`) |
| `trim.extra` | string | `""` | extra Trim Galore arguments appended verbatim |
| `bismark.align_extra` | string | `"--phred33-quals"` | extra bismark arguments appended verbatim (keep the phred encoding in sync with your data) |
| `methylation_extractor.cx_report` | bool | `false` | `true` = also run the extractor with `--CX --cytosine_report` (the full all-context cytosine report; large) |
| `methylation_extractor.merge_cpg` | bool | `true` | request the per-sample `coverage2cytosine --merge_CpG` step (`{sample}.CpG_merged.tsv.gz`) |
| `methylation_extractor.buffer_frac` | int >= 1 | `4` | the extractor `--buffer_size` (GB) = the rule's `mem_mb / 1024 / buffer_frac`, minimum 1 (§5.2) |
| `resources` | mapping of rule -> {threads, mem_mb, runtime_min} | `{}` | per-rule scheduler-resource overrides, same shape as `config/resources.yaml` (§4.4) |

Validation runs at parse time (`validate_config` in `workflow/rules/common.smk`) and aggregates all problems into a single `WorkflowError` report (missing keys, type/range checks for every key above). Reference files that are still `/path/to/` placeholders only **print a warning and do not abort** — dry-run/lint often run on machines without the reference files; in a real run a missing genome fails `bismark_genome_prep`, so confirm each warning before launching.

### 4.3 Species presets

`config/species.yaml` currently defines two presets:

```yaml
osa:
  label: "Rice Oryza sativa (IRGSP-1.0)"
  genome: "/path/to/reference/osa/genome.fa"

hsa:
  label: "Human Homo sapiens (GRCh38.p14)"
  genome: "/path/to/reference/hsa/GRCh38.p14.genome.fa"
```

Preset values fill **only the keys the project config does not set**: an unset `genome` inherits `preset.genome` after choosing the species; an explicit non-empty `genome` in the project config wins. `species: "none"` skips the preset merge entirely — nothing is inherited, so `genome` must then be set explicitly. The shipped paths are placeholders — replace them with real cluster paths at deployment time. (The `label` key is descriptive metadata only.)

### 4.4 Per-rule scheduler resources (`config/resources.yaml`)

Scheduler resources live in their own file — the single place to tune cluster requests. Each rule declares `threads`, `mem_mb`, and `runtime_min`; `runtime_sec` (used by PBS walltime and SGE h_rt) is derived automatically as `runtime_min x 60` and cannot be set independently. The legacy top-level `threads` in config.yaml still acts as a global cap on every rule's threads.

Override without touching the repository: copy the file into your working directory as `resources.yaml` and edit — `run.sh` auto-detects it. Lookup order: the `BSSEQ_RESOURCES_CONFIG` environment variable > `<project dir>/resources.yaml` > repository `config/resources.yaml`. Rules not listed in an override keep the repository defaults.

Example override (`<project dir>/resources.yaml`):

```yaml
resources:
  bismark_align:
    mem_mb: 64000
    runtime_min: 1440
  bismark_genome_prep:
    threads: 16
```

Built-in defaults per rule (identical to `config/resources.yaml`):

| Rule | threads | mem_mb | runtime_min |
|---|---:|---:|---:|
| `software_versions` | 1 | 1024 | 10 |
| `trim` | 4 | 4096 | 120 |
| `fastqc` | 2 | 2048 | 30 |
| `multiqc` | 2 | 4096 | 30 |
| `bismark_genome_prep` | 12 | 16000 | 240 |
| `bismark_align` | 12 | 32000 | 720 |
| `deduplicate` | 4 | 16000 | 240 |
| `bam2nuc_genome` | 2 | 8000 | 60 |
| `bam2nuc_sample` | 2 | 8000 | 60 |
| `methylation_extractor` | 4 | 32000 | 1440 |
| `coverage2cytosine` | 4 | 16000 | 720 |
| `bismark2report` | 1 | 4096 | 30 |
| `bismark2summary` | 1 | 4096 | 30 |

(`fastqc` has its own entry for completeness; in v0.1 FastQC runs inside the trim jobs via `trim_galore --fastqc`. The `mem_mb` of `methylation_extractor` also feeds the `--buffer_size` derivation, §5.2.)

### 4.5 software.yaml

`config/software.yaml` (or a project copy) configures the runtime, separate from analysis parameters:

```yaml
environment:
  type: system                  # system | conda
  conda_prefix: ""              # preferred for HPC, e.g. /share/conda/envs/bs-seq
  conda_name: ""                # alternative to conda_prefix
  strict: true

tools: {}          # add overrides only for binaries outside the main PATH
paths: {}
databases: {}
```

All analysis tools (the Bismark suite, bowtie2, samtools, trim_galore, fastqc, multiqc, python3) are ordinary PATH executables resolved through `workflow/scripts/runtime_config.py`; `run.sh` exports the resolved values to the workflow. Nothing here needs a hand-installed absolute path unless your site installs a tool outside PATH.

---

## 5. Run modes

### 5.1 What runs — trimming, layouts, and optional steps

Per-sample decisions, all made at DAG-building time:

| Setting / detection | what enters the DAG |
|---|---|
| PE sample detected | `trim_pe` → PE bismark alignment (`--pe --1/--2`) → PE dedup (`--pe`) → PE extraction (`-p`) |
| SE sample detected | `trim_se` → SE bismark alignment (positional FASTQ) → SE dedup → SE extraction |
| `trim.enabled: false` | no trim job at all; the raw reads in `1.rawdata/` are aligned directly (legacy behavior) |
| `methylation_extractor.merge_cpg: true` | one `coverage2cytosine` job per sample |
| `methylation_extractor.cx_report: true` | the extractor additionally writes the full CX-context cytosine report |
| always | `bismark_genome_prep` + `bam2nuc_genome` (once per run), `bam2nuc_sample` / `bismark2report` per sample, `bismark2summary` + `multiqc` + `software_versions` once |

The generated cytosine reports always cover the standard contexts (CpG/CHG/CHH) plus the gzipped bedGraph; the extra per-sample report inputs (splitting report, M-bias) require the genome, which every run has by construction.

### 5.2 The extractor `--buffer_size` is derived from the rule memory

`bismark_methylation_extractor` loads a genomic buffer into RAM; a buffer larger than the job's actual memory is the classic OOM-on-the-cluster cause. The workflow therefore derives it from the scheduler request instead of hard-coding it:

```
--buffer_size <G>   with  G = max(1, resources.methylation_extractor.mem_mb / 1024 / methylation_extractor.buffer_frac)
```

With the shipped defaults (32000 MB, `buffer_frac: 4`) that is `7G`. Lower `buffer_frac` (or a larger `mem_mb` override in a project `resources.yaml`) gives a bigger buffer; raise `buffer_frac` when jobs are killed for exceeding memory.

### 5.3 Checks, dry-run, and direct Snakemake invocation

```bash
bash run.sh -P . -n                    # dry-run: builds the DAG and prints the jobs it would run, without executing
bash run.sh -P . --validate-only       # parses workflow + config, lists rules, then exits
bash run.sh -P . --check-software      # executable preflight
```

The dry-run automatically skips the software preflight (needs only snakemake + python3/PyYAML, not the full analysis environment). Sample-table and config validation happen at parse time, so the dry-run and `--validate-only` catch the same errors listed in §3.1 / §4.2.

Direct Snakemake invocation also works without the launcher:

```bash
BSSEQ_CONFIG=/path/to/config.yaml snakemake -s workflow/Snakefile -j 10
```

with the env vars the launcher would normally provide (`BSSEQ_RESOURCES_CONFIG`, `BSSEQ_EXTRA_CONFIG`).

### 5.4 Resuming and Snakemake passthrough

- Resuming: Snakemake skips completed steps based on output timestamps; after an interruption, **simply rerun the same command**;
- profiles pin `keep-going: true` and `rerun-incomplete: true`; `latency-wait` defaults to 90 (pbs) or 60 (default/sge/slurm), overridable via `--latency-wait SEC`;
- a stale working-directory lock (e.g. a job was force-killed) is cleared with `bash run.sh -P . --unlock`;
- arguments after `--` pass through to Snakemake untouched, e.g. `bash run.sh -P . -- --rerun-triggers mtime`;
- the launcher log defaults to `snakemake.logs.txt` (change with `--log FILE`); `-q` reduces the launcher's own output.

---

## 6. Cluster submission

### 6.1 The four profiles

Profiles pin cluster parameters in the repository (`workflow/profile/<name>/config.yaml`, see `workflow/profile/README.md`):

| profile | scheduler | behavior |
|---|---|---|
| `default` | local / standalone server | no cluster command, executes locally (8 cores pinned) |
| `pbs` | PBS (Torque), snakemake 7.x classic interface | `qsub -V -N {rule} -l select=1:ncpus={threads}:mem=<mem> -l walltime=<sec> -j oe`; walltime in seconds avoids `[HH:]MM:SS` ambiguity |
| `sge` | SGE, snakemake 7.x classic interface | `qsub -V -N {rule} -l ncpus={threads} -l h_vmem=<mem> -l h_rt=<sec> -j oe` |
| `slurm` | SLURM, snakemake 7.x classic interface | `sbatch --parsable -J {rule} -c {threads} --mem=<mem> --time=<min> -o slurm-{rule}-%j.out` |

All four profiles pin `keep-going: true`, `rerun-incomplete: true`, `printshellcmds: true`. `--profile local` is equivalent to `default`.

### 6.2 Auto detection

`--profile auto` (the default; presettable via the `RUN_PROFILE` environment variable) probes in this order:

1. `sbatch` on PATH -> **slurm**;
2. `qsub` present with `SGE_ROOT` set -> **sge**; `qsub` present without `SGE_ROOT` -> **pbs**;
3. neither -> **default** (local).

An explicitly requested profile whose scheduler command is missing (e.g. pbs without qsub) errors out immediately.

### 6.3 The resource model

The resources requested from the cluster come from three layers; later layers win:

1. **Repository defaults**: `config/resources.yaml` (the table in §4.4);
2. **Project override**: a `resources.yaml` copied into the project directory (auto-detected by run.sh; or point `BSSEQ_RESOURCES_CONFIG` at any file);
3. **Global forced overrides**: `--memory VALUE` (e.g. `16G`) overrides every rule's memory request; `--runtime MIN` (minutes) overrides every rule's walltime.

The global overrides are only valid with cluster profiles (passing them with the local default profile errors out). When nothing is overridden, the cluster command uses the placeholders `{resources.mem_mb}` / `{resources.runtime_sec}` / `{resources.runtime_min}`, filled at runtime by snakemake from each rule's declaration.

### 6.4 Cluster options and examples

| Option | Profiles | Notes |
|---|---|---|
| `--queue NAME` | pbs / sge | queue name (passing it with slurm errors out) |
| `--partition NAME` | slurm | partition name |
| `--memory VALUE` | pbs / sge / slurm | global memory override (`8000`, `8G`, `16GB`, `32GiB` all accepted) |
| `--runtime MIN` | pbs / sge / slurm | global walltime override (positive integer minutes) |
| `--sge-mem-resource N` | sge | SGE memory resource name, default `h_vmem` |
| `--scheduler-extra S` | cluster profiles | raw text appended to the qsub/sbatch submit command |
| `--retries N` | all | automatic retry count for failed jobs (default 0) |
| `--latency-wait SEC` | all | wait for output visibility on shared filesystems |
| `--max-jobs-per-sec N` / `--max-status-per-sec N` | all | throttle submission/status-check rates (useful under tight cluster quotas) |
| `--unlock` | all | clear the Snakemake working-directory lock |

Examples (copy as needed):

```bash
# PBS: workq queue, all rules forced to 32G/720min, retry failed jobs once
bash run.sh -P . --profile pbs --queue workq --memory 32G --runtime 720 --retries 1

# SGE: all.q queue, submit with the per-rule resources
bash run.sh -P . --profile sge --queue all.q --runtime 720 --retries 2

# SLURM: compute partition, 40 concurrent jobs
bash run.sh -P . --profile slurm --partition compute --memory 32G -j 40
```

The concurrency cap is `-j/--jobs N` (default 10; the pbs/sge/slurm profiles pin `jobs: 20`, overridden by `-j`).

### 6.5 Environment variable presets

Common options can be preset via environment variables (explicit command-line arguments win):

| Variable | Equivalent option |
|---|---|
| `RUN_PROFILE` | `--profile` |
| `BSSEQ_JOBS` | `-j / --jobs` (default 10) |
| `BSSEQ_LOG` | `--log` (default `snakemake.logs.txt`) |
| `BSSEQ_QUEUE` | `--queue` |
| `BSSEQ_PARTITION` | `--partition` |
| `BSSEQ_MEMORY` | `--memory` |
| `BSSEQ_RUNTIME_MIN` | `--runtime` |
| `BSSEQ_SGE_MEMORY_RESOURCE` | `--sge-mem-resource` (default `h_vmem`) |
| `BSSEQ_SCHEDULER_EXTRA` | `--scheduler-extra` |
| `BSSEQ_RETRIES` | `--retries` |
| `BSSEQ_LATENCY_WAIT` | `--latency-wait` |
| `BSSEQ_MAX_JOBS_PER_SECOND` | `--max-jobs-per-sec` |
| `BSSEQ_MAX_STATUS_CHECKS_PER_SECOND` | `--max-status-per-sec` |
| `BSSEQ_SOFTWARE_CONFIG` | `--software` |
| `BSSEQ_CONFIG` | `-c / --config` |
| `BSSEQ_EXTRA_CONFIG` | `-l / --extra-config` |
| `BSSEQ_RESOURCES_CONFIG` | project `resources.yaml` lookup |

One further variable is internal but also read by direct Snakemake invocation: `BSSEQ_PYTHON` (interpreter used by the `software_versions` rule).

---

## 7. Output interpretation

### 7.1 Working-directory layout

```
workdir/
├── config.local.yaml        # optional: working-directory config overlay (auto-layered, not committed)
├── samples.csv              # sample table (filename free; pointed to by SampleListFile)
├── snakemake.logs.txt       # launcher/Snakemake main log
├── 1.rawdata/               # raw inputs (the only data directory outside results/):
│                            #   {sample}_1.fastq.gz + {sample}_2.fastq.gz (PE) or {sample}.fastq.gz (SE)
└── results/                 # ALL derived artifacts (rename via results_dir in config)
    ├── 0.index/
    │   └── bismark_genome/    # bisulfite genome (Bisulfite_Genome/{CT,GA}_conversion/) + genomic_nucleotide_totals.txt
    ├── 2.cleandata/
    │   ├── {sample}_trimmed.fq.gz / {sample}_{1_val_1,2_val_2}.fq.gz   # alignment inputs
    │   ├── {sample}_trimming_report.txt (+ per-mate for PE)            # Trim Galore statistics
    │   └── fastqc/             # {sample}_trimmed_fastqc.html / .zip (PE: per mate)
    ├── 3.align/
    │   └── {sample}.bam + {sample}_report.txt                          # Bismark alignment + report
    ├── 4.dedup/
    │   └── {sample}.bam.deduplicated.bam + {sample}.dedup_report.txt   # deduplicated alignment
    ├── 5.methylation/
    │   └── {sample}/           # cytosine report, bedGraph, M-bias, splitting report,
    │                           # nucleotide stats, CpG_merged table, {sample}_seq_context.html
    ├── 5.QC/
    │   ├── bismark2summary.html
    │   ├── multiqc/multiqc_report.html
    │   ├── software_versions.yaml
    │   └── logs/               # QC rule logs (e.g. software_versions.log.txt)
    └── logs/                   # per-rule logs (trim_pe/, bismark_align/, deduplicate/, ...)
```

### 7.2 Results quick reference

| Result | Path | Meaning |
|---|---|---|
| bisulfite genome index | `results/0.index/bismark_genome/Bisulfite_Genome/` | C→T and G→A converted bowtie2 indexes built by `bismark_genome_preparation` |
| genome nucleotide totals | `results/0.index/bismark_genome/genomic_nucleotide_totals.txt` | genome-wide nucleotide composition (bam2nuc; input to the per-sample stats) |
| trimmed reads | `results/2.cleandata/{sample}_trimmed.fq.gz` (SE) or `{sample}_1_val_1.fq.gz` + `{sample}_2_val_2.fq.gz` (PE) | quality- + adapter-trimmed, `min_len`-filtered; the alignment input when trimming is enabled |
| trimming report | `results/2.cleandata/{sample}_trimming_report.txt` | Trim Galore statistics (aggregated into MultiQC) |
| FastQC | `results/2.cleandata/fastqc/{sample}_trimmed_fastqc.html` | per-base quality of the trimmed reads |
| alignment BAM | `results/3.align/{sample}.bam` | Bismark alignment of the trim-aware reads |
| alignment report | `results/3.align/{sample}_report.txt` | alignment efficiency, conversion efficiency, per-context yields |
| deduplicated BAM | `results/4.dedup/{sample}.bam.deduplicated.bam` | PCR duplicates removed (input of every downstream step) |
| dedup report | `results/4.dedup/{sample}.dedup_report.txt` | duplication rate |
| cytosine report | `results/5.methylation/{sample}/{sample}.bam.deduplicated.bismark.cov.gz` | per-cytosine methylation (6 columns: chromosome, start, end, methylation %, methylated, unmethylated counts) |
| methylation bedGraph | `results/5.methylation/{sample}/{sample}.bam.deduplicated.bedGraph.gz` | context-level methylation percentage per region |
| splitting report | `results/5.methylation/{sample}/{sample}.bam.deduplicated.splitting_report.txt` | methylation per context (CpG/CHG/CHH) and strand |
| M-bias | `results/5.methylation/{sample}/{sample}.bam.deduplicated.M-bias.txt` | per-read-position methylation (basis for `--mbias` decisions) |
| nucleotide stats | `results/5.methylation/{sample}/{sample}.nucleotide_stats.txt` | input-DNA nucleotide composition of the sample (bam2nuc) |
| merged CpG table | `results/5.methylation/{sample}/{sample}.CpG_merged.tsv.gz` | coverage2cytosine `--merge_CpG` product (CpG sites merged across strands) |
| per-sample report | `results/5.methylation/{sample}/{sample}_seq_context.html` | alignment + dedup + splitting + M-bias + nucleotide stats in one HTML |
| run-level summary | `results/5.QC/bismark2summary.html` | one-row-per-sample overview of all Bismark numbers |
| QC report | `results/5.QC/multiqc/multiqc_report.html` | trimming + FastQC + Bismark reports in one HTML |
| version record | `results/5.QC/software_versions.yaml` | tool versions actually resolved for the run (incl. Snakemake) |
| per-rule logs | `results/logs/` | one log per rule/sample |

Note the derived-name chain: every methylation output carries the `.bam.deduplicated` infix inherited from the input BAM name (see "Bismark output naming" in the [README](../README.md)).

### 7.3 Reading the numbers

1. **Trimming yield** (`*_trimming_report.txt`, MultiQC): the trimmed library size vs raw. BS-seq libraries are adapter-heavy when fragments are shorter than the read length, so a substantial adapter-contaminated fraction is normal — what matters is how many pairs survive `--length` afterwards.
2. **Alignment efficiency** (`{sample}_report.txt`, bismark2summary): unique alignment rate. On a well-matched reference genome, bisulfite libraries typically align somewhat below ordinary DNA-seq rates (the C→T conversion collapses strand complexity); a rate near zero means wrong genome, wrong `--phred33-quals`/`--phred64-quals` encoding, or a non-converted library.
3. **Duplication rate** (`{sample}.dedup_report.txt`): the fraction of duplicate read pairs removed. Bisulfite libraries duplicate heavily (complexity loss from the conversion); judge it together with coverage.
4. **CpG methylation level** (splitting report, cytosine report): the CpG-context methylation percentage should be in the organism's expected range (rice/human endogenous CG methylation is high, tens of percent); near-zero CpG methylation with high CHH is the signature of chloroplast/organellar contamination or a conversion problem. CHG/CHH levels in rice are distinctly non-zero (plant contexts).
5. **M-bias** (`*.M-bias.txt`): methylation percentage per read position. Strong end-of-read biases are the reason people trim the first/last bases with bismark's `--mbias`; inspect before adding such trims to `bismark.align_extra`.

Everything above lands in `results/5.QC/multiqc/multiqc_report.html` (title "BS-seq QC Summary") and `results/5.QC/bismark2summary.html`; the exact tool versions used are recorded in `results/5.QC/software_versions.yaml`. Differential methylation between conditions is out of scope for v0.1 (see [docs/TODO.md](TODO.md)).

---

## 8. FAQ

**Q1: How do I use my own genome instead of the osa/hsa preset?**
Either keep `species: "osa"` (or `"hsa"`) and set an explicit `genome:` path in the project config (explicit keys win over the preset), or set `species: "none"` and set `genome:` — under `none` nothing is inherited, so the key is mandatory. The genome must be an unmodified (non-bisulfite-converted) FASTA; `bismark_genome_prep` performs the conversion itself and caches it under `results/0.index/bismark_genome/`.

**Q2: How do I skip trimming (legacy raw-read alignment)?**
Set `trim.enabled: false`. No trim job is scheduled and `bismark_align` reads straight from `1.rawdata/` — exactly the legacy behavior. Everything downstream is unchanged.

**Q3: Can I mix paired-end and single-end samples?**
Yes. The layout is detected per sample from the `1.rawdata/` file names (§3.2) and each sample's chain (trim → align → dedup → extract) uses the matching flags. The only naming consequence is the trim stage: PE samples produce `{sample}_1_val_1.fq.gz` + `{sample}_2_val_2.fq.gz`, SE samples `{sample}_trimmed.fq.gz`.

**Q4: The methylation extractor was killed for memory. What do I tune?**
The extractor's `--buffer_size` is derived from the rule's `mem_mb` (§5.2). Either lower `methylation_extractor.buffer_frac` (smaller buffer) in the project config, or override the rule's memory in a project `resources.yaml` (`resources: {methylation_extractor: {mem_mb: 64000}}`) so the derived buffer grows with the request — never hard-code a buffer bigger than the scheduler grants.

**Q5: Why are my methylation files called `{sample}.bam.deduplicated.*`?**
Bismark tools derive output names from the input file name: the extractor consumes the deduplicated BAM and inherits its full basename, so every report carries the `.bam.deduplicated` infix. The full chain is tabulated in the [README](../README.md#bismark-output-naming).

**Q6: How do I resume after an interruption? What about lock errors?**
Snakemake skips completed steps by output timestamps — **rerun the same command to resume** (profiles pin keep-going and rerun-incomplete). If you hit a working-directory lock error (`Directory cannot be locked`, usually after a force-killed job), run `bash run.sh -P . --unlock` and rerun.

**Q7: What does `[config warning] genome is still a placeholder: /path/to/...` mean?**
Parse-time validation warns when `genome` still holds the `/path/to/` placeholder path. Dry-runs work anyway; a real run fails at `bismark_genome_prep`. Fill real server paths in the project config before launching.

**Q8: A job was killed by the cluster for exceeding walltime/memory. What now?**
Three fixes, in order of preference: (1) targeted per-rule override in a project `resources.yaml` (§4.4) — `bismark_align`, `methylation_extractor`, and `bismark_genome_prep` are the usual suspects for large genomes; (2) global `--runtime` / `--memory` overrides (§6.3); (3) `--retries N` so sporadic failures retry automatically. To locate the rule: find the failed rule name in `snakemake.logs.txt`, then read `results/logs/<rule>/<sample>.log` for details.

**Q9: How do I validate a new deployment or a workflow change?**

```bash
make check                          # bash -n syntax checks (no snakemake needed)
make lint                           # static suite (missing optional tools are skipped)
bash tests/run_test.sh              # synthetic-data dry-run regression (needs snakemake + python3/PyYAML)
bash tests/run_test.sh --real-run   # end-to-end run + output assertions (needs the full analysis environment)
```

Test data is generated by `tests/make_testdata.py` with a fixed seed (2 x 20 kb chromosomes, 2 samples of paired-end 100 bp reads simulated post-bisulfite); the dry-run baseline DAG is 20 jobs (CI passes `--reads 2000`). Before changing workflow code, read the documentation-sync checklist in [CONTRIBUTING](../CONTRIBUTING.md).
