# srna-seq workflow user guide

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
| Analysis tools | bowtie + bowtie-build (1.x) / Trim Galore (brings Cutadapt + FastQC) / MultiQC |

Launcher bootstrap order: `run.sh` first uses the `python3` on the system PATH (with PyYAML) to resolve `config/software.yaml`, then injects the main environment/tools into the current process — so even on the conda_prefix reuse route, the login-node PATH must have `python3` (snakemake is supplied by the main environment after resolution).

### 1.2 Three ways to set up the environment

**Path 1: fresh server — create the environment from the all-in-one template**

```bash
mamba env create -f workflow/environment.yaml   # environment name srna-seq
conda activate srna-seq
```

This command is executed explicitly by the user; Snakemake never creates or modifies software environments on its own. Pinned versions in the template: snakemake-minimal 7.32.4, bowtie 1.3.1, trim-galore 0.6.10, fastqc 0.11.9, multiqc 1.21.

**Path 2: reuse an existing conda environment on the server**

Copy `config/software.yaml` into the project directory (or edit the repository default) and set:

```yaml
environment:
  type: conda
  conda_prefix: "/share/conda/envs/srna-seq"   # recommended on HPC
  # conda_name: "srna-seq"                     # or resolve by environment name
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

`strict: true` makes the preflight fail loudly when a required executable is missing instead of continuing silently. Ordinary tools (bowtie / bowtie-build / trim_galore / fastqc / multiqc / python3) resolve from PATH automatically; `software.yaml` needs no `tools:` overrides unless a binary lives outside PATH. Like the other v0.1 workflows there is no R section.

### 1.3 Preflight checks

```bash
# Check that all required executables resolve (bowtie, bowtie-build,
# trim_galore, fastqc, multiqc, python3); then exit
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

Six steps to a minimal project (the legacy-style rice example from [example/](../example/README.md)):

**Step 1: create the working directory**

```bash
mkdir -p ~/work/rice/1.rawdata
```

The working directory (data) and the analysis code directory (this repository) are separate; `run.sh` generates all results inside the working directory.

**Step 2: add the raw data**

Single-end FASTQ files named `{sample}.fastq.gz`, `{sample}.fq.gz`, `{sample}_R1.fastq.gz`, or `{sample}_R1.fq.gz` (exact names, checked when the workflow resolves the sample):

```bash
cp root_rep1.fq.gz root_rep2.fq.gz leaf_rep1.fq.gz leaf_rep2.fq.gz ~/work/rice/1.rawdata/
```

**Step 3: write the sample table**

```bash
cp example/samples.csv ~/work/rice/samples.csv
```

Edit the file with your own samples (a single `sample_id` column; validation rules in §3).

**Step 4: write the project config**

```bash
cp example/config.yaml ~/work/rice/config.yaml
```

Fill in the cascade class FASTAs and the genome FASTA, or keep `species: "osa"` and let the `config/species.yaml` preset supply the placeholder paths (an empty `fasta` inherits the preset; see §4.2/§4.3). `run.sh` auto-detects `<project dir>/config.yaml` when `-c` is omitted, so placing it in the working directory is enough. Full key reference in §4.

**Step 5: preflight and dry-run**

```bash
cd /path/to/repo/srna-seq
bash run.sh -P ~/work/rice --check-software   # environment preflight
bash run.sh -P ~/work/rice -n                 # dry-run: builds the DAG and prints planned jobs without executing
```

**Step 6: launch**

```bash
bash run.sh -P ~/work/rice -j 10              # local / auto-detected scheduler
# positional form:
bash run.sh ~/work/rice ~/work/rice/config.yaml 10
# cluster example:
bash run.sh -P ~/work/rice --profile pbs --queue workq --memory 16G --runtime 600
```

After a successful launch the main log is `~/work/rice/snakemake.logs.txt` (change with `--log FILE`), and per-rule logs live under `results/logs/`. See §7 for the results layout.

---

## 3. Sample table in detail

The sample table is a **single-column CSV** (template `config/samples.csv`):

```csv
sample_id
root_rep1
root_rep2
leaf_rep1
leaf_rep2
```

Every row is one sample; v0.1 has no other columns (no group/condition/batch design yet — see the [TODO backlog](TODO.md)). Each `sample_id` must match the prefix of a FASTQ file in `1.rawdata/`.

### 3.1 Validation rules (at parse time, errors carry line numbers)

Validation runs centrally at Snakemake parse time (`load_sample_table` in `workflow/rules/common.smk`); any violation aborts the workflow before any job runs:

1. The header must be **exactly one column named `sample_id`** — otherwise: `Sample table <path> must have exactly one column with header 'sample_id' (got [...]); see config/samples.csv`;
2. `sample_id` must not be empty — `Sample table line <n>: sample_id must not be empty`;
3. `sample_id` allows only **alphanumerics plus `. _ -`**, must start alphanumeric, must not start with `-`, and must not contain consecutive underscores `__` (ids become file names and bowtie command-line arguments) — `Sample table line <n>: sample_id='...' contains illegal characters; only alphanumerics and . _ - are allowed (no leading '-', no '__')`;
4. `sample_id` must be unique — `Sample table line <n>: duplicate sample_id '...'`;
5. The table must have at least one data row — `Sample table <path> has no data rows`.

### 3.2 Raw FASTQ resolution

For each declared sample the workflow looks for exactly one file, in this order: `1.rawdata/{sample}_R1.fastq.gz`, `1.rawdata/{sample}_R1.fq.gz`, `1.rawdata/{sample}.fastq.gz`, `1.rawdata/{sample}.fq.gz`. If none exists, the workflow fails with the supported-name list (`sample <s>: no raw FASTQ found under 1.rawdata/ ...`). The data is expected **single-end small-RNA reads** (typically 18-30 nt after adapter trimming).

---

## 4. Configuration in detail

### 4.1 Config chain and layering

The `workflow/Snakefile` assembles the configuration in this order; **later files override earlier keys**:

1. Repository `config/species.yaml` (species preset map; loaded first);
2. Repository `config/resources.yaml` (per-rule scheduler resources; the project copy is injected by run.sh instead when present);
3. Project config: `-c FILE` explicitly; when omitted, `<project dir>/config.yaml` is auto-detected; when neither exists, the repository `config/config.yaml` defaults apply;
4. Overlay config: `-l FILE` explicitly; when omitted, `<project dir>/config.local.yaml` is auto-layered last (an INFO line is printed when found).

Typical usage: the project config holds project-specific values (reference paths, samples, cascade layout); machine-specific differences (e.g. a lower thread cap) go into `config.local.yaml` (kept out of version control). `config/config.template.yaml` is the annotated template — copy it into the working directory and change only what differs (unlisted keys keep the repository defaults).

Scheduler resources are layered separately: see §4.4 for `config/resources.yaml` and the project-level `resources.yaml`.

### 4.2 All config.yaml keys

| Key | Type | Default | Effect |
|---|---|---|---|
| `results_dir` | string | `"results"` | root directory for all derived outputs, relative to the working directory (raw inputs `1.rawdata/` stay at the root) |
| `SampleListFile` | string | `"config/samples.csv"` | sample-table path; resolves absolute > working-directory-relative > repository-relative (§3) |
| `threads` | int >= 1 | `12` | legacy global thread cap: every rule's thread request (§4.4) is capped at this value |
| `species` | string | `"osa"` | selects a reference preset from `config/species.yaml`; supported values are `osa` (rice preset) and `none` (no preset, no fallback — anything else errors at parse time) |
| `genome.fasta` | string | from the species preset | genome FASTA for the genome bowtie index; **empty = fall back to the species preset; still empty = the genome alignment is skipped** (`/path/to/` placeholders only warn) |
| `trim.quality` | int >= 0 | `25` | 3' quality-trimming cutoff passed to Trim Galore (`-q`) |
| `trim.min_len` | int >= 1 | `15` | discard trimmed reads shorter than this (`--length`) |
| `trim.adapter` | string of `ACGTN` | `"AGATCGGAAGAGC"` | 3' adapter sequence (`-a`) |
| `trim.stringency` | int >= 1 | `3` | adapter overlap required to trigger trimming (`--stringency`) |
| `trim.error_rate` | float in (0, 1] | `0.1` | fraction of mismatches allowed in an adapter match (`-e`) |
| `trim.extra` | string | `""` | extra Trim Galore arguments appended verbatim |
| `cascade` | non-empty ordered list of `{name, fasta}` | 7 entries (rRNA, snoRNA, snRNA, tRNA, miRNA, mRNA, rhizo) | the sncRNA filter cascade; semantics below |
| `bowtie.extra` | string | `""` | extra bowtie1 arguments appended to every alignment (cascade stages and genome alike), e.g. `"-v 1 --best --strata"` |
| `resources` | mapping of rule -> {threads, mem_mb, runtime_min} | `{}` | per-rule scheduler-resource overrides, same shape as `config/resources.yaml` (§4.4) |

**Cascade entry semantics** (the most important config in this workflow):

- The list is **ordered** — classes run strictly top to bottom, and the order defines read assignment priority (a read belongs to the first class whose index it maps to).
- Each entry is `{name, fasta}`. `name` follows the sample-name rules (alphanumerics plus `. _ -`, no leading `-`, no `__`, no duplicates) because it becomes a directory name and a bowtie argument.
- An entry with an **empty `fasta` inherits the species preset** path for that class (`species.yaml` → `cascade.<name>`). If the fasta is **still empty after that fallback** (or the preset has no such class), the entry is **skipped entirely**: no bowtie index, no cascade stage, no counts, and no columns in the count/RPM matrices or the cascade summary.
- With `species: "none"` **no preset merge happens**: only entries with explicitly non-empty fastas run. This is how you run a fully custom class set.
- Validation requires **at least one configured cascade class or a configured genome** — a config with neither has nothing to align and is rejected at parse time (`neither any cascade class nor the genome has a configured fasta — nothing to align`). Keep at least one cascade class configured: the per-class quantification stage is defined over the configured classes.

Validation runs at parse time (`validate_config` in `workflow/rules/common.smk`) and aggregates all problems into a single `WorkflowError` report. Reference files that do not exist (or are still `/path/to/` placeholders) only **print a warning and do not abort** — dry-run/lint often run on machines without the reference files; in a real run a missing reference fails the corresponding rule, so confirm each warning before launching. Skipped classes warn too: `[config warning] cascade class 'x' has no fasta and is skipped`.

### 4.3 Species presets

`config/species.yaml` currently defines one preset:

```yaml
osa:
  label: "Rice Oryza sativa (IRGSP-1.0; sncRNA classes from user references)"
  cascade:
    rRNA: "/path/to/reference/osa/rRNA.fa"
    snoRNA: "/path/to/reference/osa/snoRNA.fa"
    snRNA: "/path/to/reference/osa/snRNA.fa"
    tRNA: "/path/to/reference/osa/tRNA.fa"
    miRNA: "/path/to/reference/osa/miRNA.fa"
    mRNA: "/path/to/reference/osa/mRNA.fa"
    rhizo: "/path/to/reference/osa/Rhizo.fa"
  genome:
    fasta: "/path/to/reference/osa/genome.fa"
```

Preset values fill **only the entries the project config leaves empty**: a cascade entry with an empty `fasta` inherits `preset.cascade.<name>`, and an empty `genome.fasta` inherits `preset.genome.fasta`. An explicit **non-empty** value in the project config wins. Note the asymmetry with the sibling workflows: here an empty string means "unset" and triggers the fallback, so the way to *skip* a preset class is to remove its entry from the cascade list (or switch to `species: "none"`), not to blank its fasta. The shipped paths are placeholders — replace them with real cluster paths at deployment time.

### 4.4 Per-rule scheduler resources (`config/resources.yaml`)

Scheduler resources live in their own file — the single place to tune cluster requests. Each rule declares `threads`, `mem_mb`, and `runtime_min`; `runtime_sec` (used by PBS walltime and SGE h_rt) is derived automatically as `runtime_min x 60` and cannot be set independently. The legacy top-level `threads` in config.yaml still acts as a global cap on every rule's threads.

Override without touching the repository: copy the file into your working directory as `resources.yaml` and edit — `run.sh` auto-detects it. Lookup order: the `SRNA_RESOURCES_CONFIG` environment variable > `<project dir>/resources.yaml` > repository `config/resources.yaml`. Rules not listed in an override keep the repository defaults.

Example override (`<project dir>/resources.yaml`):

```yaml
resources:
  genome_align:
    mem_mb: 16000
    runtime_min: 300
  bowtie_index:
    threads: 16
```

Built-in defaults per rule (identical to `config/resources.yaml`):

| Rule | threads | mem_mb | runtime_min |
|---|---:|---:|---:|
| `software_versions` | 1 | 1024 | 10 |
| `trim` | 4 | 4096 | 60 |
| `fastqc` | 2 | 2048 | 30 |
| `multiqc` | 2 | 4096 | 30 |
| `bowtie_index` | 8 | 8192 | 120 |
| `cascade_stage` | 4 | 8192 | 120 |
| `genome_align` | 8 | 8192 | 180 |
| `count_stage` | 1 | 2048 | 30 |
| `merge_counts` | 1 | 4096 | 30 |
| `cascade_summary` | 1 | 2048 | 30 |

(`fastqc` has its own entry for completeness; in v0.1 FastQC runs inside `rule trim` via `trim_galore --fastqc`.)

### 4.5 software.yaml

`config/software.yaml` (or a project copy) configures the runtime, separate from analysis parameters:

```yaml
environment:
  type: system                  # system | conda
  conda_prefix: ""              # preferred for HPC, e.g. /share/conda/envs/srna-seq
  conda_name: ""                # alternative to conda_prefix
  strict: true

tools: {}          # add overrides only for binaries outside the main PATH
paths: {}
databases: {}
```

All analysis tools (bowtie, bowtie-build, trim_galore, fastqc, multiqc, python3) are ordinary PATH executables resolved through `workflow/scripts/runtime_config.py`; `run.sh` exports the resolved values to the workflow. Unlike seclip-seq there is no external-legacy-tool equivalent — nothing here needs a hand-installed absolute path.

---

## 5. Run modes

### 5.1 Which classes run — and which are skipped

For every cascade entry the parse-time decision is:

| `cascade` entry (after the §4.3 preset fallback) | what enters the DAG |
|---|---|
| `fasta` non-empty | `bowtie_index` for the class + one `cascade_stage` per sample + `count_stage` per sample; the class appears in the count/RPM matrices and as `{class}_mapped`/`{class}_unmapped` columns in `cascade_summary.tsv` |
| `fasta` still empty | **skipped entirely** — no index, no stage, no counts, no matrix rows, no summary columns; a `[config warning]` is printed at parse time |

Combinations:

| Cascade classes | `genome.fasta` | DAG |
|---|---|---|
| >=1 configured | configured | full workflow: cascade filter + genome alignment + quantification + QC (the baseline configuration) |
| >=1 configured | empty (after preset fallback) | cascade-only: no `0.index/genome`, no `genome_align`, `cascade_summary.tsv` reports `genome_unmapped` as `NA` |
| all entries skipped | — | rejected at parse time when the genome is unconfigured too (`nothing to align`); keep at least one class so the quantification stage has input |

Skipped classes are reported, never silent: each one prints `[config warning] cascade class '<name>' has no fasta and is skipped` during the dry-run/preflight.

### 5.2 The genome alignment uses the TRIMMED reads — not the cascade remainder

`rule genome_align` aligns `results/2.cleandata/{sample}_trimmed.fq.gz` (the full trimmed read set), **not** the unmapped output of the last cascade stage. This is deliberate and faithful to the legacy script, and the two views answer different questions:

- the **cascade stages** classify reads against curated class references — reads assigned to `miRNA` never reach the later stages, so the later SAMs (and the cascade remainder) are progressively depleted;
- the **genome SAM** is a complete trimmed-read alignment against the reference genome, unaffected by what the cascade assigned — the basis for genome-level analyses (e.g. novel miRNA discovery, a v0.2 item) and for judging how much of the library the cascade classified at all.

Consequently `genome_unmapped` in `cascade_summary.tsv` is *not* the complement of the last class's unmapped count: it counts trimmed reads that did not match the genome.

### 5.3 Checks, dry-run, and direct Snakemake invocation

```bash
bash run.sh -P . -n                    # dry-run: builds the DAG and prints the jobs it would run, without executing
bash run.sh -P . --validate-only       # parses workflow + config, lists rules, then exits
bash run.sh -P . --check-software      # executable preflight
```

The dry-run automatically skips the software preflight (needs only snakemake + python3/PyYAML, not the full analysis environment). Sample-table and config validation happen at parse time, so the dry-run and `--validate-only` catch the same errors listed in §3.1 / §4.2.

Direct Snakemake invocation also works without the launcher:

```bash
SRNA_CONFIG=/path/to/config.yaml snakemake -s workflow/Snakefile -j 10
```

with the env vars the launcher would normally provide (`SRNA_RESOURCES_CONFIG`, `SRNA_EXTRA_CONFIG`).

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
2. **Project override**: a `resources.yaml` copied into the project directory (auto-detected by run.sh; or point `SRNA_RESOURCES_CONFIG` at any file);
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
# PBS: workq queue, all rules forced to 16G/600min, retry failed jobs once
bash run.sh -P . --profile pbs --queue workq --memory 16G --runtime 600 --retries 1

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
| `SRNA_JOBS` | `-j / --jobs` (default 10) |
| `SRNA_LOG` | `--log` (default `snakemake.logs.txt`) |
| `SRNA_QUEUE` | `--queue` |
| `SRNA_PARTITION` | `--partition` |
| `SRNA_MEMORY` | `--memory` |
| `SRNA_RUNTIME_MIN` | `--runtime` |
| `SRNA_SGE_MEMORY_RESOURCE` | `--sge-mem-resource` (default `h_vmem`) |
| `SRNA_SCHEDULER_EXTRA` | `--scheduler-extra` |
| `SRNA_RETRIES` | `--retries` |
| `SRNA_LATENCY_WAIT` | `--latency-wait` |
| `SRNA_MAX_JOBS_PER_SECOND` | `--max-jobs-per-sec` |
| `SRNA_MAX_STATUS_CHECKS_PER_SECOND` | `--max-status-per-sec` |
| `SRNA_SOFTWARE_CONFIG` | `--software` |
| `SRNA_CONFIG` | `-c / --config` |
| `SRNA_EXTRA_CONFIG` | `-l / --extra-config` |
| `SRNA_RESOURCES_CONFIG` | project `resources.yaml` lookup |

One further variable is internal but also read by direct Snakemake invocation: `SRNA_PYTHON` (interpreter used by the `software_versions` rule).

---

## 7. Output interpretation

### 7.1 Working-directory layout

```
workdir/
├── config.local.yaml        # optional: working-directory config overlay (auto-layered, not committed)
├── samples.csv              # sample table (filename free; pointed to by SampleListFile)
├── snakemake.logs.txt       # launcher/Snakemake main log
├── 1.rawdata/               # raw inputs (the only data directory outside results/):
│                            #   {sample}.fastq.gz / {sample}_R1.fq.gz (SE small-RNA reads)
└── results/                 # ALL derived artifacts (rename via results_dir in config)
    ├── 0.index/
    │   ├── {class}/            # per-class bowtie index ({class}.1.ebwt + siblings)
    │   └── genome/             # genome bowtie index (genome.1.ebwt + siblings)
    ├── 2.cleandata/
    │   ├── {sample}_trimmed.fq.gz             # Trim Galore output (alignment input)
    │   ├── {sample}_trimming_report.txt       # trimming statistics
    │   └── fastqc/             # {sample}_trimmed_fastqc.html / .zip
    ├── 3.align/
    │   ├── filter/{class}/     # {sample}.sam + {sample}_unmapped.fq per cascade stage
    │   └── genome/             # {sample}.sam + {sample}_unmapped.fq (trimmed reads)
    ├── 4.expression/
    │   ├── {class}/            # {sample}_counts.txt, {class}_counts.tsv, {class}_RPM.tsv
    │   └── all_classes_counts.tsv
    ├── 5.QC/
    │   ├── raw_counts/         # {sample}.txt (raw read counts)
    │   ├── cascade_summary.tsv / cascade_summary_mqc.tsv
    │   ├── multiqc/multiqc_report.html
    │   ├── software_versions.yaml
    │   └── logs/               # QC rule logs (e.g. software_versions.log.txt)
    └── logs/                   # per-rule logs (trim/, cascade/{class}/, genome_align/, count/{class}/, ...)
```

### 7.2 Results quick reference

| Result | Path | Meaning |
|---|---|---|
| class bowtie index | `results/0.index/{class}/{class}.1.ebwt` | auto-built from the class FASTA; `.2.ebwt`/`.rev.*` siblings are part of the same index |
| genome bowtie index | `results/0.index/genome/genome.1.ebwt` | auto-built from `genome.fasta` |
| trimmed reads | `results/2.cleandata/{sample}_trimmed.fq.gz` | quality- + adapter-trimmed, `min_len`-filtered; input to the first cascade stage and the genome alignment |
| trimming report | `results/2.cleandata/{sample}_trimming_report.txt` | Trim Galore statistics (aggregated into MultiQC) |
| FastQC | `results/2.cleandata/fastqc/{sample}_trimmed_fastqc.html` | per-base quality of the trimmed reads |
| cascade stage SAM | `results/3.align/filter/{class}/{sample}.sam` | reads assigned to this class at this stage |
| cascade unmapped reads | `results/3.align/filter/{class}/{sample}_unmapped.fq` | reads passed on to the next stage (or never assigned, for the last class) |
| genome SAM | `results/3.align/genome/{sample}.sam` | full **trimmed-read** alignment against the genome (§5.2) |
| genome unmapped reads | `results/3.align/genome/{sample}_unmapped.fq` | trimmed reads with no genome hit |
| per-sample class counts | `results/4.expression/{class}/{sample}_counts.txt` | reads per reference feature; last line `__mapped_total` (§7.3) |
| class count matrix | `results/4.expression/{class}/{class}_counts.tsv` | one row per feature, one column per sample |
| class RPM matrix | `results/4.expression/{class}/{class}_RPM.tsv` | same shape, RPM-normalized (§7.3) |
| combined matrix | `results/4.expression/all_classes_counts.tsv` | long format `class` / `feature` / one column per sample |
| raw read counts | `results/5.QC/raw_counts/{sample}.txt` | raw reads per sample (cascade summary input) |
| cascade summary | `results/5.QC/cascade_summary.tsv` | per-sample read fate across all stages (§7.3) |
| QC report | `results/5.QC/multiqc/multiqc_report.html` | trimming + FastQC + cascade bargraph in one HTML |
| version record | `results/5.QC/software_versions.yaml` | tool versions actually resolved for the run (incl. Snakemake) |
| per-rule logs | `results/logs/` | one log per rule/sample |

### 7.3 Reading the numbers

**Per-sample class counts** (`{sample}_counts.txt`): one line per reference feature of the class FASTA (`feature<TAB>reads`), counting **primary alignments only** (bowtie secondary/hits after the first are flagged `0x100` and excluded), closed by the sentinel line `__mapped_total<TAB><n>` — the total reads assigned to this class in this stage. `__mapped_total` is the denominator of the RPM normalization and the `{class}_mapped` value in the cascade summary.

**RPM definition** (`{class}_RPM.tsv`): RPM is **per class and per sample** — `count / __mapped_total * 1e6` where `__mapped_total` is that same class and sample's assigned-read total (not all reads of the library). A class that mapped no reads gets `0.000`. Cross-class comparison of RPM values therefore answers "how abundant is this feature within its class", not across classes; use the raw counts in `{class}_counts.tsv` / `all_classes_counts.tsv` for abundance comparisons.

**`cascade_summary.tsv`** — one row per sample, one column pair per configured class:

| column group | meaning |
|---|---|
| `raw` | raw reads in `1.rawdata/` (lines / 4) |
| `trimmed` | reads written passing filters, from the Trim Galore report |
| `{class}_mapped` | reads assigned to the class (`__mapped_total` of the stage counts) |
| `{class}_unmapped` | reads in the stage's unmapped FASTQ (lines / 4) — the input of the next stage |
| `genome_unmapped` | reads in `3.align/genome/{sample}_unmapped.fq` |

Special values: `-1` means the corresponding unmapped FASTQ file was **missing** (an upstream stage failed — investigate the logs, the count could not be taken); `NA` in `genome_unmapped` means the **genome is not configured** (§5.1) so the column is structurally empty, not failed.

Sanity checks on a healthy run:

1. `raw` matches your FASTQ read counts; `trimmed` is the usable library (for small RNA a large trim loss usually means wrong `trim.adapter` or heavy degradation);
2. the `{class}_unmapped` chain decreases monotonically stage by stage (each stage removes its class) — the same number appears as the next stage's stage input;
3. `{class}_mapped` sums to the fraction of the library the cascade could classify; reads in neither class nor genome stay in the last `{class}_unmapped`;
4. `genome_unmapped` + genome-mapped reads (from the genome SAM) reconstructs the `trimmed` total (§5.2 explains why the genome step starts from `trimmed`).

Everything above lands in `results/5.QC/multiqc/multiqc_report.html` (title "sRNA-seq QC Summary"), including the "sRNA cascade read fate" bargraph built from `cascade_summary_mqc.tsv`; the exact tool versions used are recorded in `results/5.QC/software_versions.yaml`.

---

## 8. FAQ

**Q1: How do I add a new sncRNA class (e.g. piRNA, or another exogenous index)?**
Append an entry to the `cascade:` list at the position where it should run — order matters, because it defines assignment priority:

```yaml
  - name: piRNA
    fasta: "/path/to/reference/osa/piRNA.fa"
```

The name follows the §4.2 naming rules. Rerun the same command: the index is built automatically and only the new class's jobs (plus the merge/summary/MultiQC refresh) execute. No other config change is needed.

**Q2: How do I drop a class (e.g. we have no rhizo data)?**
Remove the entry from the `cascade:` list — do **not** just blank its `fasta`: with a real species preset an empty fasta *inherits* the preset path (§4.3), so the class would still run. Removing the entry (or switching to `species: "none"` and listing exactly the classes you want) drops it from the DAG; the old `0.index/<class>/` and `3.align/filter/<class>/` outputs are simply no longer requested and can be deleted by hand.

**Q3: Why does the genome alignment use the trimmed reads instead of the cascade remainder?**
It reproduces the legacy script's behaviour on purpose. The cascade and the genome step are two independent views of the same trimmed library: the cascade classifies reads against curated class references (each stage depletes the remainder), while the genome SAM is a complete trimmed-read alignment usable for genome-level analyses regardless of what the cascade assigned. Chaining the genome step after the cascade would make `genome_unmapped` meaningless as a library-level number (§5.2).

**Q4: Can I run paired-end (PE) data?**
No. v0.1 is single-end only: the sample-table parser accepts exactly one `sample_id` column and the FASTQ resolver looks only for `{sample}[_R1]{.fastq,.fq}.gz` files. PE support would need a new sample-table design and is not scheduled for v0.1.

**Q5: How do I resume after an interruption? What about lock errors?**
Snakemake skips completed steps by output timestamps — **rerun the same command to resume** (profiles pin keep-going and rerun-incomplete). If you hit a working-directory lock error (`Directory cannot be locked`, usually after a force-killed job), run `bash run.sh -P . --unlock` and rerun.

**Q6: A class is missing from my results / the cascade summary has fewer columns than expected. Why?**
The class was skipped at parse time: its `fasta` was still empty after the species-preset fallback. Every skipped class prints `[config warning] cascade class '<name>' has no fasta and is skipped` during the run. Check the config's `cascade:` list and the `species.yaml` preset (§4.3); note the summary columns are generated from the *configured* classes only.

**Q7: What does `-1` (or `NA`) in `cascade_summary.tsv` mean?**
`-1` = the corresponding unmapped FASTQ file was missing when the summary was built — an upstream stage failed, so look for the failing rule in `snakemake.logs.txt` and its log under `results/logs/` before trusting any numbers. `NA` in `genome_unmapped` = the genome is unconfigured (structural, not a failure) (§7.3).

**Q8: A job was killed by the cluster for exceeding walltime/memory. What now?**
Three fixes, in order of preference: (1) targeted per-rule override in a project `resources.yaml` (§4.4) — `genome_align` and `bowtie_index` are the usual suspects for large genomes; (2) global `--runtime` / `--memory` overrides (§6.3); (3) `--retries N` so sporadic failures retry automatically. To locate the rule: find the failed rule name in `snakemake.logs.txt`, then read `results/logs/<rule>/<sample>.log` for details.

**Q9: What does `[config warning] ... is still a placeholder: /path/to/...` mean?**
Parse-time validation warns when a configured class fasta or `genome.fasta` still holds the `/path/to/` placeholder path. Dry-runs work anyway; a real run fails at the `bowtie_index` rule. Fill real server paths in the project config before launching.

**Q10: How do I validate a new deployment or a workflow change?**

```bash
make check                          # bash -n syntax checks (no snakemake needed)
make lint                           # static suite (missing optional tools are skipped)
bash tests/run_test.sh              # synthetic-data dry-run regression (needs snakemake + python3/PyYAML)
bash tests/run_test.sh --real-run   # end-to-end run + output assertions (needs the full analysis environment)
```

Test data is generated by `tests/make_testdata.py` with a fixed seed (2 x 20 kb chromosomes, rRNA/tRNA/miRNA references, 2 samples); the dry-run baseline DAG is 27 jobs (CI passes `--reads 2000`). Before changing workflow code, read the documentation-sync checklist in [CONTRIBUTING](../CONTRIBUTING.md).
