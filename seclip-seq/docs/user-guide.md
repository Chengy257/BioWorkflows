# seclip-seq workflow user guide

> Updated: 2026-09-08 (v0.2.0)
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
| Analysis tools | STAR / samtools / bedtools / umi-tools / cutadapt / seqkit / bgzip / FastQC / MultiQC / PureCLIP; plus CLIPper when peak clusters are wanted (external install, see §4.5). bedtools (2.31.0) is only needed when the optional `reproducible_peaks` / `annotate_peaks` stages are enabled (§5.5) |

Launcher bootstrap order: `run.sh` first uses the `python3` on the system PATH (with PyYAML) to resolve `config/software.yaml`, then injects the main environment/tools into the current process — so even on the conda_prefix reuse route, the login-node PATH must have `python3` (snakemake is supplied by the main environment after resolution).

### 1.2 Three ways to set up the environment

**Path 1: fresh server — create the environment from the all-in-one template**

```bash
mamba env create -f workflow/environment.yaml   # environment name seclip-seq
conda activate seclip-seq
```

This command is executed explicitly by the user; Snakemake never creates or modifies software environments on its own. Pinned versions in the template: snakemake-minimal 7.32.4, star 2.7.10b, samtools 1.17, bedtools 2.31.0, cutadapt 4.6, umi-tools 1.1.5 (pip; python 3.9 — see the notes in `workflow/environment.yaml`), seqkit 2.13.0, fastqc 0.11.9, multiqc 1.21, pureclip 1.3.1.

**Path 2: reuse an existing conda environment on the server**

Copy `config/software.yaml` into the project directory (or edit the repository default) and set:

```yaml
environment:
  type: conda
  conda_prefix: "/share/conda/envs/seclip-seq"   # recommended on HPC
  # conda_name: "seclip-seq"                     # or resolve by environment name
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

`strict: true` makes the preflight fail loudly when a required executable is missing instead of continuing silently. Unlike the chip workflow there is no R section — seclip-seq has no R dependency.

### 1.3 Preflight checks

```bash
# Check that all required executables resolve (STAR, umi-tools, cutadapt,
# seqkit, bgzip, samtools, fastqc, multiqc, pureclip) and report the
# configured CLIPper; then exit
bash run.sh -P /path/to/workdir --check-software

# Parse the workflow + all configs (sample-table/config validation runs at parse time),
# list the rules, then exit without starting jobs
bash run.sh -P /path/to/workdir -c /path/to/workdir/config.yaml --validate-only
```

Notes:

- `--check-software` requires `-P` to name the project directory (created automatically when missing);
- before a real run, `run.sh` executes a full preflight **automatically**; `--skip-software-check` skips it, and `-n` (dry-run) skips it automatically;
- CLIPper absence is a **warning, not an error**: the preflight reports `[runtime] [WARN] CLIPper not configured (software.yaml paths.clipper); CLIPper peak calling will be skipped` and the workflow runs PureCLIP only (§5.2);
- `software.yaml` lookup order: `--software FILE` > `<project dir>/software.yaml` > repository `config/software.yaml`.

### 1.4 Snakemake version matrix

| snakemake version | support | notes |
|---|---|---|
| **7.32.4** | reference version | pinned in `workflow/environment.yaml`; cluster profiles (pbs/sge/slurm) target the 7.x classic `--cluster` interface |
| **8.x** | cluster semantics unverified | parsing/dry-run work; 8.x replaced the classic `--cluster` interface with the executor plugin system — test on a small scale before real cluster runs |
| <7 or >=9 | unverified | test before use |

---

## 2. Quick start

Six steps to a minimal project (the legacy-style human FBL example from [example/](../example/README.md)):

**Step 1: create the working directory**

```bash
mkdir -p ~/work/fbl/1.rawdata
```

The working directory (data) and the analysis code directory (this repository) are separate; `run.sh` generates all results inside the working directory.

**Step 2: add the raw data**

Single-end FASTQ files named `{sample}.fastq.gz`, `{sample}.fq.gz`, `{sample}_R1.fastq.gz`, or `{sample}_R1.fq.gz` (exact names, checked at parse time):

```bash
cp FC_rep1.fq.gz FC_rep2.fq.gz ~/work/fbl/1.rawdata/
```

**Step 3: write the sample table**

```bash
cp example/samples.csv ~/work/fbl/samples.csv
```

Edit the file with your own samples (a single `sample_id` column is enough; the optional `condition`/`role` grouping columns are described in §3; validation rules in §3.1).

**Step 4: write the project config**

```bash
cp example/config.yaml ~/work/fbl/config.yaml
```

Fill in the three reference paths (`genome`, `gtf`, `repeats_fa`) or keep `species: "hsa"` and let the `config/species.yaml` preset supply them (explicit keys in your config win over the preset; see §4.2/§4.3). `run.sh` auto-detects `<project dir>/config.yaml` when `-c` is omitted, so placing it in the working directory is enough. Full key reference in §4.

**Step 5: preflight and dry-run**

```bash
cd /path/to/repo/seclip-seq
bash run.sh -P ~/work/fbl --check-software   # environment preflight
bash run.sh -P ~/work/fbl -n                 # dry-run: builds the DAG and prints planned jobs without executing
```

**Step 6: launch**

```bash
bash run.sh -P ~/work/fbl -j 10              # local / auto-detected scheduler
# positional form:
bash run.sh ~/work/fbl ~/work/fbl/config.yaml 10
# cluster example:
bash run.sh -P ~/work/fbl --profile pbs --queue workq --memory 16G --runtime 600
```

After a successful launch the main log is `~/work/fbl/snakemake.logs.txt` (change with `--log FILE`), and per-rule logs live under `results/logs/`. See §7 for the results layout.

---

## 3. Sample table in detail

The sample table is a CSV with one of two accepted headers (template `config/samples.csv`):

The classic **single-column** form (all you need for per-sample peak calling):

```csv
sample_id
FC_rep1
FC_rep2
```

The **three-column** form (optional since v0.2; drives the `reproducible_peaks` consensus stage, §5.5):

```csv
sample_id,condition,role
FC_rep1,FBL,ip
FC_rep2,FBL,ip
```

`condition` groups replicates that measure the same target (the value becomes a file name, so it obeys the same naming rules as `sample_id`); `role` marks the channel and accepts exactly `ip` or `input` — input controls pair with the ip replicates of their condition, and they join the peak-calling target set only when `reproducible_peaks.input_control: true` (§5.5; otherwise only the ip samples are peak-called once the consensus stage is enabled). Each `sample_id` must match the prefix of a FASTQ file in `1.rawdata/`. The grouping columns only take effect when `reproducible_peaks.enabled: true` (§5.5); a single-column table stays fully valid and the two v0.2 stages error out with a clear message if enabled without them.

### 3.1 Validation rules (at parse time, errors carry line numbers)

Validation runs centrally at Snakemake parse time (`load_sample_table` / `_read_sample_table` in `workflow/rules/common.smk`); any violation aborts the workflow before any job runs:

1. The header must be **exactly `sample_id`** or **exactly `sample_id,condition,role`** — otherwise: `Sample table <path> header must be exactly ['sample_id'] or exactly ['sample_id', 'condition', 'role'] (got [...]); see config/samples.csv` (a partially extended header, extra columns, or reordered columns are all rejected);
2. `sample_id` must not be empty — `Sample table line <n>: sample_id must not be empty`;
3. `sample_id` allows only **alphanumerics plus `. _ -`**, must start alphanumeric, must not start with `-`, and must not contain consecutive underscores `__` (ids become file names and STAR command-line arguments) — `Sample table line <n>: sample_id='...' contains illegal characters; only alphanumerics and . _ - are allowed (no leading '-', no '__')`;
4. `sample_id` must be unique — `Sample table line <n>: duplicate sample_id '...'`;
5. With the grouping columns present: `condition` must be non-empty and follow the same character rules as `sample_id` (`condition='...' contains illegal characters ...`), and `role` must be exactly `ip` or `input` (`role='...' must be 'ip' or 'input'`);
6. The table must have at least one data row — `Sample table <path> has no data rows`.

### 3.2 Raw FASTQ resolution

For each declared sample the workflow looks for exactly one file, in this order: `1.rawdata/{sample}_R1.fastq.gz`, `1.rawdata/{sample}_R1.fq.gz`, `1.rawdata/{sample}.fastq.gz`, `1.rawdata/{sample}.fq.gz`. If none exists, parsing fails with the supported-name list (`sample <s>: no raw FASTQ found under 1.rawdata/ ...`). The data is expected **single-end with a 10-base random UMI at the start of each read** (pattern configurable via `umi.pattern`).

---

## 4. Configuration in detail

### 4.1 Config chain and layering

The `workflow/Snakefile` assembles the configuration in this order; **later files override earlier keys**:

1. Repository `config/species.yaml` (species preset map; loaded first);
2. Repository `config/resources.yaml` (per-rule scheduler resources; the project copy is injected by run.sh instead when present);
3. Project config: `-c FILE` explicitly; when omitted, `<project dir>/config.yaml` is auto-detected; when neither exists, the repository `config/config.yaml` defaults apply;
4. Overlay config: `-l FILE` explicitly; when omitted, `<project dir>/config.local.yaml` is auto-layered last (an INFO line is printed when found).

Typical usage: the project config holds project-specific values (reference paths, samples, switches); machine-specific differences (e.g. a lower thread cap) go into `config.local.yaml` (kept out of version control). `config/config.template.yaml` is the annotated template — copy it into the working directory and change only what differs (unlisted keys keep the repository defaults).

Scheduler resources are layered separately: see §4.4 for `config/resources.yaml` and the project-level `resources.yaml`.

### 4.2 All config.yaml keys

| Key | Type | Default | Effect |
|---|---|---|---|
| `results_dir` | string | `"results"` | root directory for all derived outputs, relative to the working directory (raw inputs `1.rawdata/` stay at the root) |
| `SampleListFile` | string | `"config/samples.csv"` | sample-table path; resolves absolute > working-directory-relative > repository-relative (§3) |
| `threads` | int >= 1 | `12` | legacy global thread cap: every rule's thread request (§4.4) is capped at this value |
| `species` | string | `"hsa"` | selects a reference preset from `config/species.yaml`; currently only `hsa` is defined (anything else errors at parse time) |
| `genome` | string | from the species preset | genome FASTA for the STAR genome index and PureCLIP (`-g`); `/path/to/` placeholders only warn |
| `gtf` | string | from the species preset | gene annotation for the STAR genome index junction database |
| `repeats_fa` | string | from the species preset | sncRNA/repeats FASTA for the pre-filter index; **required (non-empty) when `filter_repeats: true`** |
| `filter_repeats` | bool | `true` | align against the sncRNA/repeats index first and keep only unmapped reads for genome alignment (see §5.1) |
| `umi.pattern` | string of `ACGTN` | `"NNNNNNNNNN"` | UMI pattern for `umi_tools extract` (10 random bases at the read start in the legacy library) |
| `cutadapt.min_len` | int >= 1 | `18` | discard trimmed reads shorter than this |
| `cutadapt.quality_cutoff` | int >= 0 | `6` | 3' quality-trimming cutoff (`-q`) |
| `cutadapt.error_rate` | float in (0, 1] | `0.1` | fraction of mismatches allowed in an adapter match (`-e`) |
| `cutadapt.adapters` | non-empty list of `ACGT` strings | 20 shifted 3' adapter variants (legacy set) | adapter sequences passed as `-a` to both cutadapt passes |
| `star.sjdb_overhang` | int >= 1 | `139` | junction-database overhang for the genome index (`ReadLength - 1`) |
| `star.genome_sa_index_nbases` | int >= 1 | `14` | `--genomeSAindexNbases` for the genome index (small genomes need a smaller value, see FAQ Q3) |
| `star.repeats_sa_index_nbases` | int >= 1 | `12` | `--genomeSAindexNbases` for the repeats index (small reference) |
| `star.repeats_limit_ram` | int > 0 | `500000000000` | `--limitGenomeGenerateRAM` for the repeats index build |
| `star.filter_multimap_nmax` | int >= 1 | `30` | `--outFilterMultimapNmax` when filtering against repeats |
| `star.align_multimap_nmax` | int >= 1 | `1` | `--outFilterMultimapNmax` for the genome alignment (1 = unique reads only) |
| `callpeak.pureclip` | bool | `true` | run PureCLIP on the deduplicated BAM; `false` removes it from the target list |
| `callpeak.clipper` | bool | `true` | request CLIPper peak clusters; effective only when a CLIPper executable is configured (§5.2) |
| `callpeak.clipper_species` | string | `"GRCh38_v40"` | CLIPper `--species` value; must be non-empty when `callpeak.clipper: true` |
| `reproducible_peaks.enabled` | bool | `false` | optional v0.2 stage (§5.5): build a per-condition consensus of the ip-role PureCLIP beds (`bedtools multiinter`); requires the `condition`/`role` sample-table columns and `callpeak.pureclip: true` |
| `reproducible_peaks.min_replicates` | int >= 1 | `2` | a consensus site must be present in at least this many ip samples of the condition; every condition must actually have that many ip samples (parse-time error otherwise) |
| `reproducible_peaks.input_control` | bool | `false` | with the consensus stage enabled: also peak-call the `role: input` samples and build a per-condition background union from their PureCLIP beds; the consensus of conditions with inputs gains the binary `in_input_background` flag column (§5.5) |
| `reproducible_peaks.filter_by_input` | bool | `false` | requires `input_control: true`: additionally write `{condition}.consensus.filtered.bed` with the input-background-flagged sites removed and feed the filtered BED to the consensus annotation (§5.5) |
| `annotate_peaks.enabled` | bool | `false` | optional v0.2 stage (§5.5): annotate every peak set (per-sample PureCLIP + consensus) with nearest gene, distance, and biotype from the configured `gtf`; requires `callpeak.pureclip: true` |
| `resources` | mapping of rule -> {threads, mem_mb, runtime_min} | `{}` | per-rule scheduler-resource overrides, same shape as `config/resources.yaml` (§4.4) |

Validation runs at parse time (`validate_config` in `workflow/rules/common.smk`) and aggregates all problems into a single `WorkflowError` report. Reference files that do not exist (or are still `/path/to/` placeholders) only **print a warning and do not abort** — dry-run/lint often run on machines without the reference files; in a real run a missing reference fails the corresponding rule, so confirm each warning before launching.

### 4.3 Species presets

`config/species.yaml` currently defines one preset:

```yaml
hsa:
  label: "Human Homo sapiens (GRCh38.p14 / GENCODE v44)"
  genome: "/path/to/reference/hsa/GRCh38.p14.genome.onlychr.fa"
  gtf: "/path/to/reference/hsa/gencode.v44.annotation.gtf"
  repeats_fa: "/path/to/reference/hsa/sncRNA_gencode_Rfam.fa"
```

Preset values fill only the keys the project config does **not** set; an explicit `genome`/`gtf`/`repeats_fa` in the project config wins. The shipped paths are placeholders — replace them with real cluster paths at deployment time.

### 4.4 Per-rule scheduler resources (`config/resources.yaml`)

Scheduler resources live in their own file — the single place to tune cluster requests. Each rule declares `threads`, `mem_mb`, and `runtime_min`; `runtime_sec` (used by PBS walltime and SGE h_rt) is derived automatically as `runtime_min x 60` and cannot be set independently. The legacy top-level `threads` in config.yaml still acts as a global cap on every rule's threads.

Override without touching the repository: copy the file into your working directory as `resources.yaml` and edit — `run.sh` auto-detects it. Lookup order: the `SECLIP_RESOURCES_CONFIG` environment variable > `<project dir>/resources.yaml` > repository `config/resources.yaml`. Rules not listed in an override keep the repository defaults.

Example override (`<project dir>/resources.yaml`):

```yaml
resources:
  star_align:
    mem_mb: 48000
    runtime_min: 600
  callpeak_pureclip:
    threads: 16
```

Built-in defaults per rule (identical to `config/resources.yaml` for the always-on rules; the three optional-stage rules at the bottom ship as workflow defaults only — a project `resources.yaml` can override them by name like any other rule):

| Rule | threads | mem_mb | runtime_min |
|---|---:|---:|---:|
| `software_versions` | 1 | 1024 | 10 |
| `fastqc` | 2 | 2048 | 30 |
| `multiqc` | 2 | 4096 | 60 |
| `star_index_genome` | 12 | 48000 | 480 |
| `star_index_repeats` | 12 | 32000 | 240 |
| `umi_extract` | 4 | 8000 | 180 |
| `cutadapt_trim` | 8 | 8000 | 180 |
| `fastq_sort` | 4 | 8000 | 120 |
| `star_filter_repeats` | 12 | 32000 | 360 |
| `star_align` | 12 | 32000 | 360 |
| `umi_dedup` | 4 | 16000 | 360 |
| `read_count` | 1 | 2000 | 10 |
| `callpeak_clipper` | 4 | 16000 | 240 |
| `callpeak_pureclip` | 8 | 16000 | 360 |
| `consensus_peaks` (optional stage) | 1 | 2048 | 30 |
| `input_background` (optional stage) | 1 | 2048 | 30 |
| `flag_input_background` (optional stage) | 1 | 2048 | 30 |
| `filter_input_background` (optional stage) | 1 | 2048 | 30 |
| `gtf_gene_regions` (optional stage) | 1 | 4096 | 60 |
| `annotate_peaks` (optional stage, both annotate rules) | 1 | 4096 | 60 |

### 4.5 software.yaml

`config/software.yaml` (or a project copy) configures the runtime, separate from analysis parameters:

```yaml
environment:
  type: system                  # system | conda
  conda_prefix: ""              # preferred for HPC, e.g. /share/conda/envs/seclip-seq
  conda_name: ""                # alternative to conda_prefix
  strict: true

tools:
  clipper: ""        # absolute path to an existing CLIPper install, e.g.
                     # /home/chengyu/soft/miniconda3/envs/clipper/bin/clipper
                     # empty = CLIPper peak calling is skipped (PureCLIP only)

paths: {}
databases: {}
```

`paths.clipper` semantics: CLIPper is an external legacy install (Python2-era upstream, deliberately not packaged). `run.sh` resolves it through `workflow/scripts/runtime_config.py` and exports the absolute path as `SECLIP_TOOL_CLIPPER`; the workflow defines the `callpeak_clipper` rule only when the path is non-empty. Ordinary tools (STAR, umi-tools, cutadapt, seqkit, FastQC, MultiQC, PureCLIP, samtools, bgzip) are resolved from the main environment/PATH automatically. When invoking Snakemake directly without run.sh (§5.3), you can export `SECLIP_TOOL_CLIPPER` yourself.

---

## 5. Run modes

### 5.1 `filter_repeats: true / false` — what changes in the DAG

| | `filter_repeats: true` (default) | `filter_repeats: false` |
|---|---|---|
| `0.index` | genome index **and** repeats index (`star_index_repeats`) | genome index only |
| `3.align` | `star_filter_repeats` aligns the sorted reads to the repeats index; only its unmapped reads (`3.align/repeats/{sample}_Unmapped.out.mate1`) enter `star_align` | `star_align` consumes the sorted trimmed reads directly; the repeats filter rule does not exist |
| Effect | sncRNA/repeats-derived reads never reach the genome alignment or the peaks | all trimmed reads map against the genome; sncRNA-derived reads appear in the BAM |

This switch replaces the legacy repo's second Snakefile variant (`seCLIP_nofilter.smk`). Validation: `filter_repeats: true` requires a non-empty `repeats_fa` (species preset or project config); with `false`, no repeats reference is needed at all.

### 5.2 CLIPper present / absent

`callpeak.clipper: true` (default) only *requests* CLIPper peaks. Whether they are produced depends on `paths.clipper` in software.yaml:

- **Configured**: `callpeak_clipper` runs per sample and `results/5.callpeak/{sample}.clipper.peakClusters.bed` joins the targets;
- **Empty (default)**: the rule is excluded from the DAG at parse time, a `[config warning] callpeak.clipper is enabled but no CLIPper executable is configured ...` line is printed, and PureCLIP runs alone;
- `callpeak.pureclip: false` removes PureCLIP targets the same way; with both callers off the workflow stops after `4.rmdup` (plus QC).

### 5.3 Checks, dry-run, and direct Snakemake invocation

```bash
bash run.sh -P . -n                    # dry-run: builds the DAG and prints the jobs it would run, without executing
bash run.sh -P . --validate-only       # parses workflow + config, lists rules, then exits
bash run.sh -P . --check-software      # executable + CLIPper + database preflight
```

The dry-run automatically skips the software preflight (needs only snakemake + python3/PyYAML, not the full analysis environment). Sample-table and config validation happen at parse time, so the dry-run and `--validate-only` catch the same errors listed in §3.1 / §4.2.

Direct Snakemake invocation also works without the launcher:

```bash
SECLIP_CONFIG=/path/to/config.yaml snakemake -s workflow/Snakefile -j 10
```

with the env vars the launcher would normally provide (`SECLIP_RESOURCES_CONFIG`, `SECLIP_EXTRA_CONFIG`, `SECLIP_TOOL_CLIPPER`).

### 5.4 Resuming and Snakemake passthrough

- Resuming: Snakemake skips completed steps based on output timestamps; after an interruption, **simply rerun the same command**;
- profiles pin `keep-going: true` and `rerun-incomplete: true`; `latency-wait` defaults to 90 (pbs) or 60 (default/sge/slurm), overridable via `--latency-wait SEC`;
- a stale working-directory lock (e.g. a job was force-killed) is cleared with `bash run.sh -P . --unlock`;
- arguments after `--` pass through to Snakemake untouched, e.g. `bash run.sh -P . -- --rerun-triggers mtime`;
- the launcher log defaults to `snakemake.logs.txt` (change with `--log FILE`); `-q` reduces the launcher's own output.

### 5.5 Optional v0.2 stages: `reproducible_peaks` and `annotate_peaks`

Both stages are **default off** and drop out of the DAG entirely while disabled — a v0.1-style config (single-column sample table, no new sections) produces exactly the v0.1 pipeline.

**Cross-sample reproducible peaks** (`reproducible_peaks.enabled: true`): v0.1/v0.2 call peaks per sample; this stage adds the cross-replicate view. Requires the `condition`/`role` sample-table columns (§3) and `callpeak.pureclip: true`.

```yaml
reproducible_peaks:
  enabled: true
  min_replicates: 2
```

One `consensus_peaks` job runs per condition: the PureCLIP beds of that condition's ip-role samples are coordinate-sorted and merged with `bedtools multiinter`; sites present in >= `min_replicates` input beds are kept and the support count (how many replicates carry the site) is written to column 4 of `results/6.reproducible_peaks/{condition}.consensus.bed`. Parse-time validation errors out when the table has no `condition`/`role` columns, when `callpeak.pureclip` is off, or when a condition has fewer ip samples than `min_replicates` (the support threshold could never be reached).

**Input controls** (`reproducible_peaks.input_control: true`, needs `enabled: true` and the `condition`/`role` columns; closes backlog item 3):

```yaml
reproducible_peaks:
  enabled: true
  min_replicates: 2
  input_control: true       # peak-call the role=input samples too
  filter_by_input: false    # optionally drop flagged consensus sites
```

The upstream chain never changes — trim/align/dedup always run for every declared sample — but the PureCLIP target set does: without `input_control` only the ip samples of each condition are peak-called, with `input_control` the role=input samples join as well. Their beds feed one `input_background` job per condition with inputs: a `bedtools multiinter` union (every row has support >= 1 by construction, so the whole union is kept; column 4 reports how many input controls cover a feature) written to `results/6.reproducible_peaks/{condition}.input_background.bed`.

For every condition that declares input samples, the final consensus keeps the W7 columns 1-4 (chrom, start, end, support) and appends the binary `in_input_background` flag as column 5 (1 = the site overlaps the condition's input background, 0 = ip-specific). The background is **condition-scoped**: an input control only flags the consensus of the ip replicates declared under the same `condition`. Conditions without input samples keep the plain unflagged BED4 — byte-identical to the no-`input_control` output (parse-time warning when `input_control` is on). With `filter_by_input: true`, an additional `results/6.reproducible_peaks/{condition}.consensus.filtered.bed` drops the flagged sites (BED4 again) and the annotation stage consumes the filtered BED for those conditions; without filtering it annotates the flagged BED5 and drops the flag column per its column contract (`score` stays the consensus support).

Design rationale (recorded when backlog item 3 was closed): the background is condition-scoped because an input control only pairs meaningfully with the ip replicates of the same condition; the union uses the lowest possible threshold (support >= 1 — every interval ever seen in any input counts as background); and flagging is the default over filtering so both the annotated-with-flag view (all reproducible sites, background membership recorded) and the cleaned view (opt-in `filter_by_input`) remain available.

Validation additions (parse-time, aggregated): `input_control` requires the `condition`/`role` columns and `reproducible_peaks.enabled: true` (a warning while the stage is off); `filter_by_input` requires `input_control`; a condition with input samples but no ip samples is a hard error (there is no ip consensus to build or flag); an ip condition without inputs only warns (the background is simply absent).

**GTF-based peak annotation** (`annotate_peaks.enabled: true`): annotates every existing peak set — each per-sample PureCLIP bed plus each consensus bed. Requires `callpeak.pureclip: true` (it annotates PureCLIP output; CLIPper beds are not annotated).

```yaml
annotate_peaks:
  enabled: true
```

Three rule kinds join the DAG: `gtf_gene_regions` parses the configured `gtf` (the same file the STAR index uses) into `results/6.annotation/_ref/genes.bed` + `exons.bed` + a `gene_id/gene_name/gene_biotype` table (stdlib parser, no extra dependency); then one job per peak set runs `bedtools intersect -u` against exons and gene bodies (feature classification: `exon` = overlaps an exon, `gene` = inside a gene body but not an exon, `intergenic` = neither) and `bedtools closest -d -t first` (nearest gene + signed distance), merged by the stdlib script `workflow/scripts/annotate_peaks.py` into `results/6.annotation/{set}.annotation.tsv` with the columns:

```
chrom  start  end  score  nearest_gene  nearest_gene_id  distance  feature_class  gene_biotype
```

`score` is the PureCLIP crosslink-site score for sample sets and the consensus support for consensus sets (`{set} = {condition}.consensus`); `distance` is 0 for overlapping peaks, negative when the nearest gene lies upstream of the peak, and `NA` when the peak's contig carries no gene at all. Peak coordinates stay in the BED system (0-based start) so results round-trip against the peak files.

Typical combined setup (the legacy five-replicate FBL example in `example/`):

```yaml
# samples.csv
sample_id,condition,role
FC_rep1,FBL,ip
FC_rep2,FBL,ip
...
```

```yaml
reproducible_peaks:
  enabled: true
  min_replicates: 3     # consensus sites must replicate in >= 3 of the 5 FBL replicates
annotate_peaks:
  enabled: true
```

Dry-run first (`bash run.sh -P . -n`): the consensus and annotation rules appear in the DAG as `consensus_peaks`, `gtf_gene_regions`, `annotate_sample_peaks`, and `annotate_consensus_peaks`; with `input_control` the raw consensus becomes the `consensus_peaks_raw` intermediate and `input_background`, `flag_input_background`, and (with `filter_by_input`) `filter_input_background` join the DAG. Scheduler resources for the new rules default to 1 thread / 2048 MB / 30 min (`consensus_peaks` and the three input-control rules) and 1 / 4096 / 60 (GTF prep and annotation); override them per rule via a project `resources.yaml` (§4.4). The regression script exercises both scenarios with synthetic data: `bash tests/run_test.sh --consensus` and `bash tests/run_test.sh --input-control`.

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
2. **Project override**: a `resources.yaml` copied into the project directory (auto-detected by run.sh; or point `SECLIP_RESOURCES_CONFIG` at any file);
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
| `SECLIP_JOBS` | `-j / --jobs` (default 10) |
| `SECLIP_LOG` | `--log` (default `snakemake.logs.txt`) |
| `SECLIP_QUEUE` | `--queue` |
| `SECLIP_PARTITION` | `--partition` |
| `SECLIP_MEMORY` | `--memory` |
| `SECLIP_RUNTIME_MIN` | `--runtime` |
| `SECLIP_SGE_MEMORY_RESOURCE` | `--sge-mem-resource` (default `h_vmem`) |
| `SECLIP_SCHEDULER_EXTRA` | `--scheduler-extra` |
| `SECLIP_RETRIES` | `--retries` |
| `SECLIP_LATENCY_WAIT` | `--latency-wait` |
| `SECLIP_MAX_JOBS_PER_SECOND` | `--max-jobs-per-sec` |
| `SECLIP_MAX_STATUS_CHECKS_PER_SECOND` | `--max-status-per-sec` |
| `SECLIP_SOFTWARE_CONFIG` | `--software` |
| `SECLIP_CONFIG` | `-c / --config` |
| `SECLIP_EXTRA_CONFIG` | `-l / --extra-config` |
| `SECLIP_RESOURCES_CONFIG` | project `resources.yaml` lookup |

Two further variables are internal but also read by direct Snakemake invocation: `SECLIP_TOOL_CLIPPER` (exported by run.sh from `software.yaml` `paths.clipper`) and `SECLIP_PYTHON` (interpreter used by the version-collection rule).

---

## 7. Output interpretation

### 7.1 Working-directory layout

```
workdir/
├── config.local.yaml        # optional: working-directory config overlay (auto-layered, not committed)
├── samples.csv              # sample table (filename free; pointed to by SampleListFile)
├── snakemake.logs.txt       # launcher/Snakemake main log
├── 1.rawdata/               # raw inputs (the only data directory outside results/):
│                            #   {sample}.fastq.gz / {sample}_R1.fq.gz (SE + 10N UMI)
└── results/                 # ALL derived artifacts (rename via results_dir in config)
    ├── 0.index/
    │   ├── genome_STARindex/    # STAR genome index (built from genome + gtf)
    │   └── repeats_STARindex/   # STAR sncRNA/repeats index (filter_repeats: true)
    ├── 2.cleandata/
    │   ├── {sample}_umi.fq.gz                  # after UMI extraction
    │   ├── {sample}_clean.fqTrTr.fq.gz         # after both cutadapt passes
    │   ├── {sample}_clean.fqTrTr.sorted.fq.gz  # alignment input
    │   ├── logs/               # umi_extract / cutadapt metrics files
    │   └── fastqc/             # {sample}_fastqc.html / .zip
    ├── 3.align/
    │   ├── repeats/            # {sample}_Unmapped.out.mate1 (filter_repeats: true)
    │   └── genome/             # {sample}_Aligned.out.bam, {sample}_Log.final.out (+ STAR siblings)
    ├── 4.rmdup/
    │   ├── {sample}.rmDupSo.bam(.bai)          # UMI-deduplicated, coordinate-sorted BAM
    │   ├── {sample}_readnum.txt                # mapped-read count
    │   └── {sample}_stats/                     # umi_tools dedup statistics
    ├── 5.callpeak/             # {sample}.pureclip.bed, {sample}.clipper.peakClusters.bed
    ├── 6.reproducible_peaks/   # optional: {condition}.consensus.bed (support in column 4; with
    │                           #   input_control conditions with inputs carry the in_input_background
    │                           #   flag in column 5), plus {condition}.input_background.bed and
    │                           #   {condition}.consensus.filtered.bed (filter_by_input)
    ├── 6.annotation/           # optional: {set}.annotation.tsv + _ref/ (genes.bed, exons.bed, genes.tsv)
    ├── 5.QC/
    │   ├── multiqc/multiqc_report.html
    │   ├── software_versions.yaml
    │   └── logs/               # QC rule logs (e.g. software_versions.log.txt)
    └── logs/                   # per-rule logs ({rule}/{sample}.log, star_index_*.log, multiqc.log)
```

### 7.2 Results quick reference

| Result | Path | Meaning |
|---|---|---|
| UMI-extracted reads | `results/2.cleandata/{sample}_umi.fq.gz` | raw read with the leading UMI removed by `umi_tools extract` |
| extract metrics | `results/2.cleandata/logs/{sample}_umi_extract.metrics` | how many reads matched the UMI pattern |
| trimmed reads | `results/2.cleandata/{sample}_clean.fqTrTr.fq.gz` | output of the strict second cutadapt pass |
| cutadapt report | `results/2.cleandata/logs/{sample}_cutadapt.metrics` | adapter-trimming statistics (aggregated into MultiQC) |
| sorted reads | `results/2.cleandata/{sample}_clean.fqTrTr.sorted.fq.gz` | name-sorted alignment input |
| FastQC | `results/2.cleandata/fastqc/{sample}_fastqc.html` | per-base quality of the trimmed reads |
| repeats-unmapped reads | `results/3.align/repeats/{sample}_Unmapped.out.mate1` | reads that did NOT match sncRNA/repeats |
| genome BAM | `results/3.align/genome/{sample}_Aligned.out.bam` | unique end-to-end genome alignment |
| STAR stats | `results/3.align/genome/{sample}_Log.final.out` | mapping-rate summary (aggregated into MultiQC) |
| dedup BAM | `results/4.rmdup/{sample}.rmDupSo.bam` (+ `.bai`) | UMI-collapsed, coordinate-sorted alignments |
| dedup stats | `results/4.rmdup/{sample}_stats/{sample}_edit_distance.tsv` | umi_tools dedup statistics (aggregated into MultiQC) |
| read count | `results/4.rmdup/{sample}_readnum.txt` | mapped reads after dedup (one integer) |
| PureCLIP peaks | `results/5.callpeak/{sample}.pureclip.bed` | crosslink sites, BED6 (chromosome, start, end, site name, crosslink-site score, strand) |
| CLIPper peaks | `results/5.callpeak/{sample}.clipper.peakClusters.bed` | peak clusters (only with a configured CLIPper) |
| Consensus peaks (`reproducible_peaks.enabled`) | `results/6.reproducible_peaks/{condition}.consensus.bed` | per-condition cross-sample consensus; BED4 with the replicate support count in column 4 (§5.5); with `input_control`, conditions with input samples carry the binary `in_input_background` flag in column 5 |
| Input background (`reproducible_peaks.input_control`) | `results/6.reproducible_peaks/{condition}.input_background.bed` | union of the condition's input-control PureCLIP beds, BED4 (column 4 = number of inputs covering the feature) |
| Filtered consensus (`reproducible_peaks.filter_by_input`) | `results/6.reproducible_peaks/{condition}.consensus.filtered.bed` | consensus minus the input-background-flagged sites (BED4) |
| Peak annotation (`annotate_peaks.enabled`) | `results/6.annotation/{set}.annotation.tsv` | nearest gene / distance / biotype per peak set, one row per peak (§5.5) |
| Annotation reference | `results/6.annotation/_ref/` | GTF-derived `genes.bed` / `exons.bed` / `genes.tsv` behind the annotation |
| QC report | `results/5.QC/multiqc/multiqc_report.html` | everything above in one HTML |
| version record | `results/5.QC/software_versions.yaml` | tool versions actually resolved for the run (incl. Snakemake) |
| per-rule logs | `results/logs/` | one log per rule/sample |

### 7.3 Reading the QC (no Bismark modules here)

This workflow intentionally ships no Bismark/RNA-seq-specific QC modules; judge a run by this chain:

1. **UMI extraction** (`*_umi_extract.metrics`, MultiQC umi_tools module): the fraction of raw reads carrying the expected UMI pattern should be near 100% for the legacy library design; a low value means wrong `umi.pattern` or a non-UMI library;
2. **Trimming** (`*_cutadapt.metrics`, MultiQC cutadapt module): the adapter-contaminated fraction reflects library quality; the final length distribution must clear `cutadapt.min_len` (reads shorter than 18 bp by default are discarded);
3. **FastQC** (`*_fastqc.html`): per-base quality after trimming — residual adapter signal or a quality collapse at the 3' end points at wrong adapter sequences;
4. **Alignment** (`*_Log.final.out`, MultiQC STAR module): input reads vs uniquely mapped reads; with `align_multimap_nmax: 1` the unique fraction is the usable signal;
5. **Deduplication** (`*_edit_distance.tsv`, MultiQC umi_tools module) and **`{sample}_readnum.txt`**: how many unique molecules survived; the ratio of dedup output to alignment input is the effective library complexity;
6. **Peaks** (`results/5.callpeak/`): compare per-sample BED sizes; wildly different peak counts across replicates of the same IP suggest depth or quality problems. With `reproducible_peaks.enabled` the consensus support column (§5.5) adds the cross-replicate view directly.

Everything above lands in `results/5.QC/multiqc/multiqc_report.html` (title "seCLIP QC Summary"); the exact tool versions used are recorded in `results/5.QC/software_versions.yaml`.

---

## 8. FAQ

**Q1: Can I run paired-end (PE) data?**
No. The workflow is single-end only: the sample-table parser accepts exactly one `sample_id` column (plus the optional `condition`/`role` grouping columns, §3) and the FASTQ resolver looks only for `{sample}[_R1]{.fastq,.fq}.gz` files. The 10-base UMI sits at the start of the single read (`umi.pattern`, configurable if your library layout differs). PE support would need a new sample-table design and is not scheduled.

**Q2: CLIPper is not installed on my server — is the workflow unusable?**
No. CLIPper is optional by design: leave `paths.clipper` empty in software.yaml and the workflow auto-skips the `callpeak_clipper` rule with a `[config warning]` line, keeping PureCLIP peaks. To enable it, install CLIPper anywhere and put its absolute path in `paths.clipper` (a project software.yaml wins over the repository default), then set `callpeak.clipper_species` to the value your CLIPper build expects (default `GRCh38_v40`).

**Q3: Why do the `genome_sa_index_nbases` / `repeats_sa_index_nbases` settings matter?**
STAR's suffix-array seed length must satisfy roughly `min(14, log2(genomeSize)/2 - 1)`: the human default 14 is right for GRCh38 but far too large for small references — index building crashes or hangs (that is why the synthetic test config uses 4/2). The repeats index defaults to 12 because the sncRNA/repeats FASTA is small. When you point the workflow at a small genome, shrink `star.genome_sa_index_nbases` accordingly (the FAQ-level rule of thumb: halve it for every ~4x of genome-size reduction and check the STAR manual).

**Q4: How do I resume after an interruption? What about lock errors?**
Snakemake skips completed steps by output timestamps — **rerun the same command to resume** (profiles pin keep-going and rerun-incomplete). If you hit a working-directory lock error (`Directory cannot be locked`, usually after a force-killed job), run `bash run.sh -P . --unlock` and rerun.

**Q5: Can I add samples to a project that already ran?**
Yes. Append the new ids to the sample table, drop their FASTQ files into `1.rawdata/`, and rerun the same command: finished outputs are reused and only the new per-sample jobs (plus a fresh MultiQC/versions run) execute. Two caveats: ids must obey the §3.1 naming rules, and the DAG rebuilds from scratch on the first run after the change (Snakemake re-checks every target), which is normal.

**Q6: What does `[config warning] ... is still a placeholder: /path/to/...` mean?**
Parse-time validation warns when `genome`/`gtf`/`repeats_fa` still hold the species-preset placeholder paths or are missing entirely. Dry-runs work anyway; a real run fails at the index/alignment rules. Fill real paths in the project config (explicit keys override the preset) before launching.

**Q7: A job was killed by the cluster for exceeding walltime/memory. What now?**
Three fixes, in order of preference: (1) targeted per-rule override in a project `resources.yaml` (§4.4); (2) global `--runtime` / `--memory` overrides (§6.3); (3) `--retries N` so sporadic failures retry automatically. To locate the rule: find the failed rule name in `snakemake.logs.txt`, then read `results/logs/<rule>/<sample>.log` for details.

**Q8: How do I validate a new deployment or a workflow change?**

```bash
make check                          # bash -n syntax checks (no snakemake needed)
make lint                           # static suite (missing optional tools are skipped)
make unit                           # pytest suite for the annotation scripts
bash tests/run_test.sh              # synthetic-data dry-run regression (needs snakemake + python3/PyYAML)
bash tests/run_test.sh --consensus  # dry-run with the optional v0.2 stages enabled
bash tests/run_test.sh --input-control  # dry-run with the input-control consensus scenario
bash tests/run_test.sh --real-run   # end-to-end run + output assertions (needs the full analysis environment)
```

Test data is generated by `tests/make_testdata.py` with a fixed seed (2 x 20 kb chromosomes, 2 samples; 4 with `--with-inputs`); the default dry-run DAG is 23 jobs (CI passes `--reads 2000`), the `--consensus` scenario dry-runs 28 jobs, and the `--input-control` scenario 45. Before changing workflow code, read the documentation-sync checklist in [CONTRIBUTING](../CONTRIBUTING.md).

**Q9: How do I get reproducible peaks across replicates, and what do the annotation columns mean?**
Enable the two optional stages (§5.5): extend the sample table to `sample_id,condition,role` (§3), then set `reproducible_peaks.enabled: true` (with `min_replicates`) and `annotate_peaks.enabled: true` in the project config. Each condition gets `results/6.reproducible_peaks/{condition}.consensus.bed` (BED4, column 4 = number of replicates carrying the site) and every peak set gets `results/6.annotation/{set}.annotation.tsv` (`feature_class` = exon/gene/intergenic, `distance` = signed distance to the nearest gene with 0 = overlapping). Both stages stay out of the DAG while disabled, so existing projects are unaffected. With `input_control: true` the consensus of conditions with input samples additionally carries the binary `in_input_background` flag (column 5), and `filter_by_input: true` writes a filtered BED without those background sites (§5.5).
