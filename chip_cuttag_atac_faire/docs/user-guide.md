# ChIP/CUT&Tag/ATAC/FAIRE workflow user guide

> Updated: 2026-09-08 (v0.5: replicate-aware peak stage with IDR/consensus, TSS enrichment, organelle fraction, blacklist filtering — all default-off; environment + run.sh ops CLI as of v0.4.0)
> Intended audience: analysts running this workflow on their cluster/server
> Since v0.4.0 the workflow no longer creates per-rule Conda environments: the runtime environment is created from `workflow/environment.yaml`, or reuses the server's existing environment and R libraries via `config/software.yaml`; cluster scheduling resources are declared per rule in `config/resources.yaml` and support project-level overrides. All derived outputs are consolidated under the project's `results/` directory (`results_dir` in config).

---

## 1. Installation and environment setup

### 1.1 Prerequisites

| Dependency | Notes |
|---|---|
| Linux (PBS/SGE/SLURM cluster or a single machine) | the scheduler is auto-detected by `run.sh` or set via `--profile` |
| Snakemake | reference version **7.32.4** (pinned in `workflow/environment.yaml`); see the version matrix in §1.4 |
| Python 3 + PyYAML | the bootstrap dependency `run.sh` needs to resolve `software.yaml` and run preflight checks; must be on the main PATH |
| Analysis tools | fastqc / trim_galore / multiqc / bowtie2 / samtools / picard / macs2 / bedtools / deeptools / phantompeakqualtools / R (ChIPseeker, GenomicFeatures) |

Launcher bootstrap order: `run.sh` first uses the `python3` on the system PATH (with PyYAML) to resolve `config/software.yaml`, then injects the main environment/tools into the current process — so even on the conda_prefix reuse route, the login node PATH must have `python3` and `snakemake` (the latter supplied by the main environment after resolution).

### 1.2 Three ways to set up the environment

**Path 1: fresh server — create the environment from the all-in-one template**

```bash
mamba env create -f workflow/environment.yaml   # environment name chip-cuttag-atac-faire
conda activate chip-cuttag-atac-faire
```

This command is executed explicitly by the user; Snakemake never creates or modifies software environments on its own. Pinned versions in the template: snakemake-minimal 7.32.4, fastqc 0.12.1, trim-galore 0.6.10, multiqc 1.14, bowtie2 2.5.1, samtools 1.17, picard 3.0.0, macs2 2.2.7.1, bedtools 2.31.0, deeptools 3.5.1, phantompeakqualtools 1.2.2, r-base 4.3, bioconductor-chipseeker 1.38.0, bioconductor-genomicfeatures 1.54.0.

**Path 2: reuse an existing conda environment on the server**

Copy `config/software.yaml` into the project directory (or edit the repository default) and set:

```yaml
environment:
  type: conda
  conda_prefix: "/share/conda/envs/chip-cuttag-atac-faire"   # recommended on HPC
  # conda_name: "chip-cuttag-atac-faire"                     # or resolve by environment name
  strict: true
```

No `conda activate` needed: after resolution, `run.sh` injects that prefix's `bin` into PATH. An environment created via path 1 can also be connected this way (`environment.yaml` then only serves to "create the environment"; the runtime always goes through `software.yaml`).

**Path 3: system mode + reuse the server's existing R libraries**

When all tools come from the system PATH and R is installed separately:

```yaml
environment:
  type: system
  strict: true

r:
  rscript: "/opt/R/4.3/bin/Rscript"   # or "Rscript" to use PATH
  version: "4.3"                       # optional check
  version_check: major_minor           # major_minor | exact | warn | off
  lib_paths:                           # reuse existing server libraries such as ChIPseeker
    - "/share/Rlibs/4.3"
  lib_mode: prepend                    # prepend | append | replace
```

`r.lib_paths` reaches all R steps via `R_LIBS_USER`; `prepend` does not shadow R's own base/site libraries. Skip this section when the conda environment already ships ChIPseeker/GenomicFeatures.

### 1.3 Preflight checks

```bash
# Check that all executables, R, and R packages resolve (then exit)
bash run.sh -P /path/to/workdir --check-software

# Check only the configured Rscript, R version, and R packages
bash run.sh -P /path/to/workdir --check-r

# Parse the workflow + all configs (sample-table/config validation runs at parse time),
# list the rules, then exit without starting jobs
bash run.sh -P /path/to/workdir -c /path/to/workdir/config.yaml --validate-only
```

Notes:

- `--check-software` / `--check-r` require `-P` to name the project directory (created automatically when missing);
- before a real run, `run.sh` executes a full preflight **automatically**; `--skip-software-check` skips it, and `-n` (dry-run) skips it automatically;
- `software.yaml` lookup order: `--software FILE` > `<project dir>/software.yaml` > repository `config/software.yaml`. The `tools:` section only needs to override tools that differ (the deeptools suite is represented by the `bamCoverage` binary and the SPP entry point is `run_spp.R`; both are pre-filled with the bioconda package command names).

### 1.4 Snakemake version matrix

| snakemake version | support | notes |
|---|---|---|
| **7.32.4** | ✅ reference version | pinned in `workflow/environment.yaml`; cluster profiles (pbs/sge/slurm) target the 7.x classic `--cluster` interface |
| **8.x** | ⚠️ cluster semantics unverified | parsing/lint verified under the unified environment; 8.x moved to the executor plugin system (`--cluster` semantics removed), the launcher prints a warning at startup — test on a small scale before real cluster runs |
| <7 or ≥9 | ⛔ unverified | test before use |

---

## 2. Quick start

Five steps to a minimal project (using the rice example as the reference; switch species by changing the reference files):

**Step 1: create the working directory**

```bash
mkdir -p ~/work/chip_demo/1.rawdata
```

The working directory (data) and the analysis code directory (this repository) are separate; `run.sh` generates all results inside the working directory.

**Step 2: add the raw data**

Place the paired raw fastq files in `1.rawdata/`. Each sample's pair is auto-detected by priority: `{sample}_1.fastq.gz` + `{sample}_2.fastq.gz`, `{sample}_1.fq.gz` + `{sample}_2.fq.gz`, `{sample}_R1.fastq.gz` + `{sample}_R2.fastq.gz`, or `{sample}_R1.fq.gz` + `{sample}_R2.fq.gz` (the first complete pair wins; a sample without any complete pair fails validation listing all supported variants). The optional batch rename below is only needed to normalize mixed deliveries:

```bash
cp my_*.fastq.gz ~/work/chip_demo/1.rawdata/
bash run.sh -P ~/work/chip_demo -r        # -r only takes effect when 1.rawdata/ exists; perl rename syntax
```

**Step 3: write the sample table**

```bash
cp config/samples.csv ~/work/chip_demo/sample_info.csv
```

Edit the file with your own samples (6-column schema and validation rules in §3). Note that the `grouplist` default in `config.yaml` is `config/samples.csv` (the repository template); **relative paths resolve against the working directory first** — the sample-table filename is free, just point `grouplist` at it.

**Step 4: write the project config**

```bash
cp config/config.yaml ~/work/chip_demo/config.yaml
```

Set the reference genome file set (`genome_fa` / `gtf` / `bed` / `chromsize`) and `genome_size` (MACS2 effective genome size) for your species — or simply keep `species: "osa"` / switch to `"hsa"` and let the `config/species.yaml` preset supply the reference keys (explicit keys in your config win over the preset; see §4.2). `run.sh` auto-detects `<project dir>/config.yaml` when `-c` is omitted, so placing it in the working directory is enough. Full key reference in §4.

**Step 5: preflight and launch**

```bash
cd /path/to/repo
bash run.sh -P ~/work/chip_demo --check-software   # environment preflight
bash run.sh -P ~/work/chip_demo -n                 # dry-run: builds the DAG and prints planned jobs without executing
bash run.sh -P ~/work/chip_demo -j 10              # real run (local / auto-detected scheduler)
```

After a successful launch the main log is `~/work/chip_demo/snakemake.logs.txt` (change with `--log FILE`), and per-step logs live under `results/logs/`. See §7 for the results layout.

---

## 3. Sample table in detail

The sample table is a 6-column CSV (template `config/samples.csv`; a mixed-assay example ready to copy):

```csv
sample_id,role,group,seqtype,layout,peak_type
myc,treat,myc_vs_IgG,chip,PE,narrow
IgG,control,myc_vs_IgG,chip,PE,narrow
H3K27me3_rep1,treat,H3K27me3_vs_IgG,cuttag,PE,broad
IgG_cuta,control,H3K27me3_vs_IgG,cuttag,PE,broad
atac_leaf_1,treat,atac_leaf,atac,PE,none
atac_leaf_2,treat,atac_leaf,atac,PE,none
faire_root,treat,faire_root_vs_Input,faire,PE,none
Input_faire,control,faire_root_vs_Input,faire,PE,none
```

### 3.1 Column definitions

| Column | Values | Notes |
|---|---|---|
| `sample_id` | non-empty | sample name; must match the fastq prefix in `1.rawdata/` |
| `role` | `treat` / `control` | treatment / control; **controls (IgG/Input) must be listed too** — they are aligned as well |
| `group` | non-empty | peak-calling group name; a group may contain multiple treats and multiple controls (MACS2 multi-file `-c` input) |
| `seqtype` | `chip` / `cuttag` / `atac` / `faire` | assay type; mixing rows is what makes a mixed-assay project |
| `layout` | `PE` | only paired-end is supported in this version |
| `peak_type` | `narrow` / `broad` / `none` | treat rows of chip/cuttag must be `narrow` or `broad`; atac/faire are fixed to `none` (the peak type is decided by the workflow) |
| `condition` | optional free label | differential binding (§5.6): labels the contrast factor; a group's treats must share one value; falls back to the group name when absent |
| `batch` | optional free label | differential blocking factor (sequencing batch); ignored unless `diffbind.batch_correction` is on and it has ≥ 2 levels |

### 3.2 Validation rules (row by row at parse time, errors carry line numbers)

Validation runs centrally at Snakemake parse time (`load_sample_table` in `workflow/rules/common.smk`); any violation aborts with the offending line number:

1. All 6 columns must be present; missing columns error out immediately and list the missing names;
2. `sample_id` / `group` are non-empty and allow only **alphanumeric plus `. _ -`**: must start alphanumeric, must not start with `-`, and **must not contain consecutive underscores `__`** (`__` is the FRiP output filename separator; commas/spaces would break peak-list joining and shell expansion);
3. `role` accepts only `treat`/`control`; `seqtype` only the four assays; `layout` only `PE`; `peak_type` only `narrow`/`broad`/`none`;
4. atac/faire rows must have `peak_type` `none`; **treat** rows of chip/cuttag must be `narrow` or `broad`;
5. The same `sample_id` must not appear in groups with different `seqtype` values;
6. Within one `group`, every row must agree on `seqtype` / `peak_type` / `layout`;
7. Every group contains at least 1 `role=treat` sample; the table has at least 1 data row.

### 3.3 Mixed-assay notes

- The four assays can be mixed row by row with no launcher options; the `seqtype` column drives the whole workflow's routing (dedup strategy, peak-calling rules; see §5.1);
- Can chip and cuttag share one `group` name — **no**: `seqtype` must be consistent within a group (rule 6); use different groups for different assays;
- The sample-table path comes from `grouplist` in the config: absolute path > working-directory-relative > repository-relative. If none of the three resolves, parsing **errors out** (no silent fallback to the development example since v0.4.0); see `example/samples.csv` for a real-project example.

---

## 4. Configuration in detail

### 4.1 Config chain and layering

`run.sh` assembles the `--configfile` chain in this order; **later files override earlier keys**:

1. Repository default `config/config.yaml` (always first in the chain, the fallback);
2. Project config: `-c FILE` explicitly; when omitted, `<project dir>/config.yaml` is auto-detected;
3. Overlay config: `-l FILE` explicitly; when omitted, `<project dir>/config.local.yaml` is auto-detected (an INFO line is printed when found).

Typical usage: the project config holds only project-specific differences; machine-specific differences (e.g. thread caps) go into `config.local.yaml` (kept out of version control); the repository `config/config.template.yaml` is the annotated override template — copy it to `config.local.yaml` and change what differs (every key is optional; unlisted keys keep the defaults).

Scheduler resources are layered separately: see §4.3 for `config/resources.yaml` and the project-level `resources.yaml`.

### 4.2 All config.yaml keys

| Key | Default | Notes |
|---|---|---|
| `species` | `"osa"` | selects a species preset from `config/species.yaml` (`osa` = rice IRGSP-1.0, `hsa` = human GRCh38/GENCODE) |
| `results_dir` | `"results"` | root directory for all derived outputs (raw inputs `1.rawdata/` stay at the working-directory root) |
| `genome_fa` | from the species preset | bowtie2 index input (**must change per species unless a preset fits**) |
| `gtf` | from the species preset | for peak annotation (ChIPseeker makeTxDbFromGFF) |
| `bed` | from the species preset | for deeptools gene-region signal |
| `chromsize` | from the species preset | for bigWig generation (chromosome length table) |
| `genome_size` | from the species preset (`osa` = `"3.7e8"`) | MACS2 effective genome size (rice ~3.7e8, human 2.7e9, mouse 1.87e9; **must change per species**) |
| `grouplist` | `config/samples.csv` | sample-table path (resolves relative to the working directory first, see §3.3) |
| `threads` | `12` | legacy global thread cap applied to every rule (integer, unquoted); per-rule requests live in `config/resources.yaml` |
| `trim.quality` | `25` | trim_galore `-q` terminal low-quality base threshold |
| `trim.stringency` | `3` | minimum adapter overlap (bp) |
| `trim.error_rate` | `0.1` | `-e` max adapter mismatch rate |
| `trim.extra` | `""` | extra arguments, e.g. `"--clip_r1 5 --clip_r2 5"` |
| `bowtie2_extra` | `--end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700` | bowtie2 argument string |
| `min_mapq` | `30` | MAPQ filter threshold (common ENCODE value) |
| `dedup.chip` / `dedup.atac` / `dedup.faire` | `true` | picard MarkDuplicates deduplication switch (per assay) |
| `dedup.cuttag` | `false` | CUT&Tag keeps PCR duplicates (Kaya-Okur et al. 2019) |
| `peak.keepdup` | `all` | MACS2 keeps all duplicates (deduplication is decided upstream) |
| `peak.qvalue` | `0.05` | narrow peak q-value cutoff |
| `peak.broad_cutoff` | `0.05` | broad peak cutoff |
| `peak.replicate.enabled` | `false` | replicate-aware peak stage: per-replicate calling + IDR (narrow) / overlap consensus (broad); see §5.3 |
| `peak.replicate.qvalue` | `0.01` | relaxed per-replicate narrow cutoff feeding IDR (ENCODE style) |
| `peak.replicate.idr_threshold` | `0.05` | global IDR cutoff |
| `peak.replicate.idr_rank` | `"p.value"` | narrowPeak rank column handed to idr: `p.value` or `signal.value` |
| `peak.replicate.consensus_min_replicates` | `2` | support cutoff (broad consensus and >2-replicate IDR unions) |
| `peak.replicate.frip_on` | `"pooled"` | peak set FRiP is computed against: `pooled` or `consensus` (requires the stage enabled) |
| `peak.atac.mode` | `bampe` | `bampe` = ENCODE ATAC v2 recipe (pile up real fragment lengths); `shifted` = classic Tn5 offset recipe (`--nomodel --shift -100 --extsize 200`) |
| `peak.atac.shift` / `peak.atac.extsize` | `-100` / `200` | effective in `shifted` mode only |
| `peak.bigwig_measure` | `"FE"` | group signal-track measure: `FE` (fold enrichment) or `logFE`; the filename stays `{group}_FE.bw` |
| `bigwig.per_sample` | `false` | per-sample normalized coverage bigWigs under `results/4.peak/samples/` (browser-level replicate comparison) |
| `bigwig.normalize` / `bigwig.bin` | `"RPGC"` / `25` | bamCoverage normalization (RPGC uses `genome_size`) and bin size |
| `motif.enabled` | `false` | HOMER motif enrichment on the final peak sets (§5.5); needs an external HOMER install |
| `motif.homer_genome` | `""` | required when the stage is on: HOMER genome tag (`hg38`, `mm10`, …) or `custom:/path/to/genome` |
| `motif.size` / `motif.background` / `motif.extra` | `"given"` / `""` / `""` | `-size`, optional `-bg` BED, extra findMotifsGenome.pl arguments |
| `diffbind.enabled` | `false` | DiffBind differential binding between sample-table groups (§5.6) |
| `diffbind.contrasts` | `[]` | list of `[groupA, groupB]` pairs; each arm needs ≥ 2 treats |
| `diffbind.analysis` | `"DESeq2"` | backend: `DESeq2` or `edgeR` |
| `diffbind.summit_flank` | `250` | count regions = summit ± this many bp (`0` = full peak regions) |
| `diffbind.use_controls` / `fdr` / `foldchange` / `batch_correction` | `false` / `0.05` / `1.0` / `true` | attach single controls as background; significant-table cutoffs; blocking on the optional `batch` column |
| `blacklist` | `""` | optional BED of artifact regions; empty disables. Filtered peak copies feed FRiP/annotation (§5.4) |
| `region_flank` | `3000` | peak-annotation flank distance and deeptools signal window up/downstream length (bp); one key controls both |
| `qc.tss` | `false` | TSS enrichment for atac/faire treat samples (per-sample CPM coverage ±2kb around TSS; §5.2) |
| `qc.organelle` | `false` | organelle (chloroplast/mitochondrion) mapped-read fraction from idxstats (§5.2) |
| `qc.organelle_patterns` | `[chrc, chrm, pt, mt, pltd, chloroplast, mitochondr, plastid]` | contig-name patterns; ≤3-char patterns match contig names exactly (case-insensitive), longer ones as substrings |

Reference keys (`genome_fa`/`gtf`/`bed`/`chromsize`/`genome_size`) that the project config leaves unset fall back to the selected species preset; an explicit key in the project config wins over the preset.

### 4.3 Per-rule scheduler resources (`config/resources.yaml`)

Scheduler resources live in their own file, `config/resources.yaml` — the single place to tune cluster requests. Each rule declares `threads`, `mem_mb`, and `runtime_min`; `runtime_sec` (used by PBS walltime and SGE h_rt) is derived automatically as `runtime_min × 60` and cannot be set independently.

Override without touching the repository: copy the file into your project directory as `resources.yaml` (next to `config.yaml`) and edit — `run.sh` auto-detects it. Lookup order: the `CHIP_RESOURCES_CONFIG` environment variable > `<project dir>/resources.yaml` > repository `config/resources.yaml`. Rules not listed in an override file keep the repository defaults, and the legacy top-level `threads` in config.yaml still acts as a global cap on every rule's threads.

Example override (`<project dir>/resources.yaml`):

```yaml
resources:
  bowtie2_mapping:
    mem_mb: 32768
    runtime_min: 480
  callpeak_narrow:
    runtime_min: 360
```

Built-in defaults per rule (identical to `config/resources.yaml`; analysis parameters stay in config.yaml):

| Rule | threads | mem_mb | runtime_min |
|---|---:|---:|---:|
| `software_versions` | 1 | 1024 | 10 |
| `trim_adapter` | 4 | 4096 | 60 |
| `fastqc` | 2 | 2048 | 30 |
| `multiqc` | 1 | 4096 | 30 |
| `bowtie2_index` | 8 | 8192 | 120 |
| `bowtie2_mapping` | 8 | 16384 | 240 |
| `dedup` | 2 | 8192 | 120 |
| `callpeak_narrow` | 1 | 8192 | 180 |
| `callpeak_broad` | 1 | 8192 | 180 |
| `callpeak_atac` | 1 | 8192 | 180 |
| `bigwig` | 1 | 4096 | 60 |
| `peak_annotation` | 1 | 8192 | 120 |
| `frip` | 1 | 4096 | 60 |
| `frip_summary` | 1 | 1024 | 10 |
| `deeptools_multibamsummary` | 2 | 8192 | 120 |
| `deeptools_correlation` | 1 | 4096 | 30 |
| `deeptools_pca` | 1 | 4096 | 30 |
| `deeptools_fingerprint` | 2 | 8192 | 60 |
| `deeptools_fragmentsize` | 2 | 8192 | 60 |
| `deeptools_profile` | 2 | 8192 | 120 |
| `spp_crosscorr` | 2 | 8192 | 180 |
| `spp_summary` | 1 | 1024 | 10 |

### 4.4 Validation behavior

Config validation runs centrally at workflow parse time (`validate_config` in `workflow/rules/common.smk`): required keys, sub-keys, types, and value ranges are **aggregated into a single report** (e.g. a string `dedup.cuttag` is listed alongside any other errors). Missing reference files only **print a warning and do not abort** — dry-run/lint often run on machines without the reference files; in a real run the missing reference will fail the corresponding rule, so confirm each warning before launching.

---

## 5. Run modes and assay routing

### 5.1 seqtype automatic routing

The four assays are driven by the sample table's `seqtype` column; one project can mix them, no launcher options needed:

| seqtype | dedup (`dedup.<assay>`) | peak-calling rule | peak format |
|---|---|---|---|
| `chip` | picard deduplication (on by default) | `callpeak_narrow` (peak_type=narrow) or `callpeak_broad` (broad) | narrowPeak / broadPeak |
| `cuttag` | **no deduplication** (off by default, PCR duplicates kept) | same as chip, routed per group peak type | narrowPeak / broadPeak |
| `atac` | picard deduplication (on by default) | `callpeak_atac` (`peak.atac.mode` = bampe/shifted) | narrowPeak |
| `faire` | picard deduplication (on by default) | `callpeak_atac` (same rule as atac) | narrowPeak |

Peak calling runs in parallel per `group`; a group without a control automatically omits MACS2 `-c` (falling back to local lambda estimation). Each group's signal track `{group}_FE.bw` is produced by bdgcmp → bedClip → bedtools sort -g → bedGraphToBigWig (chromosome order consistent with chrom.sizes).

### 5.2 The QC switches

The `qc:` section of config.yaml controls the optional QC modules (their rule sets are conditionally included by the Snakefile; when off they stay out of the DAG):

| Switch | Default | Outputs |
|---|---|---|
| `qc.frip` | `true` | `results/5.QC/frip/FRiP_summary.tsv` (also injected into the MultiQC report) |
| `qc.deeptools` | `true` | `results/5.QC/deeptools/`: correlation heatmaps, PCA, fingerprint plots, fragment-size distribution, gene-region signal profiles |
| `qc.nsc_rsc` | `false` | `results/5.QC/spp/NSC_RSC_mqc.tsv` (SPP cross-correlation, slow; also injected into MultiQC) |
| `qc.tss` | `false` | `results/5.QC/tss/`: per atac/faire treat sample a TSS profile plot + enrichment score (`{sample}_TSSE.txt`), plus `TSSE_summary.tsv` (injected into MultiQC). The score is the max of the ±2kb profile normalized by the outer-flank baseline (ENCODE-flavored definition on CPM coverage) — healthy ATAC libraries show a clear TSS spike |
| `qc.organelle` | `false` | `results/5.QC/organelle/Organelle_summary.tsv` (injected into MultiQC): per-sample chloroplast/mitochondrial mapped-read fraction from `samtools idxstats`. Plant ATAC libraries frequently lose a large fraction of reads to the chloroplast; the number is diagnostic for library quality and for whether the reference should be nuclear-only |

Regardless of the switches, `results/5.QC/software_versions.yaml` (record of the tool versions actually used, including the Snakemake version) and the MultiQC summary report are always generated; FastQC, bowtie2 alignment stats, and picard dedup metrics are pulled into MultiQC automatically.

### 5.3 Replicate-aware peak stage (v0.5, `peak.replicate`)

By default the workflow calls peaks once per group on the pooled replicates (`-t treat1,treat2`). With `peak.replicate.enabled: true` the replicate structure is additionally analyzed:

- **Per-replicate calling**: every treat sample gets its own MACS2 call (`results/4.peak/replicates/{group}/{sample}_peaks.*`) at the relaxed narrow cutoff `peak.replicate.qvalue` (default 0.01), always against the group's pooled control. Broad groups use the same `broad_cutoff` as the pooled call.
- **IDR (narrow groups, ≥2 treats)**: all replicate pairs run through the classic `idr` tool (`--rank p.value --idr-threshold 0.05` by default; both configurable). A group with exactly 2 treats uses the single pair result as its final reproducible peak set (`results/4.peak/{group}_IDR_peaks.narrowPeak`). With >2 treats the union of all pairwise IDR results is kept where the pairwise support ≥ `consensus_min_replicates` — a deliberate simplification of ENCODE's rescue/self-consistency scheme, chosen for interpretability; the per-peak support is in `{group}_IDR_support.bed`.
- **Overlap consensus (broad groups, ≥2 treats)**: IDR is not applicable to broad peaks; the per-replicate broadPeak files go through `bedtools multiinter` and intervals carried by ≥ `consensus_min_replicates` replicates form `results/4.peak/{group}_consensus_peaks.broadPeak` (+ support bed).
- **Summary**: `results/5.QC/replicate_peaks/Replicate_summary.tsv` lists per group the replicate peak counts, the final (IDR/consensus) count, and the retained fraction; it is injected into MultiQC.
- **Downstream wiring**: `peak.replicate.frip_on: consensus` computes FRiP against the reproducible set (default `pooled` keeps today's semantics; single-treat groups always fall back to pooled), and the ChIPseeker annotation covers the final reproducible set. The pooled peak files keep their names and remain in `results/4.peak/`.

**Environment note**: the `idr` tool (Liu et al., 2.0.4.x) is python2-based and deliberately not part of the conda template. Install it separately (e.g. `conda create -n idr -c bioconda idr=2.0.4`) and either put `idr` on PATH or point the software.yaml `paths:` entry at the binary (exported to the rules as `CHIP_IDR`). It is only required when the stage is enabled and narrow groups with ≥2 treats exist; dry-runs never execute it.

### 5.4 Blacklist filtering (v0.5, `blacklist`)

Set the top-level `blacklist` key to a BED file of known artifact regions and the workflow writes filtered copies (bedtools `intersect -v`) of the pooled and final peak sets into `results/4.peak/blacklist_filtered/`, with a before/after count table in `results/5.QC/blacklist/blacklist_summary.tsv`. FRiP and peak annotation read the filtered copies automatically. Filtering happens at the peak level only (BAMs are untouched). No ENCODE blacklist exists for rice — build or borrow one appropriate for your genome, or leave the key empty (default).

### 5.5 HOMER motif enrichment (v0.5, `motif`)

With `motif.enabled: true`, every group's final peak set (the same deliverable annotation uses — IDR/consensus when the replicate stage is on, blacklist-filtered when a blacklist is set) goes through `findMotifsGenome.pl`; results land in `results/6.motif/{group}/` (de novo + known motif tables and logos). Requirements and knobs:

- `motif.homer_genome` is mandatory: a HOMER genome tag installed via `configureHomer`, or a `custom:/path/to/genome` directory. Rice has no stock HOMER genome — configure one for your assembly.
- HOMER is an external distribution (not in the conda template): the `findMotifsGenome.pl` entry point resolves via the software.yaml `paths:` section (`homer_findmotifs` → `CHIP_HOMER_FINDMOTIFS`) or from PATH.
- `motif.size` (`given` = peak widths), an optional matched `motif.background` BED, and free-form `motif.extra` (e.g. `"-len 8,10,12 -nmotifs 12"`) cover the common recipes.
- De novo discovery is slow (hours on large peak sets); tune the `motif_enrichment` entry in `config/resources.yaml` before queueing.

### 5.6 Differential binding (v0.5, `diffbind`)

Differential enrichment between two conditions, DiffBind (DESeq2 or edgeR backend), driven by explicit contrasts:

```yaml
diffbind:
  enabled: true
  contrasts: [["H3K27ac_WT_vs_IgG", "H3K27ac_mut_vs_IgG"]]   # two sample-table groups
```

Each contrast arm must carry ≥ 2 treat replicates (parse-time error otherwise). The sample table's optional `condition` column labels the two factor levels (falling back to the group names); the optional `batch` column becomes a blocking factor when `diffbind.batch_correction` is on and has ≥ 2 levels. Per contrast, `results/6.diffbind/{A}__vs__{B}/` receives the generated sample sheet, `DB_results.tsv` (every consensus region with Fold/FDR/p), `DB_significant.tsv` (passing `fdr` and `|Fold| >= foldchange`), MA/volcano/PCA plots, and sessionInfo.

Two design notes: DiffBind counts reads over the consensus peak set derived from the **per-sample peak files** — with `peak.replicate.enabled: true` those are the per-replicate calls (recommended); without it every sample of a group maps to the identical pooled set, which still runs but loses replicate-level peak structure. `diffbind.summit_flank: 250` recenters counting on summits (narrow peaks); set `0` to count full peak regions (usually better for broad marks). `diffbind.use_controls: true` attaches a group's exactly-one control as the DiffBind background sample.

**Environment note**: `bioconductor-diffbind` is part of the conda template (pulls DESeq2/edgeR transitively). When reusing server R libraries instead, add DiffBind (e.g. `BiocManager::install("DiffBind")`) to the `r.lib_paths` libraries.

### 5.7 Checks and dry-run

```bash
bash run.sh -P . -n                    # dry-run: builds the DAG and prints the jobs it would run, without executing
bash run.sh -P . --validate-only       # parses workflow + config, lists rules, then exits
bash run.sh -P . --check-software      # tools + R + R packages preflight
bash run.sh -P . --check-r             # R-side preflight only
```

The dry-run automatically skips the software preflight (no tools executed; needs only snakemake + python3/PyYAML, not the full analysis environment). Sample-table and config validation happen at parse time, so the dry-run and `--validate-only` catch the same errors listed in §3.2 / §4.4.

### 5.8 Resuming and Snakemake passthrough

- Resuming: Snakemake skips completed steps based on output timestamps; after an interruption, **simply rerun the same command**;
- profiles pin `keep-going: true` and `rerun-incomplete: true`; `latency-wait` defaults to 90 (default/pbs) or 60 (sge/slurm), overridable via `--latency-wait SEC`;
- a stale working-directory lock (e.g. a job was force-killed) is cleared with `bash run.sh -P . --unlock`;
- arguments after `--` pass through to Snakemake untouched, e.g. `bash run.sh -P . -- --rerun-triggers mtime`;
- the launcher log defaults to `snakemake.logs.txt` (change with `--log FILE`); `-q` reduces the launcher's own output.

---

## 6. Cluster submission

### 6.1 The four profiles

Profiles pin cluster parameters in the repository (`workflow/profile/<name>/config.yaml`, see `workflow/profile/README.md`):

| profile | scheduler | behavior |
|---|---|---|
| `default` | local / standalone server | no cluster command, executes locally |
| `pbs` | PBS (Torque), snakemake 7.x classic interface | `qsub -V -N {rule} -l select=1:ncpus={threads}:mem=<mem> -l walltime=<sec> -j oe`; walltime in seconds avoids `[HH:]MM:]SS` ambiguity |
| `sge` | SGE, snakemake 7.x classic interface | `qsub -V -N {rule} -l ncpus={threads} -l h_vmem=<mem> -l h_rt=<sec> -j oe` |
| `slurm` | SLURM, snakemake 7.x classic interface | `sbatch --parsable -J {rule} -c {threads} --mem=<mem> --time=<min> -o slurm-{rule}-%j.out` |

All four profiles pin `keep-going: true`, `rerun-incomplete: true`, `printshellcmds: true`.

### 6.2 Auto detection

`--profile auto` (the default; presettable via the `RUN_PROFILE` environment variable) probes in this order:

1. `sbatch` on PATH → **slurm**;
2. `qsub` present with `SGE_ROOT` set → **sge**; `qsub` present without `SGE_ROOT` → **pbs**;
3. neither → **default** (local).

An explicitly requested profile whose scheduler command is missing (e.g. pbs without qsub) errors out immediately. `--profile local` is equivalent to `default`.

### 6.3 The resource model

The resources requested from the cluster come from three layers; later layers win:

1. **Repository defaults**: `config/resources.yaml` (the table in §4.3);
2. **Project override**: a `resources.yaml` copied into the project directory (auto-detected by run.sh; or point `CHIP_RESOURCES_CONFIG` at any file);
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

The concurrency cap is `-j N` (default 3; the pbs/sge/slurm profiles pin `jobs: 20`, overridden by `-j`).

### 6.5 Environment variable presets

Common options can be preset via environment variables (explicit command-line arguments win): `RUN_PROFILE`, `CHIP_JOBS`, `CHIP_QUEUE`, `CHIP_PARTITION`, `CHIP_MEMORY`, `CHIP_RUNTIME_MIN`, `CHIP_SGE_MEMORY_RESOURCE`, `CHIP_SCHEDULER_EXTRA`, `CHIP_RETRIES`, `CHIP_LOG`, `CHIP_SOFTWARE_CONFIG`, `CHIP_RESOURCES_CONFIG`.

### 6.6 PBS output collection and the snakemake 8 note

- In PBS mode, qsub `.o*` job logs land at the project-directory root; **after a successful run the launcher collects them into `results/logs/`** (they stay in place on failure for debugging);
- with snakemake ≥ 8 the launcher prints a warning: 8.x moved to the executor plugin system (`--cluster` semantics changed), and this repository's four profiles target the 7.x classic interface — **test on a small scale before real 8.x cluster runs** (a `-n` dry-run first, then a single-sample trial, is recommended).

---

## 7. Output interpretation

### 7.1 Working-directory layout

```
workdir/
├── config.local.yaml        # optional: working-directory config overlay (auto-layered, not committed)
├── sample_info.csv          # sample table (filename free; pointed to by config.grouplist)
├── snakemake.logs.txt       # launcher/Snakemake main log
├── 1.rawdata/               # raw inputs (the only data directory outside results/):
│                            #   {sample}_1.fq.gz / {sample}_2.fq.gz
└── results/                 # ALL derived artifacts (rename via results_dir in config)
    ├── 0.index/             # bowtie2 index (bowtie2*.bt2; reusable across projects)
    ├── 2.cleandata/         # trimmed fastq: {sample}_1_val_1.fq.gz / {sample}_2_val_2.fq.gz
    │   └── fastqc/          # per-sample FastQC + multiqc/multiqc_report.html
    ├── 3.align/bowtie2/     # {sample}_sorted.bam(.bai), {sample}_rmdup.bam(.bai), {sample}_dup_metrics.txt
    ├── 4.peak/              # {group}_peaks.{narrowPeak,broadPeak}, {group}_summits.bed, {group}_FE.bw
    │   ├── replicates/      # peak.replicate stage: {group}/{sample}_peaks.{narrowPeak,broadPeak}
    │   ├── idr/             # peak.replicate stage: {group}/{a}__vs__{b}.narrowPeak pairwise IDR
    │   ├── blacklist_filtered/   # blacklist stage: filtered copies of the pooled/final peak sets
    │   ├── samples/         # bigwig.per_sample stage: {sample}.bw normalized coverage tracks
    │   └── anno_result/     # {group}.Anno.xls, Peakanno_PeakDistributions.pdf
    ├── 5.QC/
    │   ├── frip/            # {group}__{sample}.frip.tsv, FRiP_summary.tsv
    │   ├── replicate_peaks/ # peak.replicate stage: Replicate_summary.tsv (+ MultiQC table)
    │   ├── tss/             # qc.tss stage: {sample}_TSSE.txt, profile plots, TSSE_summary.tsv
    │   ├── organelle/       # qc.organelle stage: {sample}_idxstats.tsv, Organelle_summary.tsv
    │   ├── blacklist/       # blacklist stage: blacklist_summary.tsv (before/after counts)
    │   ├── spp/             # optional: {sample}_NSC.txt / _RSC.txt / _fragment_len.txt, NSC_RSC_mqc.tsv
    │   ├── deeptools/       # correlation heatmap / PCA / fingerprint / fragment size / gene-region signal profile
    │   ├── software_versions.yaml   # tool versions actually resolved for this run (incl. Snakemake)
    │   └── logs/            # QC rule logs (e.g. software_versions.log.txt)
    ├── 6.motif/             # motif stage: {group}/ HOMER results (known + de novo)
    ├── 6.diffbind/          # diffbind stage: {A}__vs__{B}/ sample sheet, DB tables, plots
    └── logs/                # per-rule logs (PBS .o job logs are collected here after success)
```

### 7.2 Results quick reference

| Result | Path |
|---|---|
| QC summary (fastqc + bowtie2 + picard + FRiP + NSC/RSC) | `results/2.cleandata/fastqc/multiqc/multiqc_report.html` |
| Trimmed fastq | `results/2.cleandata/{sample}_1_val_1.fq.gz`, `{sample}_2_val_2.fq.gz` |
| Aligned BAM / deduplicated BAM / dedup metrics | `results/3.align/bowtie2/{sample}_sorted.bam`, `{sample}_rmdup.bam`, `{sample}_dup_metrics.txt` |
| Peak files / summits | `results/4.peak/{group}_peaks.{narrowPeak,broadPeak}`, `results/4.peak/{group}_summits.bed` |
| Per-replicate peaks / pairwise IDR / final reproducible set (`peak.replicate.enabled`) | `results/4.peak/replicates/{group}/{sample}_peaks.*`, `results/4.peak/idr/{group}/{a}__vs__{b}.narrowPeak`, `results/4.peak/{group}_IDR_peaks.narrowPeak` or `{group}_consensus_peaks.broadPeak` (+ `_support.bed`) |
| Replicate summary table | `results/5.QC/replicate_peaks/Replicate_summary.tsv` |
| TSS enrichment (`qc.tss: true`) | `results/5.QC/tss/TSSE_summary.tsv` |
| Organelle fraction (`qc.organelle: true`) | `results/5.QC/organelle/Organelle_summary.tsv` |
| Blacklist-filtered peaks / counts (`blacklist` set) | `results/4.peak/blacklist_filtered/`, `results/5.QC/blacklist/blacklist_summary.tsv` |
| Per-sample normalized bigWigs (`bigwig.per_sample`) | `results/4.peak/samples/{sample}.bw` |
| HOMER motif results (`motif.enabled`) | `results/6.motif/{group}/` |
| Differential binding tables/plots (`diffbind.enabled`) | `results/6.diffbind/{A}__vs__{B}/DB_results.tsv` (+ `DB_significant.tsv`, plots) |
| Signal-track bigWig (fold enrichment) | `results/4.peak/{group}_FE.bw` |
| Peak annotation tables and plots | `results/4.peak/anno_result/{group}.Anno.xls`, `Peakanno_PeakDistributions.pdf` |
| bowtie2 index | `results/0.index/bowtie2*.bt2` |
| FRiP summary | `results/5.QC/frip/FRiP_summary.tsv` |
| NSC/RSC summary (`qc.nsc_rsc: true`) | `results/5.QC/spp/NSC_RSC_mqc.tsv` |
| deeptools QC (correlation/PCA/fingerprint/fragment size/gene-region signal) | `results/5.QC/deeptools/` |
| Software version record | `results/5.QC/software_versions.yaml` |
| Per-rule logs | `results/logs/` |

### 7.3 QC reference thresholds (ENCODE)

| Metric | Reference standard | Source |
|---|---|---|
| FRiP | TF ≥ 1% (ideally 5%+); relax as appropriate for histone marks | ENCODE |
| NSC | ≥ 1.05, ideally ≥ 1.1 (needs `qc.nsc_rsc: true`) | ENCODE |
| RSC | ≥ 0.8, ideally ≥ 1 (needs `qc.nsc_rsc: true`) | ENCODE |
| Alignment rate | typically ≥ 70% | empirical |

Per-sample FRiP and NSC/RSC values are injected into the MultiQC report (custom content) and can be viewed directly in `multiqc_report.html`; the per-sample NSC/RSC raw text lives in `results/5.QC/spp/`.

---

## 8. FAQ

**Q1: The sample table errors with "missing columns / illegal characters / inconsistent group"?**
Walk through §3.2 item by item: is the 6-column header complete; do `sample_id`/`group` contain commas, spaces, consecutive `__`, or start with `-`; are `seqtype`/`role`/`layout`/`peak_type` spelled within the whitelists; do all rows in one `group` agree on the three columns; does every `group` have a `treat`. Error messages carry line numbers — jump straight to the offending line.

**Q2: The launcher prints "[config warning] reference file does not exist (verify before running)"?**
This is a parse-time **warning, not an abort** (so dry-run/lint work on machines without the reference files). Before a real run, point `genome_fa`/`gtf`/`bed`/`chromsize` at real server paths (explicit keys override the species preset) and double-check `genome_size`; with references missing, rules such as bowtie2 indexing will fail at execution time.

**Q3: What happens if I run without a project config?**
The config chain falls back to the repository `config/config.yaml` (rice example), and the sample-table `grouplist` default is `config/samples.csv` (resolving against the working directory first, ultimately falling back to the repository template — the template holds an 8-row mixed-assay example that parses fine). The workflow then errors at the missing fastq (the `1.rawdata` data for the template samples does not exist). For real projects, follow §2 steps 3/4 and put the sample table and config in the working directory.

**Q4: A job was killed by the cluster for exceeding walltime/memory?**
Insufficient walltime or memory. Three fixes, in order of preference: (1) targeted per-rule override in a project `resources.yaml` (§4.3); (2) global `--runtime` / `--memory` overrides (§6.3); (3) `--retries N` so sporadic failures retry automatically. To locate the rule: find the failed rule name in the main log `snakemake.logs.txt`, then read `results/logs/<rule>.log` for details.

**Q5: Can snakemake 8.x use the cluster profiles?**
Parsing and dry-run work; but 8.x removed the classic `--cluster` semantics (replaced by executor plugins), and the four profiles' cluster commands target 7.x. The launcher warns when it detects 8.x. The safe route: use snakemake-minimal 7.32.4 from `workflow/environment.yaml` (the reference version) on cluster environments, or with 8.x verify with `--validate-only`/`-n` first and test on a small scale.

**Q6: How do I slow down under cluster quota/concurrency limits?**
Lower the concurrency cap with `-j N`; throttle submission and status checks with `--max-jobs-per-sec` / `--max-status-per-sec`; on PBS/SGE point `--queue` at the queue holding the quota; if still limited, `--scheduler-extra` can append site-required qsub/sbatch arguments (e.g. account/project numbers).

**Q7: How do I resume after an interruption? What about lock errors?**
Snakemake skips completed steps by output timestamps — **rerun the same command to resume** (profiles pin keep-going and rerun-incomplete). If you hit a working-directory lock error (`Directory cannot be locked`, usually from a force-killed job), run `bash run.sh -P . --unlock` and rerun.

**Q8: Why does CUT&Tag skip deduplication by default?**
The CUT&Tag literature recommends keeping PCR duplicates (Kaya-Okur et al. 2019): libraries undergo few amplification cycles and duplicates carry real signal. Flip it with `dedup.cuttag: true` (per-assay switches; chip/atac/faire deduplicate by default).

**Q9: Is a control required?**
No: for a group with no control rows, MACS2 automatically omits `-c` and falls back to local lambda estimation, at the cost of peak quality. Controls are recommended for strong-background histone marks (e.g. H3K27me3 CUT&Tag) and for FRiP computation.

**Q10: Can I run single-end (SE) data?**
This version supports PE only (the sample table's `layout` is fixed to `PE`; anything else fails validation).

**Q11: How do I switch species?**
Either set `species: "hsa"` (or another preset in `config/species.yaml`), or override the reference file set `genome_fa`/`gtf`/`bed`/`chromsize` + `genome_size` explicitly (§4.2); everything else stays the same. The GTF must contain gene_id/transcript_id for ChIPseeker to build its TxDb.

**Q12: How do I validate a new deployment or a workflow change?**
```bash
make check                       # unit tests + bash -n syntax checks (no snakemake needed)
make lint                        # static suite (bash/shellcheck/py/R/yaml/snakemake --lint; missing tools skipped)
bash tests/run_test.sh           # synthetic-data dry-run regression (needs snakemake + python3/PyYAML)
bash tests/run_test.sh --replicate   # replicate/IDR scenario dry-run (adds a 2-treat broad group)
bash tests/run_test.sh --qc-full     # extended QC scenario dry-run (tss + organelle + blacklist)
bash tests/run_test.sh --motif       # motif stage scenario dry-run (dummy genome tag)
bash tests/run_test.sh --diffbind    # differential binding scenario dry-run (adds a contrast group)
bash tests/run_test.sh --real-run    # end-to-end run + output assertions (server validation; needs the full analysis environment)
```
Test data is generated by `tests/make_testdata.py` with a fixed seed (2 × 100kb chromosomes, 3 chip + 2 atac samples); the dry-run defaults to 50000 read pairs per sample (CI passes `--reads 2000`). Before changing workflow code, read the documentation-sync checklist in [CONTRIBUTING](../CONTRIBUTING.md).

**Q13: How do I add IDR / replicate-aware peak analysis?**
Set `peak.replicate.enabled: true` (§5.3). Everything else is automatic per the sample table's replicate structure (a group's treat rows are its replicates). You need the external `idr` tool installed (python2; e.g. `conda create -n idr -c bioconda idr=2.0.4`, see the environment note in §5.3). Single-replicate groups keep using their pooled peaks.

**Q14: Where does the blacklist come from?**
ENCODE maintains blacklists for human/mouse; for rice there is none, so the workflow leaves `blacklist` empty by default. If you have one (self-built from repeated-artifact evidence, or from a closely related assembly), point the key at the BED file and FRiP/annotation switch to the filtered peak copies (§5.4).

**Q15: My ATAC library shows a huge organelle fraction — what do I do?**
Enable `qc.organelle` and check `results/5.QC/organelle/Organelle_summary.tsv`. High chloroplast fractions (common in plant ATAC from green tissues) waste sequencing depth; aligning against a nuclear-only reference, or in silico removing organelle-mapped reads upstream, are the standard remedies. The metric is informational — the workflow never filters BAMs on it.

**Q16: How do I run a differential binding comparison?**
Give the two conditions their own sample-table groups (≥ 2 treat replicates each), optionally label them with the `condition` column (plus `batch` for sequencing-batch blocking), then set `diffbind.enabled: true` and `diffbind.contrasts: [["groupA", "groupB"]]` (§5.6). Best combined with `peak.replicate.enabled: true` so DiffBind counts over per-replicate peak sets.

**Q17: Do I need HOMER for the motif stage?**
Yes — HOMER (findMotifsGenome.pl + a configured genome) is an external distribution, deliberately not in the conda template. Point the software.yaml `paths: homer_findmotifs` entry at the binary, install a genome for your assembly via `configureHomer`, and set `motif.homer_genome` (§5.5).
