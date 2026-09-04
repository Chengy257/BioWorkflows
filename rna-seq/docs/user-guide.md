# RNA-seq Workflow User Guide

> Updated: 2026-09-04 (unified software/R runtime; flattened output layout; dedicated resources.yaml; config layering)
> Audience: analysts running this pipeline on their cluster/server
> Since v0.8.0 the workflow no longer creates multiple Conda environments per rule; it uniformly reuses the user's existing software runtime. `config/software.yaml` manages the main Conda environment, Rscript/R library, external tools, and databases. The v0.7.0 per-rule `threads / mem_mb / runtime` scheduling is retained, now in the dedicated `config/resources.yaml`.

---

## 1. Feature overview

This pipeline targets bulk RNA-seq (reference-genome based) data and offers four run modes via the `pipeline` parameter:

| pipeline | Content | Final outputs |
|---|---|---|
| `upstream` | QC -> alignment -> gene quantification | cleaned data, BAM, count matrix, TPM table |
| `deg` | upstream + differential expression + enrichment + between-group comparison | DEG results, GO/KEGG/GSEA, between-group comparison |
| `as` | upstream + StringTie transcript assembly/quantification (STAR uses assembly-optimized parameters) | per-sample assembly GTFs, merged GTF, isoform quantification |
| `lncrna` | upstream + de novo lncRNA identification | candidate lncRNA set and its expression matrix |

Main analysis line: `FastQC -> Trim Galore (incl. FastQC) -> MultiQC whole-pipeline summary -> (strandedness inference) -> STAR two-pass alignment -> featureCounts quantification -> DESeq2 differential analysis -> clusterProfiler GO/KEGG enrichment + GSEA`.

---

## 2. System requirements

| Dependency | Notes |
|---|---|
| Linux (SGE/SLURM cluster or single machine) | `qsub` detected -> SGE, `sbatch` -> SLURM, otherwise local |
| Snakemake >= 7 | Must be findable from the main runtime or the current PATH; development verified with 7.32.4 |
| Python 3 + PyYAML | Used by `run.sh` to parse `software.yaml` and run preflight checks |
| Main software environment | Reusing an existing `rna-seq` Conda environment is recommended; `environment.type: system` (use the PATH directly) is also supported |
| R | Can come from the main Conda env, or an independent Rscript via `r.rscript`; R libraries controlled by `r.lib_paths` |
| lncRNA external dependencies | CPC2, CNCI, pfam_scan.pl, Pfam DB, NR DIAMOND DB; CNCI can specify its own Python2 |

Without a ready environment, `workflow/environment.yaml` provides an **optional all-in-one environment template**:

```bash
mamba env create -f workflow/environment.yaml
Rscript -e 'install.packages("aPEAR", repos="https://cloud.r-project.org")'
```

This command is executed explicitly by the user; Snakemake itself never creates or modifies software environments.

---

## 3. Project structure

```
rna-seq/
├── run.sh                    # unified launcher (repo subdirectory root)
├── Makefile                  # make check / lint / test development checks
├── workflow/                 # pipeline code
│   ├── Snakefile             # single entry point
│   ├── rules/                # common / meta / align / quant / deg / as / lncrna
│   ├── scripts/              # R/Python/Shell scripts (stable names)
│   ├── environment.yaml      # optional all-in-one environment template (not auto-created)
│   ├── profile/              # scheduler profiles (sge / slurm / default)
│   └── multiqc_config.yaml   # custom MultiQC
├── config/                   # configuration templates (copy into the project directory before running)
│   ├── config.yaml           # analysis parameters and reference resources
│   ├── config.template.yaml  # minimal project config template (or copy as config.local.yaml)
│   ├── software.yaml         # software/R runtime/database configuration
│   ├── resources.yaml        # per-rule scheduler resources (threads/mem_mb/runtime_min)
│   ├── species.yaml          # osa/hsa species resource presets
│   └── samples.csv           # sample table template
├── example/                  # example project (includes a real sample table reference)
├── tests/                    # regression tests (make_testdata.py / run_test.sh / lint.sh / check_outputs.py)
├── docs/
├── CONTRIBUTING.md           # contribution guide (doc-sync checklist / CHANGELOG requirements)
└── results/                  # run outputs (gitignored; root configurable via results_dir)
```

---

## 4. Preparation

### 4.1 Create a new analysis project

```bash
mkdir -p myproject && cd myproject
mkdir -p 1.rawdata
# Place raw data: paired-end {sample}_1.fastq.gz + {sample}_2.fastq.gz (or .fq.gz)
#                 single-end {sample}.fastq.gz (or .fq.gz)
cp /path/to/repo/config/config.yaml .
cp /path/to/repo/config/software.yaml .
cp /path/to/repo/example/samples_231107XTL.csv samples.csv   # or build your own sample table
```

⚠️ **Sample ID naming rules**: must not contain `-`, spaces, or `/`; avoid names that are prefixes of each other (the validator warns about prefix conflicts).

### 4.2 Sample group table

Columns: `id,group`, optionally `layout` (`PE`/`SE`/`auto`; the current version detects from the actual files, so `auto` is fine):

```csv
id,group,layout
WT_1,control,auto
drugA_1,drugA,auto
```

The control group name is set by `control_group` in the config (default `control`). The sample table is validated automatically at startup; it can also be run manually:
`python3 workflow/scripts/validate_samples.py samples.csv control`

### 4.3 Configuration files

Copy `config/config.yaml` into the project directory and **replace the `/path/to/...` placeholders with real paths**:

| Key | Description |
|---|---|
| `pipeline` | Run mode (the run.sh argument overrides it) |
| `results_dir` | Run outputs root directory (default `results`, relative to the run directory) |
| `SampleListFile` | Sample table (relative to the run directory) |
| `control_group` / `threads` | Control group / legacy default thread cap (explicit `resources.<rule>.threads` can override) |
| `FoldChange` / `padj` / `pca_ntop` | Differential thresholds / number of genes used for PCA |
| `batch_correction` | DESeq2 batch correction (default `F`; with `T` the design becomes `~ batch + group`, and **the sample table needs a third `batch` column**) |
| `species` | `osa` / `hsa`; reference resources not set explicitly fall back to the species presets in `config/species.yaml` |
| `bed` / `genome` / `gtf` / `annotation_tsv` | Reference resources (explicit top-level settings take priority over species.yaml) |
| `star_extra_args` | Extra STAR arguments |
| `lncrna.*` | lncRNA analysis parameters, only `gtf`, `gtf_PcGs`, `threads` remain; software and databases moved to `software.yaml` |

At parse time the workflow hard-validates the merged configuration (`validate_config` in `workflow/rules/common.smk`): required keys, numeric checks (`FoldChange`/`padj`), `batch_correction` must be `T`/`F`, and unknown `resources` fields are rejected — all issues are reported in one aggregated error. Reference paths still containing `/path/to/` placeholders produce warnings.

### 4.3.1 Software and R runtime: `software.yaml`

A project-level `software.yaml` takes priority over the repository default `config/software.yaml`; it can also be given explicitly with `run.sh --software FILE`. Recommended configuration example:

```yaml
environment:
  type: conda
  conda_prefix: "/share/miniconda/envs/rna-seq"   # HPC deployments should prefer prefix
  # conda_name: "rna-seq"                         # or resolve by environment name
  strict: true

r:
  rscript: "/opt/R/4.3.3/bin/Rscript"             # or just Rscript, using the main environment PATH
  version: "4.3"
  version_check: major_minor
  lib_paths:
    - "/share/Rlibs/4.3"
    - "/home/user/R/4.3/library"
  lib_mode: prepend
  package_sources:
    org.Osativa.eg.db: "/share/db/org.Osativa.eg.db.tar.gz"  # for installation hints only

tools:
  # Ordinary software does not need to be listed one by one; all of the following are optional overrides
  # star: "/opt/STAR/2.7.11b/bin/STAR"
  cpc2: "/share/software/CPC2/bin/CPC2.py"
  pfam_scan: "/share/software/PfamScan/pfam_scan.pl"
  cnci_python: "/share/software/python2/bin/python2"

paths:
  cnci_dir: "/share/software/CNCI"

databases:
  pfam: "/share/database/Pfam"
  nr_diamond: "/share/database/NR/nr.dmnd"
```

Resolution priority is: **explicit tool path > main Conda `bin` > original system PATH**. `r.rscript` is resolved independently, so R can differ from the main Conda environment. `r.lib_paths` is propagated uniformly via `R_LIBS_USER` to Snakemake rules, `enrich.sh`, `DEGgroupCompare.sh`, and their child Rscript processes; `prepend` does not override R's own base/site libraries.

Recommended before a production run:

```bash
bash run.sh -p deg -P /path/to/project --software /path/to/project/software.yaml --check-software
bash run.sh -p deg -P /path/to/project --software /path/to/project/software.yaml --check-r
```

`--check-software` checks executables, the R version/R packages, and the lncRNA databases per pipeline; `--check-r` checks R only. The pipeline **never** calls `install.packages()` during a run. When a custom OrgDb is missing, install it into a directory listed in `r.lib_paths` first.

### 4.3.2 Configuration layering

Configurations merge in the following order (later layers win):

1. Repository defaults: `config/config.yaml` (species presets from `config/species.yaml` are also loaded);
2. Project configuration: given with `-c FILE` or the positional argument;
3. Extra layer: `config.local.yaml` in the project directory is picked up automatically, or any file via `-l/--extra-config FILE`. This layer is intended for machine-local overrides that stay out of version control.

For a minimal starting point, copy `config/config.template.yaml` into the project and set only the keys you need — everything else inherits the repository defaults.

### 4.3.3 Per-rule scheduler resources: `resources.yaml`

Scheduler requests (threads / mem_mb / runtime_min per rule) live in the dedicated `config/resources.yaml` — the single place to tune cluster requests. The `resources:` block was removed from `config.yaml`. To override, copy the file into the project directory as `resources.yaml`; `run.sh` auto-detects it (explicit override: `RNASEQ_RESOURCES_CONFIG`). The legacy top-level `threads` value in `config.yaml` still acts as a global cap on every rule's threads.

| resource key | threads | mem_mb | runtime_min |
|---|---:|---:|---:|
| `software_versions` | 1 | 1024 | 10 |
| `trim` | 4 | 8000 | 180 |
| `multiqc` | 2 | 4000 | 60 |
| `star_index` | 12 | 48000 | 480 |
| `star_align` | 12 | 32000 | 360 |
| `mapping_stat` | 1 | 2000 | 30 |
| `strandedness` | 1 | 4000 | 60 |
| `featurecounts` | 8 | 16000 | 240 |
| `count_merge` | 1 | 4000 | 60 |
| `deseq2` | 4 | 16000 | 240 |
| `enrichment` | 2 | 12000 | 360 |
| `deg_compare` | 4 | 12000 | 360 |
| `stringtie` | 8 | 16000 | 360 |
| `gtf_merge` | 4 | 12000 | 240 |
| `isoform_expr` | 8 | 16000 | 360 |
| `lnc_coding` | 8 | 16000 | 720 |
| `pfam` | 16 | 24000 | 720 |
| `nr` | 16 | 32000 | 1440 |
| `lnc_final` | 2 | 8000 | 120 |
| `lnc_featurecounts` | 8 | 16000 | 240 |
| `lnc_count_merge` | 1 | 4000 | 60 |

Thread compatibility rule: when `resources.<rule>.threads` is not set, a normal rule's built-in threads are capped by the top-level `threads`; `pfam` / `nr` keep falling back to the legacy `lncrna.threads` without an explicit override. Explicit `resources.<rule>.threads` always wins.

### 4.4 The `run.sh` launcher

The old positional-argument invocation stays compatible:

```bash
bash run.sh <upstream|deg|as|lncrna> <project_dir> [config.yaml] [jobs]
# e.g.: bash run.sh deg /path/to/myproject /path/to/myproject/config.yaml 10
```

Explicit arguments are recommended; run `bash run.sh --help` for full help:

```bash
# local / automatic scheduling
bash run.sh -p deg -P /path/to/myproject -c /path/to/myproject/config.yaml -j 20

# check the unified software/R runtime
bash run.sh -p deg -P /path/to/myproject --software /path/to/myproject/software.yaml --check-software

# validate samples + software/R, without starting jobs
bash run.sh -p deg -P /path/to/myproject -c /path/to/myproject/config.yaml --validate-only

# dry-run
bash run.sh -p deg -P /path/to/myproject -c /path/to/myproject/config.yaml --dry-run

# SLURM
bash run.sh -p deg -P /path/to/myproject --profile slurm --partition compute --memory 32G -j 40

# SGE
bash run.sh -p deg -P /path/to/myproject --profile sge --queue all.q --runtime 720 --retries 2
```

Environment arguments: `--software FILE`, `--check-software`, `--check-r`, `--skip-software-check`. Scheduling arguments: `--memory VALUE` (globally overrides per-rule memory), `--runtime MIN` (globally overrides per-rule walltime), `--latency-wait SEC`, `--max-jobs-per-sec N`, `--max-status-per-sec N`, `--scheduler-extra "..."`, `--unlock`, `--log FILE`; arguments after `--` are passed to Snakemake verbatim. Most long options also support the `--option=value` form.

Script behavior:
1. Parses `software.yaml`, injecting the main Conda prefix, tool overrides, `RNASEQ_RSCRIPT`, and `R_LIBS_USER` into the current process, inherited by SGE/SLURM jobs;
2. The pipeline and analysis configuration are injected into the Snakefile via `RNASEQ_PIPELINE` / `RNASEQ_CONFIG`; the resources file is injected via `RNASEQ_RESOURCES_CONFIG` (project-local `resources.yaml` wins, otherwise the repository default);
3. Before a production run it performs the software/R preflight per pipeline and validates the sample table automatically; with `batch_correction: "T"` the `batch` column is mandatory and every row must be non-empty;
4. Automatic scheduling order is SGE (`qsub`) -> SLURM (`sbatch`) -> local; can be forced via `--profile` / `RUN_PROFILE`;
5. By default SGE/SLURM request per-rule memory/walltime dynamically from `{resources.mem_mb}` and `{resources.runtime_*}`; `--memory` / `--runtime` force global overrides; the SGE memory resource name can be overridden with `--sge-mem-resource`;
6. Logs go to `$project_dir/snakemake.logs.txt` by default; on failure the exit code, duration, and log location are reported.

Calling Snakemake directly skips the automatic `software.yaml` resolution, so it is only recommended when the required environment is already activated and PATH/R_LIBS are set up manually:

```bash
RNASEQ_PIPELINE=deg RNASEQ_CONFIG=$PWD/config.yaml \
  snakemake -s /path/to/repo/workflow/Snakefile --profile /path/to/repo/workflow/profile/sge -j 10
```

---

## 5. Running and monitoring

```bash
# 1. Dry-run preview (nothing is executed)
cd myproject
RNASEQ_PIPELINE=deg RNASEQ_CONFIG=$PWD/config.yaml \
  snakemake -s ../workflow/Snakefile --profile ../workflow/profile/default -n --quiet

# 2. Production run
bash ../run.sh deg . config.yaml 10

# 3. Progress and logs
tail -f snakemake.logs.txt
ls results/logs/                  # per-step logs
qstat                             # SGE job queue
```

Resume behavior: Snakemake automatically skips completed steps based on output files. Since v0.8 no Conda environments are created and no R packages are installed during a run; the runtime should be ready before running.

---

## 6. Results layout (under results/)

| Directory/file | Content |
|---|---|
| `multiqc/multiqc_report.html` | **Whole-pipeline** QC summary (FastQC/Trim Galore/STAR/featureCounts) |
| `software_versions.yaml` | The runtime type/Conda prefix, Rscript/R library, tool paths, and database paths actually resolved for this run |
| `2.cleandata/trim/` | Trimmed fastq; `{sample}(_1/_2)?_trimming_report.txt` canonicalized trimming reports |
| `2.cleandata/trim/fastqc/` | Post-trim FastQC reports (html/zip) |
| `3.align/{sample}_Aligned.sortedByCoord.out.bam(.bai)` | Alignment results (STAR direct sorted BAM) |
| `3.align/{sample}.strandedness` (and `_infer_experiment.out`) | Strandedness inference results |
| `3.align/mapping_stat.xls` | Alignment statistics summary (parsed per sample, complete IDs; Total_Reads is the PE read-pair count) |
| `4.expression/{sample}.count` / `.log` | Per-sample featureCounts quantification and assignment statistics |
| `4.expression/count.matrix.tsv` / `GeneExpression_TPM.xls` / `_FPKM.xls` / `GeneCount_Assigned_logs.xls` | Expression matrices |
| `5.DEG/Diff_Expr_Analysis_Results/` | Differential analysis: clustering heatmap, PCA (ntop configurable), per `{treat}_vs_{control}_DESeq2.output.tsv`, volcano/MA plots, sessionInfo.txt |
| `5.DEG/DEGs/` / `GO_KEGG_enrich/` / `GSEA_enrich_GO/` | DEG lists and annotation / enrichment results / GSEA |
| `6.DEGcompare/` | Between-group DEG set intersections (venn/upset) and cross-enrichment comparison |
| `4.assembly/stringtie/` + `4.assembly/isoform/` (as pipeline) | StringTie assembly, `merged.gtf`, gffcompare reports, isoform FPKM/TPM |
| `4.lncrna/` (lncrna pipeline) | `assembly/` (stringtie assembly + classcode_u candidates), `coding_predict/` (CPC2/CNCI/length filtering), `pfam_search/` + `nr_search/` (domain/protein hits), `final/` (`final_lncRNA.fa/gtf`), `expression/` (lncRNA quantification matrices) |
| `logs/` | Per-rule logs |

The differential result `*_DESeq2.output.tsv` columns: `gene_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, FoldChange, per-sample normalized counts, type(Up/Down/Unsig)`; thresholds are controlled by `padj`/`FoldChange` in the config.

---

## 7. Known limitations and caveats (v0.8.0)

| # | Item | Impact | Recommendation |
|---|---|---|---|
| 1 | KEGG enrichment depends on KEGG REST API internet access | On nodes without internet the step is skipped automatically (GO results are unaffected) | Re-run KEGG on an internet-connected node; an offline gmt-based alternative is possible via clusterProfiler's local gmt workflow |
| 2 | Differential comparison naming `{treat}_vs_{control}` | Underscores in the control group name may break downstream prefix parsing | Avoid underscores in the control group name |
| 3 | Batch correction is off by default (`batch_correction: "F"`) | Must be enabled explicitly when needed | Set `batch_correction: "T"` in the config and add a third `batch` column to the sample table (see §8 Q2) |
| 4 | `mapping_stat.xls` Total_Reads takes the trim report R1 count (i.e. PE read pairs) | Different notion of "reads" than for single-end samples | For detailed counts rely on the original `*_trimming_report.txt` / `*_Log.final.out` |
| 5 | CNCI still requires Python2, and CPC2/pfam_scan/Pfam/NR are often external server installs | lncRNA steps fail when the environment is missing | Configure `cnci_python`, `cnci_dir`, CPC2/PfamScan, and database paths in `software.yaml`, then run `--check-software` first |
| 6 | End-to-end regression requires an already prepared unified RNA-seq software environment | This repo does not install test dependencies automatically | `bash tests/run_test.sh --pipeline deg` does it in one shot (data generated deterministically; see §8 FAQ) |

The historical review/roadmap documents were archived out of the repository (not in version control); retrieve them from a local archive or git history if needed.

---

## 8. FAQ

**Q1: Can single-end sequencing run?**
Yes. Put `{sample}.fastq.gz` (or `.fq.gz`) in rawdata; the pipeline auto-detects SE/PE from exact file names, and both naming schemes work in all pipelines.

**Q2: How are batch effects handled?**
Set `batch_correction: "T"` in the config and add a third `batch` column to the sample table; DESeq2 then runs with the `~ batch + group` design (default `"F"` disables it).

**Q3: How do I re-run only the differential analysis?**
When the count matrix already exists, run directly with `pipeline=deg`; the upstream steps are skipped automatically.

**Q4: Can I switch species?**
Yes (osa/hsa). Swap the reference resources at the alignment/quantification level; set `species: hsa` at the enrichment level to switch automatically to org.Hs.eg.db + ENSEMBL keyType + KEGG hsa (gene ids are stripped of version suffixes automatically). Other species need a self-built OrgDb and an extension following `workflow/scripts/run_enrichment.R`.

**Q5: Why do my differential results only contain vs control?**
That is the current design: DESeq2 only outputs "each treatment group vs control group" (`control_group` config). Arbitrary pairwise comparisons are handled by `6.DEGcompare`, which intersects/enrichment-compares DEG sets.

**Q6: How do I reuse an existing server Conda environment?**
Set `environment.type: conda` plus `conda_prefix` (recommended) or `conda_name` in the project `software.yaml`. No `conda activate` needed; `run.sh` adds that prefix's `bin` to PATH. If R is installed separately, also configure `r.rscript` and `r.lib_paths`. Use `--check-software` to verify the resolution result.

**Q7: Do the pre-v0.2.0 entry points still work?**
No. The old entry points `rna-seq-workflow/RNA-seq_*.smk` were deprecated and archived out of the repository; use `workflow/Snakefile` (launched via `run.sh`). Retrieve old code from a local archive or git history if needed.

**Q8: How do I verify a new environment deployment or pipeline change?**
```bash
make check                                # runtime resolver + R runtime propagation tests + shell syntax
make lint                                 # static checks (missing tools skip their section automatically)
make test                                 # everything CI runs (= check + lint)
bash tests/run_test.sh --pipeline deg    # end-to-end regression: synthetic data -> dry-run -> run -> assertions
```
Test data is generated by `tests/make_testdata.py` with a fixed seed (2 x 100 kb chromosomes, 60 simulated genes, 2 groups x 2 samples) plus an expected-differential-gene truth table; assertions include mapping_stat consistency and numeric checks, count matrix shape, and DEG direction checks. The `lncrna` pipeline depends on external tools and is not covered by automated tests. Before changing pipeline code, read the documentation sync checklist in [CONTRIBUTING](../CONTRIBUTING.md).
