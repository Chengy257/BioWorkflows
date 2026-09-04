# Changelog

All notable changes to this project are documented here. Format based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

Architecture consolidation release: flattened output layout, dedicated scheduler-resource configuration, an extra config layering mechanism, parse-time validation, and a shared cross-project layer. All text unified in English; the project is now licensed under Apache-2.0.

### Added
- Flattened the lncRNA/assembly output layout under the `results_dir` root: `4.assembly/4.1.Assembly_stringtie/` -> `4.assembly/stringtie/`, `4.assembly/4.2.IsoformExpr/` -> `4.assembly/isoform/`, `4.LncRNA/4.1.Assembly_stringtie/` -> `4.lncrna/assembly/`, `4.LncRNA/4.2.Coding_predict/` -> `4.lncrna/coding_predict/`, `4.LncRNA/4.3.Pfam_search/` -> `4.lncrna/pfam_search/`, `4.LncRNA/4.4.Nr_search/` -> `4.lncrna/nr_search/`, `4.LncRNA/4.5.Final_lncRNA/` -> `4.lncrna/final/`, and `5.expression/lncRNA/` -> `4.lncrna/expression/`.
- New dedicated `config/resources.yaml` for per-rule scheduler resources (threads / mem_mb / runtime_min). A project-local `resources.yaml` is auto-detected by `run.sh` (also settable via `RNASEQ_RESOURCES_CONFIG`); the legacy top-level `threads` still caps all rules globally.
- New config layering: repository defaults -> project config (`-c` or positional) -> `config.local.yaml` in the project directory (auto-detected) or `-l/--extra-config FILE`. Added `config/config.template.yaml` as the minimal project-config template.
- Parse-time hard validation in `workflow/rules/common.smk` (`validate_config`): required keys, numeric checks (`FoldChange`/`padj`), `batch_correction` T/F, a resources field whitelist, and an aggregated error report; placeholder `/path/to/` reference paths produce warnings.
- New `Makefile` with `make check` / `make lint` / `make test` development checks.

### Changed
- The `resources:` block was removed from `config/config.yaml` (which now points to `resources.yaml`); scheduler resources are tuned in one dedicated file.
- `run.sh` now sources `../shared/lib/launcher.sh`; `workflow/scripts/runtime_config.py` and `collect_versions.py` are thin wrappers over `shared/python/bioworkflows_runtime.py` / `bioworkflows_versions.py` (cross-project shared layer).
- CI moved to the repository root `.github/workflows/ci.yml` (rna-seq job: lint + runtime resolver tests + deg end-to-end regression); the sub-project-level `.github/` was removed.
- Documentation updated to the new architecture (output tree, resources section, config layering); the Chinese user guide was renamed to `docs/user-guide.md` and translated. `rna-seq-workflow/` and the archived review/roadmap documents remain out of the repository.
- All remaining Chinese comments, messages, and documentation unified in English.
- License changed to Apache-2.0.

### Fixed
- `workflow/environment.yaml` solves again (validated with a dry-run Conda solve): `python` is pinned to 3.10 because `rseqc=5.0.1` ships no newer Python builds, `r-base` is set to the 4.3 series, and `aPEAR` moved from a post-install CRAN step into the environment as `bioconductor-aPEAR` (its CRAN release requires a newer R than the solved series).

## [0.8.0] - 2026-09-03

Unified software and R runtime: the workflow reuses the user's existing server environment by default and no longer creates a dedicated Conda environment per rule.

### Added
- New `config/software.yaml`: centralizes the main runtime (`system` / Conda prefix / Conda name), Rscript, R version checks, R library paths, tool overrides, CNCI Python2, and the Pfam/NR databases.
- New `workflow/environment.yaml`, an optional all-in-one Conda environment template for new servers; created explicitly by the user only, never deployed automatically by Snakemake.
- New `workflow/scripts/runtime_config.py`: resolves the software configuration uniformly, builds PATH / `R_LIBS_USER`, exports `RNASEQ_TOOL_*` runtime variables, and runs executable / R / R package / database preflight per pipeline.
- `run.sh` gained `--software`, `--check-software`, `--check-r`, `--skip-software-check`; a `software.yaml` in the project directory is picked up automatically.
- New `tests/test_runtime_config.sh` and `tests/test_r_runtime_propagation.sh`: verify Conda prefix / Rscript / R library resolution, and that nested R sub-jobs of `enrich.sh` / `DEGgroupCompare.sh` inherit the same Rscript/R_LIBS.

### Changed
- All rules reuse the unified runtime; STAR, samtools, Trim Galore, MultiQC, StringTie, DIAMOND, bedtools, Rscript/Python etc. all support explicit executable overrides and are resolved from the main environment/PATH otherwise.
- The R runtime became a first-class standalone configuration: Rscript may differ from the main Conda environment, and `r.lib_paths` is propagated via `R_LIBS_USER` to direct rules and to the nested R calls of `enrich.sh` / `DEGgroupCompare.sh`.
- lncRNA's CPC2, CNCI, pfam_scan, Pfam DB, and NR DIAMOND DB moved from the analysis `config.yaml` to `software.yaml`; CNCI's Python2 interpreter can be specified independently.
- `software_versions.yaml` now records the actual runtime — Conda prefix, Rscript/R version/R library, resolved tool paths, and large database locations — instead of a list of environment recipes.
- CI and regression tests prepare one unified environment before running the workflow; Snakemake lint still checks other issues but deliberately ignores the "each rule should specify a conda/container" hint that conflicts with the current architecture.

### Removed
- Deleted the `workflow/envs/*.yaml` multi-environment setup from the production workflow, all rule `conda:` directives, and `use-conda: true` from the profiles.
- Removed the side effect of R scripts auto-running `install.packages()` at runtime; required R packages must be installed in the configured R library before running.
- `orgdb_tarball` is no longer passed as an analysis/run parameter; a local OrgDb tarball, if needed, can go into `r.package_sources` for preflight installation hints only.

### Migration
- Existing server environment: copy `config/software.yaml` into the project and set `environment.conda_prefix` (recommended) or `conda_name`; if the current PATH is already fully configured, `environment.type: system` works directly.
- No ready environment: explicitly run `mamba env create -f workflow/environment.yaml`, then point `software.yaml` at that environment.
- Before production runs, `run.sh ... --check-software` is recommended; the R configuration can be verified separately with `--check-r`.

## [0.7.0] - 2026-09-03

Introduced the per-rule resource model and completed Snakemake engineering normalization, letting SGE/SLURM request CPU, memory, and walltime dynamically per job, while making repository-level lint actually pass.

### Added
- All 24 executable rules declare `threads` and `resources.mem_mb/runtime_min/runtime_sec`; the defaults cover workloads from trim, STAR, featureCounts, DESeq2, and StringTie to Pfam and DIAMOND.
- `config/config.yaml` gained a `resources:` section to override `mem_mb` and `runtime_min` per rule, with optional `threads`; projects without the section automatically use the built-in defaults.
- `run.sh` gained `--runtime MIN` / `RNASEQ_RUNTIME_MIN`; the SGE/SLURM default submit templates read `{resources.mem_mb}` and the rule walltime respectively, with `--memory` / `--runtime` as global forced overrides.
- New `workflow/envs/utils.yaml` and `workflow/rules/meta.smk` for Python aggregation and software version recording.

### Changed
- Actual tool threads inside rules use Snakemake `{threads}` uniformly instead of reading the single `config[threads]`; the legacy top-level `threads` and `lncrna.threads` fallbacks are retained.
- Helper functions were consolidated into `common.smk` and metadata rules into `meta.smk`; scripts/paths/configs needed by shells are derived via `params` or input/output, reducing hidden global dependencies.
- Conda environments grew from 7 to 8; pure Python aggregation and metadata use the lightweight `utils` environment.
- `collect_versions.py` can accept the Snakemake version currently in use, avoiding the metadata environment changing the recorded value.

### Fixed
- `snakemake --lint` passes for all four pipelines; `tests/lint.sh` also passes bash, ShellCheck, Python, and Snakemake lint (Rscript was not installed on the machine at the time, so the R parse step was skipped by design).
- Fixed the unhandled-`cd` ShellCheck SC2164 issues in `enrich.sh`, `lncRNA_functions.sh`, and `tests/lint.sh`.

## [0.6.0] - 2026-09-03

Enhanced the run entry point and cluster scheduling, and fixed rule-parsing compatibility issues with Snakemake 7.32.4.

### Added
- `run.sh` gained a full English CLI: `--help` / `--version` / `--dry-run` / `--validate-only` / `--skip-validation` / `--unlock` / `--log`, while keeping the original four-positional-argument invocation.
- Scheduler auto-selection upgraded to SGE (`qsub`) -> SLURM (`sbatch`) -> local; added `--profile`, `--queue`, `--partition`, `--memory`, `--sge-mem-resource`, `--scheduler-extra`.
- Added failure retries and scheduler throttling: `--retries`, `--latency-wait`, `--max-jobs-per-sec`, `--max-status-per-sec`; arguments after `--` are forwarded to Snakemake; common long options support `--option=value`.
- With `batch_correction: "T"`, startup automatically requires the sample table to contain a non-empty `batch` column.

### Changed
- `validate_samples.py` moved to an English `argparse` CLI with a structured summary, `--require-batch`, and `--strict-warnings`, staying compatible with the old `samples.csv control` invocation.
- Runtime messages and core comments in `workflow/Snakefile`, `workflow/rules/common.smk`, and the three profiles were unified in English.
- README and the user guide were synchronized with the new scheduling logic and launcher arguments.

### Fixed
- Fixed Snakemake 7.32.4 treating the named input/output field `count` as a reserved name, which broke workflow parsing; the relevant fields were renamed and all four pipelines parse normally.
- Fixed PE sample R1 files (e.g. `ctrl_1_1.fastq.gz`) being misread by the SE trim rule as independent samples and triggering `AmbiguousRuleException`; all `{sample}` wildcards are now strictly restricted to real IDs from the sample table.

## [0.4.1] - 2026-09-03

Review retrospective fixes: closed the last configuration gaps plus several engineering details. All three behavior fixes are backward compatible (without the new config keys, behavior is identical to v0.4.0).

### Fixed
- **DEGgroupCompare species argument breakage**: `rules/deg.smk` previously did not pass `species`/`orgdb_tarball` into `DEGgroupCompare.sh`, so `run_deg_compare.R` always ran the default osa branch — with `species=hsa`, the between-group GO enrichment silently failed or used the wrong OrgDb, while osa depended on the side effect of upstream enrichment steps installing the OrgDb first. Both arguments are now propagated through rule -> script -> job list (passed explicitly per job).
- `DEGgroupCompare.sh`: the temp file for intersection annotation moved from the run-directory `tmp.gene` to `mktemp` + `trap` cleanup (no longer left behind when fgrep has no hits); partial failures of parallel jobs now emit an explicit WARN (no longer fully silent; fault-tolerance semantics unchanged).

### Added
- `tests/test_deggroupcompare.sh`: regression test for DEGgroupCompare argument propagation (an Rscript shim records the arguments; runs with bash + python, no R/snakemake needed).
- New config key `batch_correction` (default `"F"`, matching the old behavior; with `"T"` DESeq2 runs the `~ batch + group` design and the sample table needs a third `batch` column).
- New config key `lncrna.pfam_scan` (default `"pfam_scan.pl"` via PATH lookup; pfam_scan.pl became the last external tool with a configurable path).

### Changed
- Comment fixes: stale v0.2-era file-name references left in `lncRNA_functions.sh` and `envs/{lncrna,qc,assembly}.yaml` (`config_lncRNA.yaml`, `RNA-seq_lncRNA_DenovoIdenti.smk`, `multiQC_cleaned`, etc.) were updated to the current `config/config.yaml` lncrna section and `rules/*.smk` naming.

## [0.4.0] - 2026-09-03

Stage 3 (sustainability: tests, CI, docs moving with code) completed; see the optimization roadmap document (archived since). The three-stage overhaul plan is now fully landed.

### Added
- `tests/make_testdata.py`: deterministic miniature test-data generator (reference/GTF/BED12/annotation for 2 x 100 kb chromosomes and 60 simulated genes; 2 groups x 2 samples of simulated PE reads — with sequencing errors, 10% adapter, 5% noise reads, plus an expected-differential-gene truth table); pure standard library, data not committed and generated on the fly on the test machine (3.1).
- `tests/check_outputs.py`: end-to-end output asserter — key file existence, mapping_stat consistency/numeric checks (complete IDs, Total_Reads = input read pairs, parseable alignment ratios), count matrix shape, DEG result direction checks (against the truth up/down genes), and plot/sessionInfo outputs (3.1).
- `tests/run_test.sh`: one-shot regression (generate data -> dry-run -> end-to-end -> assertions -> DAG regeneration, configurable via `--pipeline`/`--reads`/`--keep`) (3.1).
- `tests/lint.sh`: aggregated static checks — bash -n, shellcheck, Python compilation, R parsing, snakemake --lint (4 pipelines), auto-skipping missing tools (3.2).
- `.github/workflows/ci.yml`: GitHub Actions CI (lint job + miniconda end-to-end regression job + DAG artifact) (3.2).
- `workflow/profile/slurm/`: SLURM scheduling profile (3.3).
- `CONTRIBUTING.md`: contribution guide — branch/commit-message conventions, CHANGELOG requirements, documentation sync checklist, testing requirements (3.4).

### Changed
- `DEGgroupCompare.sh`: xargs gained `-r`, so fewer than 2 treatment groups (empty job list) no longer mis-executes.
- `envs/enrich.yaml`: added `bioconductor-org.hs.eg.db`; enrichment with species=hsa works out of the box.

## [0.3.0] - 2026-09-03

Stage 2 (structure normalization and quality improvement) completed; see the optimization roadmap document (archived since). **Note: the old entry points `rna-seq-workflow/RNA-seq_*.smk` are deprecated as of this version** (a deprecation stub remains); use `run.sh` / `workflow/Snakefile` instead.

### Added
- `workflow/Snakefile` as the single entry point: `pipeline=upstream|deg|as|lncrna` injected via the `RNASEQ_PIPELINE` environment variable (fixed at parse time), project configuration injected via `RNASEQ_CONFIG`; falls back to the repository default configuration when unset (for --lint and examples).
- Modularized `rules/`: `common` (helpers and STAR arguments) / `align` (trim/QC/STAR/strandedness) / `quant` (quantification merge) / `deg` / `as` / `lncrna`; the output root `results_dir` is configurable (default `results/`).
- Whole-pipeline MultiQC summary (FastQC/Trim Galore/STAR/featureCounts) + `workflow/multiqc_config.yaml` (P2-9).
- `software_versions.yaml`: each run automatically records environment definitions, the snakemake version, and the git commit (`collect_versions.py`); DESeq2 output includes `sessionInfo.txt` (P2-9).
- `config/` configuration system (P2-4/P1-11): `config.yaml` main configuration template (all keys commented, nested lncRNA section), `species.yaml` (osa/hsa species resource preset map), `samples.csv` template.
- `example/` example project: configuration template + a real project sample table (231107XTL, 27 samples, migrated from `rna-seq-workflow/sample_info.csv`, P2-3).
- `run.sh` moved to the repository root.

### Changed
- Unified the two trim rules (P2-2): exact detection of the dual raw fastq naming convention (`.fastq.gz`/`.fq.gz`) works across all pipelines; trim always runs FastQC (fixes the P1-1 empty MultiQC report); trim reports are renamed to canonical names keyed by sample id and declared as rule outputs (P1-6).
- STAR (P2-2): outputs `BAM SortedByCoordinate` directly (P1-4); the index is tracked at file level (P1-5); the assembly pipeline keeps assembly-optimized parameters and other pipelines use standard parameters, all extendable via `star_extra_args`.
- Plot fixes (P1-8): volcano plot dynamic coordinate ranges (no more clipping), `fontface`, explicit color mapping, tolerance for empty categories; PCA `ntop` configurable (`pca_ntop`, default 20000).
- Differential-result directory spelling fix: `Diff_Expr_Analysis_Reults` -> `Diff_Expr_Analysis_Results` (references in `enrich.sh` updated, P2-7/P2-8).
- Stabilized script naming (date/hyphen suffixes removed): `run_deseq2.R`, `run_enrichment.R`, `run_gsea.R`, `run_deg_compare.R`, `run_featurecounts.R`.
- The KEGG organism code is passed to the enrichment scripts via `species.yaml` (`kegg_organism`); `run_deg_compare.R` outputs now follow results_dir.
- lncRNA pipeline threads and external tool paths all read from config (the `lncrna` section, P1-7); key CPC2/CNCI/Pfam/NR products are declared as rule outputs (P1-6).
- `run_deseq2.R` keeps v0.2's exact-match control detection; the sample table `group` column is read by header name.

### Removed
- The old four-entry Snakefiles, old rule files (`rules/RNA-seq_*.smk`), old `config_*.yaml` (replaced by the `config/` system), and stale DAG figures (regenerated by stage 3 CI); deprecation stubs remain at the old entry locations.

## [0.2.0] - 2026-09-03

Stage 1 (P0 fixes) completed; see the optimization roadmap document (archived since).

### Added
- `envs/*.yaml`: 7 Conda environment definitions (qc / align / quant / deseq2 / enrich / assembly / lncrna), wired to rules via `conda:` directives (P0-2).
- `profile/sge`, `profile/default`: cluster and local scheduling profiles.
- `run.sh`: unified launcher — automatic scheduler selection, configuration priority resolution, pre-run sample table validation; replaces the deleted `main_runRNA-seq.sh`.
- `scripts/merge_featurecounts.py`: self-contained quantification merge script (replaces `featureCount.R_result_merge.sh` and an empty-shell `.py` that depended on out-of-repo `njoin.sh`/`transposition.sh`, P0-3/P0-4).
- `scripts/validate_samples.py`: sample table validator (duplicate ids, `-` characters, prefix conflicts, control group presence, layout values).
- `sample_info.csv` gained a `layout` column (PE/SE/auto).

### Changed
- Path parameterization (P0-1/P0-8): the entry Snakefile's `configfile:`/`include:` are based on `workflow.basedir`; in-rule script calls go through the `SCRIPTS` variable; R scripts dropped hardcoded `.libPaths`; machine-specific paths in configs became `/path/to/...` placeholders.
- `Mapping_stat` rewritten (P0-6): parses `*_Log.final.out` and trim reports per sample; sample ids are no longer truncated and paired-end samples no longer misalign.
- Raw fastq detection switched to exact file-name matching (no more `ls {sample}*` globbing, eliminating cross-prefix matching risk); `runSTAR --readFilesIn` uses explicit files; featureCounts strandedness/paired-end detection no longer relies on `ls` counting.
- Configurable control group (P0-7): new `control_group` config key, propagated through runDESeq2 (`-r`; control detection changed from `grepl` substring matching to exact matching) / enrich / getGroups.py / DEGgroupCompare.
- Configurable species (P0-8): new `species`/`annotation_tsv`/`orgdb_tarball` config keys; the enrichment R scripts support both `osa` and `hsa`; KEGG query failures skip the step instead of aborting the pipeline.
- Between-group comparison parallelization: `ParaFly` replaced with `xargs -P`; `getGroups.py` reads the `group` column by header name.
- lncRNA pipeline thread count configurable (`lncrna_threads`, formerly hardcoded 30); CPC2/CNCI/Pfam paths injected into `lncRNA_functions.sh` via config.
- `count_merge`/`count_merge2` explicitly declare the `GeneExpression_FPKM.xls`, `GeneCount_Assigned_logs.xls`, and `{sample}.log` inputs/outputs; `check_strandedness` declares the `_infer_experiment.out` output.
- `STAR_index` fixed its broken stderr redirection (`; 2>{log}` -> `>> {log} 2>&1`); `runSTAR` dropped the never-cleaned `--outReadsUnmapped Fastx` products.

### Removed
- Dead/broken scripts: `enrich_KEGG_clusterProfiler.R` (syntax error), `check_strandness.sh` (unexpanded template), `getExpr_featureCounts.sh` (truncated), `featureCount.R_result_merge.{sh,py}` (replaced by the new merge script), `enrich_GO_KEGG_clusterProfiler_gProfilerGO_hsa.R` (merged into the main script's species branch), `test.py`/`test.groups`, and the empty `SampleListFile`.

### Known issues
- The upstream/DEG pipeline trim rules did not run FastQC (MultiQC post-trim report was empty); KEGG required internet; the volcano plot coordinates were hardcoded, etc. — P1/P2 issues left for stage 2 (see the review report document, archived since).

## [0.1.0] - 2026-09-03

### Added
- Initialized the git repository, ingesting the `rna-seq-workflow/` pipeline code as-is as the baseline (four Snakemake entry points + rules + scripts + configuration).
- Standard project scaffolding:
  - `README.md` (project description and documentation index)
  - the user guide (operations manual: data preparation/configuration/run/results interpretation/FAQ/breakpoint list)
  - the review report (full code review: P0 x8 / P1 x14 / P2 x9 with file:line locations and methodology assessment)
  - the optimization roadmap (three-stage overhaul plan, acceptance criteria, option comparisons)
  - `CHANGELOG.md`, `LICENSE`
- `.gitattributes` (LF line endings for scripts, targeting Linux clusters), `.gitignore` (ignores Snakemake run outputs).

### Unchanged
- The pipeline code itself was not modified; all fix work proceeded from stage 1 of the optimization roadmap.

### Known issues
- See the review report document (archived since): hardcoded absolute paths, missing Conda environment files, missing count-merge external dependencies, empty-shell lncRNA pipeline end-point scripts, misaligned `mapping_stat.xls` statistics, hardcoded "control" group name, hardcoded enrichment species.
