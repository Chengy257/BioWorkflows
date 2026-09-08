# Changelog

All notable changes to this project are documented in this file. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.6.0] - unreleased

### Added

- **QC gate summary stage (`qc.gates`, default off)**: one PASS/WARN/FAIL row per sample aggregating the QC metrics the workflow already computes — mapping rate (`samtools flagstat` of the analysis BAM, new per-sample `results/5.QC/gates/{sample}_flagstat.txt`), duplication rate (picard `PERCENT_DUPLICATION`), FRiP, NSC/RSC, TSS enrichment, and organelle fraction — into `results/5.QC/gates/gate_summary.tsv` (plus a MultiQC custom-content table). A metric whose source stage is off (or whose file is missing) renders `NA`; the per-sample gate is FAIL when any metric fails its threshold, WARN when nothing fails but at least one metric is NA, PASS otherwise (`failed`/`na` columns name the offending metrics). All thresholds are user-tunable numbers under `qc.gates.thresholds` (defaults: `mapping_rate_min` 0.70, `dup_rate_max` 0.50, `frip_min` 0.01, `nsc_min` 1.05, `rsc_min` 0.8, `tss_min` 6.0, `organelle_max` 0.20; validated as numbers in [0, 1] except `nsc_min`/`rsc_min`/`tss_min`, which accept any positive value) and are strictly informational — the pipeline never hard-fails on a gate. Aggregation lives in `workflow/scripts/gates_summary.py` (stdlib only, sample=group tokens baked from the sample table); the new `gates_flagstat`/`qc_gates` rules and their `config/resources.yaml` entries only enter the DAG when the stage is enabled (default dry-run baseline unchanged at 47 jobs; `tests/run_test.sh --gates` dry-runs the scenario, CI gained a matching step).
- **Spike-in normalization stage (`spike_in`, default off)**: second-pass alignment of the unmapped read pairs (the `--un-conc-gz` artifacts of the main bowtie2 step, now declared as `results/3.align/bowtie2/{sample}_unmapped.fq.{1,2}.gz`) against a spike-in genome (`spike_in.fasta`, required when enabled; a missing file warns only like every other reference) into `results/3.align/spike_in/{sample}_sorted.bam(.bai)`, with per-sample `samtools idxstats`/`flagstat` QC and a project-wide `results/5.QC/spike_in/Spikein_summary.tsv` (plus a MultiQC custom-content table under the configurable `spike_in.name` label) carrying `spike_total`, `spike_mapped`, `spike_rate`, and `scale_factor = 1e6 / max(spike_mapped, 1)`. With `spike_in.scale_bigwigs: true` the bigWig signal tracks are multiplied by the factors — per-sample tracks by their own factor, group FE tracks by the mean over the group's treat samples (with the stage off, the bamCoverage/bedGraphToBigWig commands are byte-identical to v0.5). The new `spike_bowtie2_index`/`spike_align`/`spike_summary` rules and their `config/resources.yaml` entries only enter the DAG when the stage is enabled (default dry-run baseline unchanged at 47 jobs; `tests/run_test.sh --spike-in` dry-runs the scenario with a synthetic 2-contig spike reference at 54 jobs, CI gained a matching step; summary math in `workflow/scripts/spikein_summary.py`, stdlib only).
- **Alternative pooled peak caller (`peak.caller: macs2 | seacr`, default `macs2`)**: SEACR (Sparse Enrichment Analysis for CUT&RUN, Yo et al. 2021) joins MACS2 as the pooled group caller for CUT&RUN/CUT&Tag-style data. With `peak.caller: seacr` the pooled MACS2 `callpeak_narrow`/`callpeak_broad`/`callpeak_atac` rules are replaced (parse-time include guard — the two caller modules never coexist in the DAG) by `seacr_bedgraph_treat`/`seacr_bedgraph_control` (pooled analysis BAMs → `samtools merge` → `bamCoverage --normalizeUsing None --binSize 1` raw-depth bedGraphs under `results/4.peak/seacr/`; controls only when the group has them) and `seacr_callpeak` (external `bash $CHIP_SEACR ...`, control bedGraph + `peak.seacr.normalize` when a control exists, numeric FDR threshold `peak.seacr.fdr_threshold` otherwise; `peak.seacr.mode` = stringent|relaxed). SEACR's 6-column output is converted to the standard 10-column narrowPeak/broadPeak contract by `workflow/scripts/seacr_to_narrowpeak.py` (stdlib: score = AUC capped at 0-1000, p/q = 0, `peak` = 0) at the same `{group}_peaks.*` paths MACS2 writes — FRiP, annotation, blacklist, motif, and DiffBind consumers are unchanged; `seacr_bigwig` keeps the `{group}_FE.bw` deliverable via the bedClip/sort/bedGraphToBigWig recipe (pooled raw depth, not bdgcmp FE). Limitations: the replicate/IDR stage stays MACS2-based, so `caller=seacr` + `peak.replicate.enabled=true` is an aggregated validation error, and `bigwig.per_sample` (a MACS2-module rule) is skipped with a parse-time warning under SEACR. SEACR resolves like idr/HOMER via the software.yaml `paths:` section (exported as `CHIP_SEACR`). Default dry-run baseline unchanged at 47 jobs; `tests/run_test.sh --seacr` dry-runs the scenario at 50 jobs and asserts the pooled MACS2 caller rules are absent while `seacr_callpeak` and the FRiP stage are present (real runs additionally need the external SEACR script); CI gained a matching step.
- **TOBIAS footprinting stage (`footprint`, default off)**: bias-corrected footprinting of the open-chromatin assays — the ATAC/FAIRE groups only (footprinting models Tn5 cut-site bias and is meaningless on chip/cuttag enrichment data; a table without any atac/faire group skips the stage with a parse-time warning). Per group: `footprint_ataccorrect` pools the group's treat analysis BAMs (`samtools merge`, dot-prefixed temp removed after the job), runs `TOBIAS ATACorrect` over `group_final_peak_file(group)` (IDR/consensus when the replicate stage is on, pooled otherwise) against `genome_fa`, and writes `{group}_corrected.bw` plus QC plots; `footprint_scorebigwig` computes footprint scores over the group's peak regions from the corrected track (TOBIAS 0.8 API `--signal/--regions/--output`, one `{group}_footprint_scores.bw`); `footprint_bindetect` (optional, `footprint.bindetect`, default on; `--skip-excel`, `footprint.motif_pvalue` → `--motif-pvalue`, empty = tool default) detects per-TF binding. Deliverables land under `results/7.footprint/{group}/`: the corrected bigWig is a declared file (snakemake cannot link inputs nested inside another rule's `directory()` output, and both downstream stages consume the track), while ScoreBigwig and BINDetect declare one directory output each (motif.smk precedent). TOBIAS resolves like idr/HOMER/SEACR via the software.yaml `paths:` section (exported as `CHIP_TOBIAS`; separate environment, e.g. `conda create -n chip-tobias -c bioconda tobias`). Validation: `footprint.enabled`/`footprint.bindetect` booleans, `footprint.motifs` a required non-empty string when enabled (a missing or `/path/to/`-style file warns only), `footprint.motif_pvalue` numeric in (0, 1] when set. Default dry-run baseline unchanged at 47 jobs; `tests/run_test.sh --footprint` dry-runs the scenario with a synthetic 2-motif JASPAR-style PFM at 50 jobs and asserts the three footprint rules schedule for the atac group g2 only (real runs additionally need the external TOBIAS install); CI gained a matching step.

### Fixed

- **WSL real-run validation round for the v0.5 replicate stage (2026-09-09)**: `tests/run_test.sh --replicate --real-run --reads 6000` now runs green end to end against the real `idr` 2.0.4 binary — the pairwise-IDR output column contract (idr appends per-peak IDR columns after the narrowPeak columns; the rules trim back to 10) and the broad-consensus behavior are confirmed on real outputs (13/18-row IDR final sets, 24-row broad consensus). Three findings folded back: (1) the synthetic enrichment design now seeds twenty-four 2kb peak regions on chr1 (was three) because `idr` refuses peak files with fewer than 20 peaks post-merge — the old data could never exercise a real idr run; (2) `annoPeak_batch.R` passes `xlim` to `tagHeatmap()` only when the installed ChIPseeker signature accepts it (the argument is required in ≤ 1.36 but removed in 1.38+, where the tagMatrix attribute carries the range); (3) `workflow/environment.yaml` pins `matplotlib<3.9` next to `deeptools=3.5.1` (matplotlib 3.9 removed `matplotlib.cm.register_cmap`, which deeptools imports), unpins `python` (macs2 2.2.7.1 solves on 3.10), and unpins the bioconductor-chipseeker/genomicfeatures/bedtools/ucsc builds whose exact versions left the current repodata (2026-09-09 solve round).

## [0.5.0] - unreleased

### Added

- **Replicate-aware peak stage (`peak.replicate`, default off — the v0.5 headline)**: pooled-replicate calling is no longer the only route. When enabled, every treat sample additionally gets its own MACS2 peak call at a relaxed narrow cutoff (`peak.replicate.qvalue`, default 0.01, ENCODE style) against the group's pooled control (`results/4.peak/replicates/{group}/{sample}_peaks.{narrowPeak,broadPeak}`); narrow groups with ≥2 treats get pairwise IDR via the classic `idr` tool (`results/4.peak/idr/{group}/{a}__vs__{b}.narrowPeak`) and a final reproducible peak set (`results/4.peak/{group}_IDR_peaks.narrowPeak` + `{group}_IDR_support.bed`); broad groups with ≥2 treats get a bedtools multiinter overlap consensus with a replicate-support cutoff (`results/4.peak/{group}_consensus_peaks.broadPeak` + support bed). Exactly-2-replicate groups use the single IDR pair as the final set; >2-replicate groups keep the union of pairwise IDR peaks with pairwise support ≥ `consensus_min_replicates` (a documented simplification of the ENCODE rescue/self-consistency scheme). A project-wide `results/5.QC/replicate_peaks/Replicate_summary.tsv` (per-group replicate counts, final count, retained fraction) is injected into MultiQC. `peak.replicate.frip_on` switches FRiP between the pooled set (default) and the reproducible set; peak annotation covers the final reproducible set. The `idr` binary is a python2 tool that deliberately stays OUT of the conda template: configure it via the software.yaml `paths:` section (exported as `CHIP_IDR`) or leave `idr` on PATH; it is only required when the stage is enabled.
- **TSS enrichment (`qc.tss`, default off)**: per atac/faire treat sample — strand-aware 1-bp TSS set derived from the gene-model BED (`scripts/tss_from_bed.py`), deeptools reference-point matrix over ±2kb (per-sample CPM coverage, 10bp bins), enrichment score = max of the baseline-normalized profile (`scripts/tss_score.py`); outputs under `results/5.QC/tss/` with a MultiQC table.
- **Organelle read fraction (`qc.organelle`, default off)**: `samtools idxstats` per sample plus a summary of the chloroplast/mitochondrial mapped-read fraction (`results/5.QC/organelle/`, MultiQC table) — the plant analog of the ATAC mitochondrial check. Contig patterns configurable via `qc.organelle_patterns` (patterns of ≤3 chars match contig names exactly, longer ones as substrings).
- **Artifact blacklist filtering (top-level `blacklist`, default empty = off)**: filtered copies of the pooled and final peak sets land in `results/4.peak/blacklist_filtered/` (bedtools intersect -v) with a before/after count table in `results/5.QC/blacklist/`; FRiP and peak annotation read the filtered copies. No ENCODE blacklist exists for rice — supply your own BED.
- **Differential binding (`diffbind`, default off — v0.5 Phase 3)**: DiffBind (DESeq2/edgeR) between explicit pairs of sample-table groups (`diffbind.contrasts`; each arm ≥ 2 treats, parse-time validated). The sample table gains OPTIONAL `condition`/`batch` columns (legacy 6-column tables unchanged): condition labels the contrast factor (group-name fallback), batch becomes a blocking factor. Per contrast, `scripts/diffbind_sheet.py` builds the DiffBind sample sheet (per-replicate peak files when the replicate stage is on — recommended; pooled group sets otherwise, with control attachment under `use_controls`) and `scripts/run_diffbind.R` counts over the consensus (summit-recentered via `summit_flank`, 0 = full regions), fits the contrast (optionally batch-blocked), and writes `results/6.diffbind/{A}__vs__{B}/` with `DB_results.tsv`, the `fdr`/`foldchange` filtered `DB_significant.tsv`, MA/volcano/PCA plots, and sessionInfo. `bioconductor-diffbind` joins the conda template (unpinned, R 4.3 line, mirroring the srna-seq DE shipping); the R preflight is unchanged.
- **Per-sample signal tracks (`bigwig.per_sample`, default off — v0.5 Phase 4)**: `bamCoverage` normalized coverage bigWigs per sample (`results/4.peak/samples/{sample}.bw`; RPGC via `genome_size` or CPM, configurable bin size) for browser-level replicate comparison, plus a `peak.bigwig_measure: FE | logFE` switch for the existing group track (filename stays `{group}_FE.bw`).
- **HOMER motif enrichment (`motif`, default off — v0.5 Phase 5)**: `findMotifsGenome.pl` on every group's final peak set into `results/6.motif/{group}/`; HOMER stays an external distribution resolved via the software.yaml `paths:` section (`homer_findmotifs` → `CHIP_HOMER_FINDMOTIFS`), `motif.homer_genome` (configureHomer tag or `custom:` path) is mandatory when enabled; size/background/extra pass through to the tool.
- **Scenario regressions in `tests/run_test.sh`**: `--replicate` (adds a 2-treat broad group to the synthetic data and enables the replicate stage), `--qc-full` (tss + organelle + synthetic blacklist), `--motif`, and `--diffbind` (adds a 2-treat narrow group with condition/batch columns and one contrast); the flags combine. The dry-run prints the planned job count (rule/localrule blocks): default 47 (unchanged), replicate 79, qc-full 60, motif 49, diffbind 69, all-on 134. CI runs each scenario as a separate step.
- `tests/run_tests.py`: extraction-backed checks for the new pure helpers (pair enumeration, IDR slug round-trip, narrow/broad group routing, organelle contig matcher), 12 new `validate_config` cases for the v0.5 keys, and subprocess tests for the four new stage scripts (tss_from_bed / tss_score / organelle_summary / replicate_summary).

### Changed

- FRiP and peak annotation resolve their peak inputs through `frip_peak_file()` / `annot_peak_files()` in common.smk; with every v0.5 switch off the resolved files are identical to the previous direct pooled paths, and the default dry-run baseline stays at 47 jobs.
- `workflow/scripts/runtime_config.py` (WorkflowSpec wrapper) exports `CHIP_IDR` from the software.yaml `paths:` section via the framework's `extra_exports` hook; the preflight is unchanged (idr is never demanded because it is not a pipeline tool).
- Sample-table/config validation aggregates the new keys (peak.replicate fields, blacklist string, qc.tss/qc.organelle booleans, organelle_patterns list, bigwig/motif/diffbind blocks) into the existing single error report; missing v0.5 blocks stay valid, so legacy project configs are unaffected. `load_sample_table` additionally accepts the optional condition/batch columns (consistency-checked within a group).

## [Unreleased]

### Added

- **`results_dir` output consolidation** (!): all derived artifacts now live under the configurable results root (default `results/`) in the project working directory — `results/0.index/bowtie2*` (bowtie2 index), `results/2.cleandata/`, `results/3.align/bowtie2/`, `results/4.peak/` (incl. `anno_result/`), `results/5.QC/` (incl. `spp/`, `frip/`, `deeptools/`, `software_versions.yaml`, `logs/`), and `results/logs/` (per-rule logs); raw inputs `1.rawdata/` stay at the working-directory root and the former top-level output directories (`0.index/`, `2.cleandata/`, ... in the workdir) no longer exist there.
- **Per-rule scheduler resources extracted to `config/resources.yaml`**: one dedicated file with per-rule `threads/mem_mb/runtime_min` (copy it into the project directory as `resources.yaml` and run.sh auto-detects it; the legacy top-level `threads` in config.yaml still acts as a global cap; `runtime_sec` is derived automatically). The `res()` helper is gone; rules read resources via the `rthreads`/`rmem`/`rruntime` helpers backed by centralized `RESOURCE_DEFAULTS`.
- **Species presets in `config/species.yaml`**: `osa`/`hsa` presets for genome_fa/gtf/bed/chromsize/genome_size, selected by `species: "osa"` in `config.yaml`; unset reference keys fall back to the preset, explicit config keys win.
- **Config layering completed**: repository `config/config.yaml` -> project config (`-c`) -> `config.local.yaml` in the project directory (auto-detected, or `-l/--extra-config`); `config.template.yaml` remains the annotated override template. Parse-time validation (sample table + config, aggregated errors with line numbers) is retained unchanged.
- **Shared-layer extraction**: run.sh now sources `../shared/lib/launcher.sh` (logging helpers) and auto-detects a project `resources.yaml` (exported as `CHIP_RESOURCES_CONFIG`); `workflow/scripts/collect_versions.py` became a thin wrapper over `shared/python/bioworkflows_versions.py` (the `runtime_config.py` half of the original claim only landed with the Changed entry below — it stayed standalone until then); `meta.smk` records the Snakemake version into `results/5.QC/software_versions.yaml`.
- **Root-level CI**: the CI workflow moved to the repository root `.github/workflows/ci.yml` (chip unit tests + lint + dry-run regression); the sub-level `.github/` directory was removed; the Makefile stays (`make check/lint/test`).

### Changed

- **Raw FASTQ naming variants supported**: `trim_adapter` no longer requires the literal `1.rawdata/{sample}_1.fq.gz` + `{sample}_2.fq.gz` — `raw_fastq_pair()` in `common.smk` resolves each sample's pair with priority `{id}_1/_2.fastq.gz` -> `{id}_1/_2.fq.gz` -> `{id}_R1/_R2.fastq.gz` -> `{id}_R1/_R2.fq.gz` (first complete pair wins; identical to the rna-seq workflow's resolution order; a missing pair raises a clear `WorkflowError` listing all supported patterns; the pipeline remains paired-end only), removing the need for manual renaming/symlinks at deployments.
- **R invocation routed through the configured R runtime**: `peak_annotation` now calls `{params.rscript}` (`RSCRIPT` from the `CHIP_RSCRIPT` environment, exported by `run.sh` from `software.yaml r.rscript`), instead of a bare `Rscript` that resolved to whatever R was on the cluster node's PATH (found on the enhancer_lncRNA_2026 deployment, 2026-09-05: the bare call picked up R 4.1.3 and crashed with `rlang.so: undefined symbol: EXTPTR_PROT`).
- **`runtime_config.py` migrated onto the shared WorkflowSpec framework**: the standalone 279-line resolver (ported from rna-seq v0.8.0 before the shared-layer extraction; it had never actually imported the framework despite the earlier "Shared-layer extraction" claim) is now a thin `WorkflowSpec` wrapper over `shared/python/bioworkflows_runtime.py` in the bs-seq style, closing `docs/TODO.md` §5 and completing shared-layer adoption across all five workflows. The `CHIP_*` export surface is byte-identical to the previous resolver (golden-diff verified), `check --scope all|software|r` semantics are preserved (`rscript` stays exported but is excluded from the `--scope software` preflight; the R branch stays always-on via non-empty `r_packages`; no OrgDb checks — annotation is GTF-based), and `run.sh` works unchanged (single-pipeline spec: no `--pipeline` needed).
- **`5.QC_deeptools/` renamed to `5.QC/deeptools/`**.
- **English unification**: all comments, user-facing messages, and documentation across the project were converted to English.
- **Relicensed from MIT to Apache-2.0**.

### Fixed

- P0: fixed the `res()` helper signature crash (resources now flow exclusively through the `rthreads`/`rmem`/`rruntime` helpers).
- Fixed the `config.local.yaml` overlay: the auto-detected project-local overlay is reliably appended last to the `--configfile` chain.
- `bigwig`: `bedGraphToBigWig` validates C-collation chromosome order (Chr1 < Chr10 < Chr11 < Chr12 < Chr2) regardless of the chrom.sizes line order, but the rule sorted with `bedtools sort -g` (which follows the chrom.sizes line order); whenever the two orders differ the rule fails at the Chr12 -> Chr2 boundary. Replaced with `LC_COLLATE=C sort -k1,1 -k2,2n` as the bedGraphToBigWig error message recommends (found on the enhancer_lncRNA_2026 deployment, 2026-09-05).
- `annoPeak_batch.R`: passed the required `xlim = c(-flank, flank)` to ChIPseeker `tagHeatmap()` (mirroring the preceding `plotAvgProf` call); without it the script aborted on the last PDF plot with `argument "xlim" is missing, with no default`, failing the whole `peak_annotation` rule after the first three plots had rendered (found on the enhancer_lncRNA_2026 deployment, 2026-09-05, job 304546).
- `run.sh --validate-only` works again under the pinned snakemake 7.32.4: `--list-rules` (the snakemake 8 name) was replaced with `--list`; under 7.x it failed with `unrecognized arguments: --list-rules`.

### TODO

- Minimal-sample end-to-end real run on a server (first real CI run + conda environment solving + MACS2 no-control `control_lambda` output confirmation; entry point `bash tests/run_test.sh --real-run`)
- Complete the DiffBind differential analysis (needs contrast/design-formula decisions; the stub scripts were archived out of the repository as of v0.4.0)
- bowtie2 `.bt2l` large-genome index support
- executor-plugin style profiles for snakemake 8.x (`snakemake-executor-plugin-cluster-generic`)

## [0.4.0] - 2026-09-04

Full alignment with the sister project rna-seq (v0.8.0) engineering system: directory layout, environment management, launcher experience, resource model, tests, and documentation. The design document and implementation plan for this alignment were archived out of the repository.

### Added (aligned with rna-seq v0.8.0)

- **Unified environment system trio**: `workflow/environment.yaml` (pinned all-in-one main environment merging every version constraint of the former 11 per-rule envs), `config/software.yaml` (environment.type=system/conda_prefix/conda_name + R runtime + tool path overrides), `workflow/scripts/runtime_config.py` (`export` injects `CHIP_*` environment variables / `check` preflight subcommands)
- **software_versions rule + collect_versions.py**: records the actually-used tool versions, git commit, and runtime mode into `5.QC/software_versions.yaml` at run time (the basis for reproducibility audits and methods sections)
- **Full run.sh ops CLI** (535 lines, replaces main_run.sh): `--profile auto|default|pbs|sge|slurm` (auto detection: sbatch→slurm; qsub disambiguated to PBS/SGE via SGE_ROOT), per-rule resource placeholder cluster submit strings, `--memory/--runtime/--queue/--partition` overrides, `--check-software/--check-r` preflight, `--validate-only`, `--unlock`, `--retries`, `--log` tee + trap timing, `--` passthrough, `-r` raw-data renaming, config.local.yaml auto-layering
- **Per-rule resource model**: all 22 rules declare `mem_mb/runtime_min/runtime_sec` (4 rules also gained threads:1); `config.yaml` gained a `resources:` override section (per-rule overrides); the `res()` helper supports project-level tuning
- **Four profiles**: `workflow/profile/{default,pbs,sge,slurm}/config.yaml` (unified resource placeholders; pbs walltime in seconds to avoid format ambiguity)
- **MultiQC customization**: `workflow/multiqc_config.yaml` (title/workflow identity/intro), wired into the multiqc rule via `-c`
- **Tests and CI**: `tests/lint.sh` (six-stage static checks, missing tools skipped automatically); `tests/make_testdata.py` (deterministic synthetic-data generator: 2×100kb reference genome + chr1 three-peak-region enriched chip/atac PE reads + sample table + test config, fixed seed, byte-identical reproducibility); `tests/run_test.sh` (synthetic-data dry-run regression by default + `--real-run` server switch + output assertions); CI rewritten as two jobs (lint + dry-run regression, the latter needing no conda environment)
- **Four docs**: `docs/user-guide.md` (8-chapter operations manual), `CONTRIBUTING.md`, a full README rewrite (keeping the mermaid diagram/QC threshold table/results quick-reference), and the `example/` real-project template (2-sample table + project config + one-command guide)

### Changed (breaking)

- **Directory layout moved to the Snakemake standard**: `workflow.smk` → `workflow/Snakefile` (pure orchestration); ~200 lines of shared definitions (sample-table parsing/config validation/query helpers/target aggregation) split into `workflow/rules/common.smk`; rules and the R script moved into `workflow/{rules,scripts}/`; sample-table template → `config/samples.csv`; profiles → `workflow/profile/` (all via git mv, preserving history)
- **Launcher replaced**: `main_run.sh` removed; use `bash run.sh` uniformly (positional argument = working directory, `-P/-w` both accepted); the old `-b/-e/-E` conda deployment options were dropped along with the environment-route switch
- **Sample-table parse-failure behavior**: when the `grouplist` three-level resolution (absolute > working directory > repository) finds no file, it now errors out (no more silent fallback to the repository-root example file); the default value now points at the `config/samples.csv` template
- The 5 standalone QC stub scripts not wired into the DAG (DiffBind/ChIPQC/DROMPAplus etc.) were archived out of the repository (an archive mapping table accompanied the move)
- Unit tests 45 → **55 checks** (removed 1 envs-integrity check; added 5 resource-declaration checks + 6 synthetic-data generator checks)

### Removed

- **Per-rule conda system** (!): the 11 environment files under `envs/` and the 19 `conda:` directives inside rules were all removed in favor of the unified environment route; servers need a one-time `mamba env create -f workflow/environment.yaml`, or reuse an existing environment via `config/software.yaml` (three onboarding paths in user-guide §1)

### Fixed

- Tool-name mapping fidelity: preflight/version recording map deeptools → `bamCoverage` (representative binary) and spp → `run_spp.R` (matching the actual spp_qc.smk invocation, eliminating false preflight failures)
- The meta.smk Python interpreter is injected via `CHIP_PYTHON` (aligned with rna-seq); PBS `.o` log collection gained a profile guard (no longer swallows slurm outputs)

## [0.3.0] - 2026-09-03

### Added (ops review P3 batch)

- **FRiP / NSC-RSC injected into the MultiQC report**: the frip_summary and the new spp_summary rules emit `_mqc.tsv` custom-content tables that the multiqc aggregation rule picks up automatically — QC metrics (fastqc + bowtie2 + picard + FRiP + NSC/RSC) concentrated in a single report; QC switch state is wired to conditional includes
- **Analysis window parameter unified**: the new `region_flank` (default 3000) controls the ChIPseeker flank/TSS window and the deeptools computeMatrix up/downstream length from one place (previously three hard-coded spots); `annoPeak_batch.R` gained a fourth parameter
- **`profiles/pbs/`**: snakemake 7.x PBS profile (cluster parameters pinned into the repository: `{rule}` job name, `{threads}` cores, latency-wait/rerun-incomplete defaults); 8.x users keep using main_run.sh (version auto-adaptation), with the migration direction noted in the README
- 4 new tests (region_flank validation ×2 + **mqc shell rule-body execution tests** ×2 — decoded via the snakemake-equivalent `ast.literal_eval` chain + format rendering + actual bash execution with output-format assertions), for **45 checks** all passing

### Changed

- SPP rule dropped `-savp` (the pdf side-product filename derives from the input BAM and lands in the cwd, so it cannot be declared; all metrics are already in the `-out` text)
- `envs/bigwig.yaml`: ucsc-bedclip/bedgraphtobigwig pinned to bioconda build 482 (removing drift risk from semantically versionless packages)
- `envs/trim-galore.yaml`: removed the redundant explicit `cutadapt=4.4` pin (pulled in by the trim-galore dependency itself; double pinning increases solver conflict surface)

## [0.2.2] - 2026-09-03

### Added (ops review P2 batch)

- **`main_run.sh` cluster robustness**: `--rerun-incomplete` and `--latency-wait` enabled by default (`-t` tunable, default 90s), covering the two common false-failure classes of PBS resume runs and shared-filesystem output visibility delay
- **Shared conda environment directory**: new `-e DIR` (`--conda-prefix`) to reuse one environment set across projects; new `-E` pre-build mode (`--conda-create-envs-only`) to create environments on the login node when PBS compute nodes lack internet
- **trim_galore parameters made configurable**: four new keys `trim.quality/stringency/error_rate/extra` (previously hard-coded `-q 25 --stringency 3 -e 0.1`), template updated in sync
- **Centralized config validation** (`workflow.smk` `validate_config`): required keys/sub-keys/types/value ranges aggregated into one report (including a friendly error for non-integer threads); missing reference files only warn and do not abort (keeping --lint/dry-run parseable on machines without the reference files)
- 10 new tests (validate_config executed from real extracted source), for **45 checks** all passing

### Changed

- **callpeak three rules' threads lowered to 1**: MACS2 is single-threaded; previously requesting cluster resources at config threads=12 caused over-subscription/waste; combined with v0.2.1's `ncpus={threads}`, PBS requests now match actual usage exactly

## [0.2.1] - 2026-09-03

### Fixed (2 P1 items from the ops review)

- **PBS resource parameters decoupled from workflow thread counts**: the `main_run.sh`/README cluster submit examples changed to `-l ncpus={threads}` — snakemake formats the cluster string with each job's actual thread count, aligning automatically with config `threads` (the old example's fixed `ncpus=6` contradicted the default `threads: 12`, causing over-subscription)
- **snakemake version matrix unclear**: `main_run.sh` probes the snakemake major version at startup — ≥8 automatically uses `--software-deployment-method conda` (`--use-conda` is deprecated in 8.x), 7.x keeps `--use-conda`; a clear error is emitted when snakemake is missing or its version cannot be parsed. The README gained a version support matrix (7.32.4 reference version / 8.x auto-adapted / <7 and ≥9 unverified), with notes added to the manual commands

## [0.2.0] - 2026-09-03

> The v0.1.0 original implementation was archived out of the repository in its entirety. The refactor went through an independent code review (fix-first); the 3 P1 findings it produced plus all P2/P3 findings were fixed within the same version.

### Added

- Unified entry `workflow.smk`: new sample-table schema (sample_id/role/group/seqtype/layout/peak_type), row-level validation (line-number errors, name character-set checks), per-assay automatic routing; the old 4 entry points archived
- Brand-new rule set `rules/`: upstream (trim/fastqc/multiqc/bowtie2, fastqc parallel per sample), dedup (picard, per-assay switch, CUT&Tag keeps duplicates by default), callpeak (narrow/broad/atac three rules parallel per group + bigwig), annotation (ChIPseeker batch annotation), frip (FRiP + summary table), qc_deeptools (correlation/PCA/fingerprint/fragment size/gene-region signal, the full set), spp_qc (optional NSC/RSC)
- Per-rule conda environments `envs/` (11, conda-forge+bioconda only), replacing the illegal `conda: "chip"` usage
- Config system: complete config schema (genome_size/min_mapq/per-assay dedup/peak thresholds/qc switches) + `config.template.yaml` + the `config.local.yaml` layering mechanism
- Sensible default thresholds: narrow q=0.05, broad_cutoff=0.05, min_mapq=30 (common ENCODE values), all configurable
- ATAC/FAIRE peak-calling dual mode: `bampe` (default, ENCODE ATAC v2 recipe) | `shifted` (classic Tn5 offset recipe)
- `main_run.sh` rewritten: getopts parameterization (incl. dry-run pre-check, optional PBS, config.local auto-layering, correct cluster/local -j/--cores split)
- Tests and CI: `tests/run_tests.py` (31 dependency-free unit tests executing real source extracted from workflow.smk), Makefile (check/lint/dryrun), GitHub Actions (tests + shellcheck + snakemake --lint)
- MIT License; an independent-review fix-status table and an archive-mapping README (both since archived out of the repository)

### Fixed (all v0.1.0 P0/P1 findings from the independent review)

- Peak-calling rules wired into the DAG (the include was broken and all rule-all targets were commented out)
- Dedup rule shell command mistakenly placed in a conda block, illegal `${id}` wildcard, threads as a string type
- Peak-calling loop variable unused, causing each group to be called twice; conflicting empty-rule outputs
- Sample-table schema contradicted its consumers (three mutually incompatible definitions coexisted)
- `get_samples()` deduplicated samples by the seqtype column — wrong logic
- Annotation rule referenced a non-existent script path
- All 16 machine-specific absolute paths (/home/chengyu, /opt, /share)
- `run_ChIPQC.R`: args used before definition, `basename()` missing an argument, tutorial leftover names, private paths
- bigwig chromosome lexicographic order mismatched chrom.sizes order (switched to `bedtools sort -g`)
- r-chipseeker conda environment R 4.2 conflicted with Bioc 3.18 package versions (r-base=4.3)
- ATAC shift/extsize silently ignored in BAMPE mode (switched to the dual-mode toggle)
- MACS2 loose thresholds like `-q 0.5` (changed to standard defaults and made configurable)
- Duplicate call_peak.sh at the repository root and in scripts/

### Changed

- README fully rewritten (new usage/new schema/QC thresholds/known limitations); line-ending policy enforced to LF via `.gitattributes`

## [0.1.0] - 2024-03

### Added (original snapshot, commit 03f5c0c, since archived out of the repository)

- ChIP-seq / CUT&Tag / ATAC-seq / FAIRE-seq four-entry Snakemake workflow
- Upstream rules: trim_galore, FastQC, MultiQC, bowtie2 alignment
- Peak-calling scripts: SPP fragment-length estimation + MACS2 (narrow/broad/ATAC modes), HMMRATAC, FSeq2
- MACS2 bdgcmp → bigWig conversion scripts
- ChIPseeker batch/single-sample peak annotation R scripts
- deeptools / ChIPQC / DROMPAplus QC scripts (some partial or empty stubs)
- Legacy full-environment export `chip_environment.yaml` and the PBS launcher `main_run.sh`
