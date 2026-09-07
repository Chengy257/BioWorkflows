# chip_cuttag_atac_faire follow-up items (TODO)

> Produced during the enhancer_lncRNA_2026 real-data preparation (2026-09-05).
> Status update (2026-09-05, branch fix/realdata-todos): items 1-4 below are all
> resolved at the repository level; the notes under each item record what server
> deployments can retire after pulling this branch. Item 5 was closed on
> 2026-09-08 (runtime_config.py migrated onto the shared WorkflowSpec framework);
> all five items are now resolved at the repository level.

## 1. FASTQ naming-variant support — DONE

- **Original issue**: the `workflow/rules/upstream.smk` fastq inputs only matched
  the literal `1.rawdata/{sample}_1.fq.gz` + `{sample}_2.fq.gz`; the
  `{sample}_R1/_R2.fq.gz` and `.fastq.gz` suffix variants commonly delivered by
  sequencing vendors were not recognized.
- **Fix**: `raw_fastq_pair()` in `workflow/rules/common.smk` resolves each
  sample's raw pair with priority `{id}_1/_2.fastq.gz` -> `{id}_1/_2.fq.gz` ->
  `{id}_R1/_R2.fastq.gz` -> `{id}_R1/_R2.fq.gz` (first complete pair wins;
  identical to the rna-seq workflow's resolution order);
  `trim_adapter` consumes it via input functions, and a missing pair raises a
  clear `WorkflowError` listing every supported pattern (the pipeline remains
  paired-end only).
- **Note**: the server deployment's workaround — symlinking a
  `{sample}_1.fq.gz`/`{sample}_2.fq.gz` set inside the project `1.rawdata/` —
  can be retired after pulling this branch.

## 2. `--validate-only` incompatible with the pinned Snakemake 7 — DONE

- **Original issue**: `run.sh --validate-only` invoked `snakemake --list-rules`,
  which is the Snakemake 8 rename (7.x uses `--list`); under the pinned
  snakemake-minimal=7.32.4 it failed with `unrecognized arguments: --list-rules`.
- **Fix**: `run.sh` now calls `--list` (help text and log messages updated to
  match); validation keeps working as a full parse + rule listing.

## 3. Rules calling bare `Rscript` bypass the `r.rscript` config — DONE

- **Original issue**: rules called `Rscript ...` directly (the `peak_annotation`
  rule / `annoPeak_batch.R`), never going through the `software.yaml`
  `r.rscript` / `CHIP_RSCRIPT` channel. When a project config pointed at a
  standalone R (e.g. an R 4.2.3 wrapper script), the bare call inside the job
  still resolved to the main environment's R (R 4.1.3 + newer R libraries on the
  enhancer_lncRNA_2026 deployment -> `rlang.so: undefined symbol: EXTPTR_PROT`,
  `peak_annotation` failure).
- **Fix (repo level)**: the RSCRIPT channel — `common.smk` defines
  `RSCRIPT = os.environ.get("CHIP_RSCRIPT", "Rscript")` (exported by `run.sh`
  from `software.yaml r.rscript` via `runtime_config.py`, resolved at
  orchestrator parse time and baked into the jobscript), and `annotation.smk`
  invokes `{params.rscript}` instead of a bare `Rscript`.
- **Notes**: the project `bin/env.sh` (BASH_ENV PATH-injection) workaround can
  stay in place as belt-and-suspenders; the symlink workaround used in the
  server deployment can be retired after pulling this branch.
- **Real-run validation (2026-09-05)**: the env.sh workaround proved the
  analysis itself is sound — R/ChIPseeker loaded, TxDb built, all 3 narrowPeak
  files annotated, and plotAnnoBar/plotDistToTSS/plotAvgProf all rendered.

## 4. `annoPeak_batch.R` `tagHeatmap()` call missing the required `xlim` — DONE

- **Original issue**: `workflow/scripts/annoPeak_batch.R` called
  `tagHeatmap(tagMatrixList)`. ChIPseeker's `tagHeatmap` has a required `xlim`
  argument with no default -> `Error in peakHeatmap.internal2(...): argument
  "xlim" is missing, with no default`; the script aborted on the last PDF plot
  (TSS tag heatmap) with a non-zero exit code, failing the whole
  `peak_annotation` rule. Reproduced on the enhancer_lncRNA_2026 real run
  (2026-09-05, job 304546, 13:19-13:23: the first three plots rendered fine and
  only tagHeatmap failed — a pure script bug, environment-independent).
- **Fix**: the call now passes `xlim = c(-flank, flank)`, consistent with the
  preceding `plotAvgProf(..., xlim = c(-flank, flank))` line (flank=3000, i.e.
  the TSS +/-3kb window).

## 5. Migrate `runtime_config.py` onto the shared WorkflowSpec framework — DONE

- **Original issue**: chip was the only workflow whose
  `workflow/scripts/runtime_config.py` (279 lines, ported from rna-seq v0.8.0
  before the shared-layer extraction) did not import
  `shared/python/bioworkflows_runtime.py`; the other four workflows declare a
  thin `WorkflowSpec` wrapper (43-95 lines). chip already consumed the other
  two shared pieces (`bioworkflows_versions.py` via `collect_versions.py`,
  `lib/launcher.sh` via `run.sh`).
- **Fix**: the standalone resolver was replaced by a thin `WorkflowSpec`
  wrapper in the bs-seq style (56 lines):
  - the 14 `DEFAULT_TOOLS` entries keep their exact insertion order, so the
    emitted `CHIP_TOOL_*` export lines stay byte-identical;
  - `rscript` stays in `DEFAULT_TOOLS` (the framework resolves it
    unconditionally and `export` emits `CHIP_RSCRIPT`/`CHIP_TOOL_RSCRIPT`)
    but is excluded from `pipeline_tools`, so `--scope software` still skips
    the R binary check while the r scope covers Rscript, its version, and the
    packages;
  - `r_packages = {"default": ["GenomicFeatures", "ChIPseeker"]}` keeps the
    shared framework's R branch always active under `--scope r` / `--scope all`;
  - `orgdb` stays empty: peak annotation is GTF-based (`makeTxDbFromGFF` in
    `annoPeak_batch.R`), so no OrgDb checks exist (the "orgdb/annotation
    checks" concern in the original deferral rationale was overstated).
  `run.sh` and the rules work unchanged: the single-key `pipeline_tools` means
  `check` needs no `--pipeline`, and `--analysis-config` is still accepted as
  the reserved extension point.
- **Validation**: `export --config config/software.yaml` output diffed
  byte-identical against a pre-migration golden capture; `tests/run_tests.py`
  all PASS; `tests/lint.sh` + `tests/run_test.sh` green with the 47-job
  dry-run DAG unchanged.
- **Reference**: the consumers matrix in `shared/README.md`; the root
  `AGENTS.md`.

