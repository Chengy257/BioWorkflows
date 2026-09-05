# chip_cuttag_atac_faire follow-up items (TODO)

> Produced during the enhancer_lncRNA_2026 real-data preparation (2026-09-05).
> Status update (2026-09-05, branch fix/realdata-todos): items 1-4 below are all
> resolved at the repository level; the notes under each item record what server
> deployments can retire after pulling this branch.

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
