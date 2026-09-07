# shared/

Cross-project code shared by the workflows in this monorepo. Every workflow
consumes it in place (nothing is copied), so a fix here lands everywhere at
once.

| Path | Consumed by | Purpose |
|---|---|---|
| `python/bioworkflows_runtime.py` | all five workflows — each via its thin `workflow/scripts/runtime_config.py` `WorkflowSpec` wrapper | Unified software-runtime framework: resolves the main environment (system/Conda), per-tool executables, the R runtime (Rscript / version policy / library paths), and external databases from `config/software.yaml`; implements the `export` (env injection for `run.sh`) and `check` (preflight) subcommands. Projects declare a `WorkflowSpec` (tool tables, pipelines, R packages, extra checks) in a thin wrapper. |
| `python/bioworkflows_versions.py` | all five workflows — each `workflow/scripts/collect_versions.py` | Captures the actually-resolved runtime (environment, Conda prefix, Rscript/R version/R libraries, per-tool paths, databases, git commit, Snakemake version) into `results/**/software_versions.yaml` for reproducibility audits. |
| `lib/launcher.sh` | all five `run.sh` | Launcher helpers: `timestamp` / `info` / `warn` / `die` logging and `resolve_path`. Sourced via `BIO_WORKFLOWS_SHARED` (exported by `run.sh`; the Python wrappers fall back to a path relative to the script location). |

Deliberately **not** shared (yet): scheduler auto-detection and cluster submit
strings inside each `run.sh`. All five launchers now embed the same
qsub-disambiguation policy (an installed `SGE_ROOT` means SGE, otherwise PBS),
but extracting the logic into this layer is a behavior decision — unify only
after one is deliberately made.

## Conventions

- Python modules import via `sys.path.insert` on `$BIO_WORKFLOWS_SHARED/python`
  with a relative fallback (`workflow/scripts → ../../..`), so direct
  execution keeps working when the env var is not set.
- Shell libraries are sourced, never executed.
- This layer must stay project-agnostic: anything project-specific belongs in
  the project's wrapper (tables, checks, exports), not here.
