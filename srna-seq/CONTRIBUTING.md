# Contributing guide (CONTRIBUTING)

## 1. Development workflow

1. **Branches**: keep `master`/`main` runnable; develop features and fixes on feature branches (`feat/<topic>`, `fix/<topic>`) and merge after self-testing passes.
2. **Commit messages**: `<type>: <English summary>` (scope optional, e.g. `feat(launcher): ...`; add `!` for breaking changes, e.g. `feat(launcher)!: ...`); keep types consistent with the existing repository history:
   - `feat` new feature / `fix` bug fix / `refactor` restructuring (no behavior change) / `docs` documentation / `test` tests / `ci` CI / `chore` misc
   - One commit does one thing; commit behavior changes and pure refactors separately (easier to trace later).
3. **Line endings and artifacts**: scripts use LF uniformly (enforced by `.gitattributes`); run outputs and test data are never committed (`.gitignore`; `tests/run_test.sh` cleans up after itself by default — do not commit `--keep` artifacts).

## 2. CHANGELOG requirements

- **Any change that affects user behavior** (adding/modifying/removing config keys or `software.yaml`/`resources.yaml`/`species.yaml` keys, `run.sh` options, rule inputs/outputs, script arguments, how things run, dependency pins) must be recorded in the Unreleased section of `CHANGELOG.md`, following the [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) format (Added/Changed/Fixed/Removed/Deprecated sections).
- Breaking changes (`!`) must describe the migration path in their entry; pure documentation or comment-level changes may be merged into a single entry at your discretion.
- On release, turn the Unreleased section into a version number + date and cut a semantic version tag.

## 3. Documentation sync checklist (check whenever touching rules)

Before modifying `workflow/`, `run.sh`, or `config/`, verify each item:

- [ ] Added/modified **config keys**? → sync the `config/config.yaml` comments, `config/config.template.yaml`, and the full key table in `docs/user-guide.md` §4 (plus the §4.3 defaults table when cluster resources are involved)
- [ ] Added/modified **per-rule resource defaults**? → sync `config/resources.yaml` and the defaults table in `docs/user-guide.md` §4.3
- [ ] Added/modified **output file paths or directories**? → sync the directory tree and results quick-reference in `docs/user-guide.md` §7 and the real-run assertions in `tests/run_test.sh`
- [ ] Added/modified **sample-table columns or validation rules**? → sync the `config/samples.csv` template and `docs/user-guide.md` §3
- [ ] Modified **`run.sh` options or behavior**? → sync the usage text inside `run.sh` and `docs/user-guide.md` §2/§5/§6
- [ ] Added/modified **software dependencies**? → sync `workflow/environment.yaml` (pinned) and `config/software.yaml`; do not reintroduce standalone `conda:` directives in rules (external-runtime system since v0.1.0)
- [ ] Does it break **compatibility with outputs/usage of older versions**? → describe the migration path in a prominent CHANGELOG entry; update the README version matrix if necessary

## 4. Testing requirements

Before submitting changes that affect workflow logic, run what the environment allows:

```bash
make check                            # bash -n syntax checks (no snakemake needed)
make lint                             # static suite tests/lint.sh (bash/shellcheck/py_compile/yaml/snakemake --lint; missing optional tools are skipped)
make test                             # CI-equivalent full check (= check + lint)
bash tests/run_test.sh                # synthetic-data dry-run regression (needs snakemake + python3/PyYAML; CI passes --reads 2000)
bash tests/run_test.sh --real-run     # end-to-end run + output assertions (needs a full analysis environment; server validation)
```

- The local minimum bar is `make check`; PRs touching Snakefile/rules/config must additionally run `bash tests/run_test.sh` in an environment with snakemake.
- `--real-run` needs a full environment with bowtie/trim-galore/fastqc/multiqc (i.e. the same as `workflow/environment.yaml`); use it when validating deployments on servers.
- When adding/modifying output files, add matching entries to the real-run assertions in `tests/run_test.sh`.
- CI (GitHub Actions, repository-root `.github/workflows/ci.yml`) runs the lint suite + dry-run regression (`bash tests/run_test.sh --reads 2000`); when a CI environment is not available locally, the commands above are the reference.

## 5. Code style

- **Naming**: rules, config keys, scripts, and functions use `snake_case` uniformly; sample/group names follow the sample-table character-set constraints (alphanumeric plus `. _ -`).
- **Centralized validation**: config and sample-table validation run only at parse time (`validate_config` / `load_sample_table` in `workflow/rules/common.smk`), with error messages **aggregated in one report and carrying line numbers**; do not scatter validation into rule shells.
- **Language**: comments, config templates, and user-facing messages are English throughout (project-wide English unification); keep the style consistent within each file.
- **Snakemake**: keep the same indentation and directive order as existing rules (input → output → params/log → threads/resources → shell); cluster resources go through the `rthreads()`/`rmem()`/`rruntime()` helpers backed by `config/resources.yaml` (with the project-level `resources:` override and the legacy top-level `threads` cap), and memory/duration must not be hard-coded in shell strings.
- **Shell**: start from `set -Eeuo pipefail` (maintenance scripts at least `set -euo pipefail`), quote paths, self-resolve `REPO_DIR`; prefer portable constructs (GNU coreutils first).
- **Python**: standard library first (the test generator/validators must not pull in third-party dependencies; `tests/make_testdata.py` must keep a fixed seed and byte-identical reproducibility).
- **Formatting**: 4-space indentation (smk/py/sh), 2-space (yaml/md), tabs in Makefiles (agreed in `.editorconfig`).
