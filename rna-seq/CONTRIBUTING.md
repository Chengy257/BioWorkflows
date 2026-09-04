# Contribution Guide (CONTRIBUTING)

## 1. Development process

1. **Branches**: `main` stays runnable; feature development and fixes happen on feature branches (`feat/<topic>`, `fix/<topic>`) and merge after self-testing passes.
2. **Commit messages**: `<type>: <summary>`, with types consistent with the existing repository history:
   - `feat` new feature / `fix` bug fix / `refactor` structural change (no behavior change) / `docs` documentation / `chore` miscellaneous
   - One commit does one thing; behavior changes and pure refactors are committed separately (easier to trace).
3. **Line endings**: scripts use LF uniformly (enforced by `.gitattributes`); run outputs are never committed (`.gitignore`).

## 2. CHANGELOG requirements

- **Any change that affects user behavior** (added/changed/removed config keys, rule input/output paths, script arguments, run procedure) must be recorded in the unreleased section of `CHANGELOG.md`, following the Keep a Changelog format.
- Purely documentation or comment-level changes may be recorded together at your discretion.

## 3. Documentation sync checklist (check for every rule change)

Before modifying `workflow/rules/` or `workflow/Snakefile`, verify each item:

- [ ] Added/changed **config keys**? -> sync the `config/config.yaml` comment template + the §4.3 configuration table in `docs/user-guide.md`
- [ ] Added/changed **output file paths or directories**? -> sync §6 results layout in `docs/user-guide.md` + the assertions in `tests/check_outputs.py`
- [ ] Changed **rule inputs or script call arguments**? -> sync the Usage comments of the corresponding `workflow/scripts/` script + the user guide FAQ
- [ ] Added/changed **software or R runtime dependencies**? -> sync `config/software.yaml`, `workflow/environment.yaml`, and the user guide; do not reintroduce standalone `conda:` directives in rules
- [ ] Changed **R package / Rscript / R library** requirements? -> update the runtime preflight and `tests/test_runtime_config.sh`
- [ ] Does it break **compatibility with older outputs**? -> describe the migration path in a prominent CHANGELOG entry

## 4. Testing requirements

Before submitting changes that affect pipeline logic, run on Linux:

```bash
bash tests/lint.sh                          # static checks (bash/shellcheck/py/R/snakemake --lint)
bash tests/test_runtime_config.sh           # unified software/R runtime resolution
bash tests/test_r_runtime_propagation.sh    # nested R call inheritance checks
bash tests/run_test.sh --pipeline deg       # end-to-end regression (generate data -> dry-run -> run -> assert)
```

- `deg` is the default acceptance pipeline (covers all align+quant+deg rules); when touching alignment/assembly-related rules, also run `--pipeline as`.
- When adding output files, add the corresponding assertions in `tests/check_outputs.py`.
- `make check` / `make lint` / `make test` (in this directory) wrap the same checks; CI (GitHub Actions, repository-root `.github/workflows/ci.yml`) runs lint and a fast `--reads 20000` regression automatically on push/PR.

## 5. Code style

- Snakefile/rules: keep the same indentation and directive order as existing rules (input -> output -> log/params -> threads/resources -> shell); `{}` escaping inside shell strings follows `{{}}`.
- Shell scripts: start with `set -u`, quote paths, self-resolving `SCRIPT_DIR`; prefer portable constructs (GNU coreutils first).
- Python: standard library first (test generators/validators must not add third-party dependencies); R: keep the existing `getopt` style.
