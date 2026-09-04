# chip_cuttag_atac_faire development check entry points
# make check  -- unit tests + shell syntax checks (no snakemake needed)
# make lint   -- unified static check suite tests/lint.sh (bash -n/shellcheck/
#                py_compile/R parse/yaml/snakemake --lint; missing optional tools
#                are skipped automatically; CI runs everything)
# make test   -- the CI-equivalent full check (= check + lint)
# Regression tests (need snakemake; same as the CI dry-run job):
#   bash tests/run_test.sh             # synthetic-data dry-run, validates DAG integrity
#   bash tests/run_test.sh --real-run  # end-to-end run + output assertions (server validation)

PYTHON ?= python3

SHELL_SCRIPTS := run.sh $(wildcard workflow/scripts/*.sh)

.PHONY: check lint test

check:
	$(PYTHON) tests/run_tests.py
	@for f in $(SHELL_SCRIPTS); do bash -n $$f && echo "  PASS  bash -n $$f" || exit 1; done

lint:
	bash tests/lint.sh

test: check lint
