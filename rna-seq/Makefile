# rna-seq development checks.
# make check  -- runtime resolver + R runtime propagation tests + shell syntax (no snakemake)
# make lint   -- static lint suite tests/lint.sh (bash -n/shellcheck/py_compile/R parse/
#                yaml/snakemake --lint; missing optional tools are skipped, CI runs all)
# make test   -- everything CI runs (= check + lint)
# Regression (needs snakemake, i.e. the CI dry-run job):
#   bash tests/run_test.sh --pipeline deg

PYTHON ?= python3

SHELL_SCRIPTS := run.sh $(wildcard workflow/scripts/*.sh)

.PHONY: check lint test

check:
	bash tests/test_runtime_config.sh
	bash tests/test_r_runtime_propagation.sh
	bash tests/test_deggroupcompare.sh
	@for f in $(SHELL_SCRIPTS); do bash -n $$f && echo "  PASS  bash -n $$f" || exit 1; done

lint:
	bash tests/lint.sh

test: check lint
