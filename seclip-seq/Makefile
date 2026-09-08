# seclip-seq development check entry points
# make check  -- shell syntax checks (no snakemake needed)
# make lint   -- unified static check suite tests/lint.sh (missing optional
#                tools are skipped automatically; CI runs everything)
# make unit   -- stdlib unit tests for the workflow scripts (pytest)
# make test   -- check + lint + unit
# Regression test (needs snakemake; same as the CI dry-run job):
#   bash tests/run_test.sh               # synthetic-data dry-run, validates DAG integrity
#   bash tests/run_test.sh --consensus   # dry-run with the optional v0.2 stages enabled
#   bash tests/run_test.sh --real-run    # end-to-end run + output assertions (server validation)

PYTHON ?= python3

SHELL_SCRIPTS := run.sh $(wildcard workflow/scripts/*.sh)

.PHONY: check lint unit test

check:
	@for f in $(SHELL_SCRIPTS); do bash -n $$f && echo "  PASS  bash -n $$f" || exit 1; done

lint:
	bash tests/lint.sh

unit:
	$(PYTHON) -m pytest tests -q

test: check lint unit
