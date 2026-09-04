# chip_cuttag_atac_faire 开发检查入口
# make check  —— 单元测试 + shell 语法检查（无需 snakemake）
# make lint   —— 统一静态检查套件 tests/lint.sh（bash -n/shellcheck/py_compile/
#               R parse/yaml/snakemake --lint，缺可选工具自动跳过，CI 全量执行）
# make dryrun —— dry-run DAG 预检（需 snakemake + 已备好 1.rawdata）
# make test   —— CI 同款全量检查

PYTHON ?= python3
SNAKEMAKE ?= snakemake

SHELL_SCRIPTS := run.sh $(wildcard workflow/scripts/*.sh) $(wildcard scripts/*.sh)

.PHONY: check lint dryrun test

check:
	$(PYTHON) tests/run_tests.py
	@for f in $(SHELL_SCRIPTS); do bash -n $$f && echo "  PASS  bash -n $$f" || exit 1; done

lint:
	bash tests/lint.sh

dryrun:
	$(SNAKEMAKE) -n -s workflow/Snakefile --configfile config/config.yaml \
		--profile workflow/profile/default

test: check lint
