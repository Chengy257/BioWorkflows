# chip_cuttag_atac_faire 开发检查入口
# make check  —— 单元测试 + shell 语法检查（无需 snakemake）
# make lint   —— 统一静态检查套件 tests/lint.sh（bash -n/shellcheck/py_compile/
#               R parse/yaml/snakemake --lint，缺可选工具自动跳过，CI 全量执行）
# make test   —— CI 同款全量检查（= check + lint）
# 回归测试（需 snakemake，即 CI dry-run job 同款）:
#   bash tests/run_test.sh             # 合成数据 dry-run，验证 DAG 完整性
#   bash tests/run_test.sh --real-run  # 端到端实跑 + 产物断言（服务器验证用）

PYTHON ?= python3

SHELL_SCRIPTS := run.sh $(wildcard workflow/scripts/*.sh) $(wildcard scripts/*.sh)

.PHONY: check lint test

check:
	$(PYTHON) tests/run_tests.py
	@for f in $(SHELL_SCRIPTS); do bash -n $$f && echo "  PASS  bash -n $$f" || exit 1; done

lint:
	bash tests/lint.sh

test: check lint
