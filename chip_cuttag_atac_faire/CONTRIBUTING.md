# 贡献指南（CONTRIBUTING）

## 1. 开发流程

1. **分支**：`master`/`main` 保持可运行状态；功能开发与修复在特性分支进行（`feat/<主题>`、`fix/<主题>`），自测通过后合并。
2. **提交信息**：`<类型>: <中文概述>`（可带 scope，如 `feat(launcher): ...`；破坏性变更加 `!`，如 `feat(launcher)!: ...`），类型与仓库现有历史保持一致：
   - `feat` 新功能 / `fix` 缺陷修复 / `refactor` 结构调整（不改行为）/ `docs` 文档 / `test` 测试 / `ci` CI / `chore` 杂项
   - 一个提交只做一件事；行为变更与纯重构分开提交（便于回溯定位）。
3. **行尾与产物**：脚本统一 LF（`.gitattributes` 已约束）；运行产物与测试数据一律不入库（`.gitignore`；`tests/run_test.sh` 默认自清理，`--keep` 产物勿提交）。

## 2. CHANGELOG 要求

- **任何影响用户行为的变更**（新增/修改/移除 config 键或 `software.yaml` 键、`run.sh` 选项、规则输入输出路径、脚本参数、运行方式、依赖钉版）必须记入 `CHANGELOG.md` 未发布小节，遵循 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/) 格式（Added/Changed/Fixed/Removed/Deprecated 分节）。
- 破坏性变更（`!`）需在条目中写明迁移方式；纯文档或注释级改动可酌情合并记录。
- 发版时将未发布小节落为版本号 + 日期，并打语义化版本 tag。

## 3. 文档同步 checklist（改规则必查）

修改 `workflow/`、`run.sh` 或 `config/` 前，逐项核对：

- [ ] 新增/修改了 **config 键**？→ 同步 `config/config.yaml` 注释、`config/config.template.yaml`、`docs/使用说明.md` §4 全键表（涉集群资源再加 §4.3 默认值表）
- [ ] 新增/修改了 **规则 resources 默认值**？→ 同步 `docs/使用说明.md` §4.3 默认值表
- [ ] 新增/修改了 **输出文件路径或目录**？→ 同步 `docs/使用说明.md` §7 目录树与结果速查表、`tests/run_test.sh` 的实跑断言
- [ ] 新增/修改了 **样本表列或校验规则**？→ 同步 `config/samples.csv` 模板、`docs/使用说明.md` §3、`tests/run_tests.py` 对应断言
- [ ] 修改了 **`run.sh` 选项或行为**？→ 同步 `run.sh` 内 usage 文本、`docs/使用说明.md` §2/§5/§6
- [ ] 新增/修改了 **软件或 R 依赖**？→ 同步 `workflow/environment.yaml`（钉版）与 `config/software.yaml`；不要在 rule 中重新引入独立 `conda:` 指令（v0.4.0 起为 external-runtime 体系）
- [ ] 是否破坏 **旧版本产物/用法兼容**？→ 在 CHANGELOG 用显著条目说明迁移方式，必要时更新 README 版本矩阵

## 4. 测试要求

提交影响流程逻辑的变更前，按环境可用度运行：

```bash
make check                            # 55 项单元测试 + bash -n 语法检查（无需 snakemake）
make lint                             # 静态检查套件 tests/lint.sh（bash/shellcheck/py_compile/R parse/yaml/snakemake --lint，缺可选工具自动跳过）
make test                             # CI 同款全量检查（= check + lint）
bash tests/run_test.sh                # 合成数据 dry-run 回归（需 snakemake + python3/PyYAML；CI 传 --reads 2000）
bash tests/run_test.sh --real-run     # 端到端实跑 + 产物断言（需完整分析环境，服务器验证用）
```

- 本地最低门槛为 `make check`；改了 Snakefile/rules/config 的 PR 在有 snakemake 的环境必须补 `bash tests/run_test.sh`。
- `--real-run` 需含 bowtie2/fastqc/trim_galore/macs2/deeptools/R（ChIPseeker）的完整环境（即 `workflow/environment.yaml` 同款），在服务器上验证部署时使用。
- 新增/修改输出文件时，在 `tests/run_test.sh` 实跑断言中补对应条目；新增配置校验逻辑时，在 `tests/run_tests.py` 补单元测试并保持 **55 项基线只增不减**。
- CI（GitHub Actions，`.github/workflows/ci.yaml`）两个 job：lint（make check + tests/lint.sh）与合成数据 dry-run 回归（`bash tests/run_test.sh --reads 2000`）；本地无 CI 环境时以上命令为准。

## 5. 代码风格

- **命名**：规则、config 键、脚本、函数统一 `snake_case`；样本/分组名沿用样本表字符集约束（字母数字与 `. _ -`）。
- **集中校验**：config 与样本表校验只在解析期集中执行（`workflow/rules/common.smk` 的 `validate_config` / `load_sample_table`），错误信息**带行号、一次汇总报出**；不要把校验散落到各 rule 的 shell 里。
- **注释语言**：config 模板、规则文件与校验错误信息以中文为主（面向最终用户）；`run.sh` 的 usage/help 文本保持英文（与 rna-seq 体系对齐），行内日志可中文。同一文件内风格保持一致。
- **Snakemake**：与现有规则保持相同缩进与 directive 顺序（input → output → params/log → threads/resources → shell）；集群资源经 `res()` 读取并支持 config `resources:` 段覆盖，勿在 shell 串里硬编码内存/时长。
- **Shell**：`set -Eeuo pipefail` 起步（维护脚本至少 `set -euo pipefail`）、路径加引号、`REPO_DIR` 自解析；优先可移植写法（GNU coreutils 为主）。
- **Python**：标准库优先（测试生成器/校验器不得引入第三方依赖；`tests/make_testdata.py` 需固定种子、逐字节可复现）。
- **格式**：缩进 4 空格（smk/py/sh）、yaml/md 2 空格、Makefile 用 tab（`.editorconfig` 已约定）。
