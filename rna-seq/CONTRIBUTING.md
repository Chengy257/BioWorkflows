# 贡献指南（CONTRIBUTING）

## 1. 开发流程

1. **分支**：`main` 保持可运行状态；功能开发与修复在特性分支进行（`feat/<主题>`、`fix/<主题>`），自测通过后合并。
2. **提交信息**：`<类型>: <中文概述>`，类型与仓库现有历史保持一致：
   - `feat` 新功能 / `fix` 缺陷修复 / `refactor` 结构调整（不改行为）/ `docs` 文档 / `chore` 杂项
   - 一个提交只做一件事；行为变更与纯重构分开提交（便于回溯定位）。
3. **行尾**：脚本统一 LF（`.gitattributes` 已约束）；运行产物一律不入库（`.gitignore`）。

## 2. CHANGELOG 要求

- **任何影响用户行为的变更**（新增/修改/移除 config 键、规则输入输出路径、脚本参数、运行方式）必须记入 `CHANGELOG.md` 未发布小节，遵循 Keep a Changelog 格式。
- 纯文档或注释级改动可酌情合并记录。

## 3. 文档同步 checklist（改规则必查）

修改 `workflow/rules/` 或 `workflow/Snakefile` 前，逐项核对：

- [ ] 新增/修改了 **config 键**？→ 同步 `config/config.yaml` 注释模板 + `docs/使用说明.md` §4.3 配置表
- [ ] 新增/修改了 **输出文件路径或目录**？→ 同步 `docs/使用说明.md` §6 结果解读 + `tests/check_outputs.py` 断言
- [ ] 修改了 **规则入参或脚本调用参数**？→ 同步对应 `workflow/scripts/` 脚本的 Usage 注释 + 使用说明 FAQ
- [ ] 新增了 **外部工具依赖**（conda 之外）？→ 更新 `docs/使用说明.md` §2 系统要求与 §7 已知限制
- [ ] 修改了 **conda 环境**（`workflow/envs/*.yaml`）？→ 确认版本 pin 合理，并在 CHANGELOG 记录
- [ ] 是否破坏 **旧版本产物兼容**？→ 在 CHANGELOG 用显著条目说明迁移方式

## 4. 测试要求

提交影响流程逻辑的变更前，在 Linux 环境运行：

```bash
bash tests/lint.sh                          # 静态检查（bash/shellcheck/py/R/snakemake --lint）
bash tests/run_test.sh --pipeline deg       # 端到端回归（生成数据 → dry-run → 运行 → 断言）
```

- `deg` 为默认验收管线（覆盖 align+quant+deg 全部规则）；改动比对/组装相关规则时另跑 `--pipeline as`。
- 新增输出文件时，请在 `tests/check_outputs.py` 中补对应断言。
- CI（GitHub Actions）会在 push/PR 时自动执行 lint 与 `--reads 20000` 的快速回归；本地无 CI 环境时以上两条命令为准。

## 5. 代码风格

- Snakefile/规则：与现有规则保持相同的缩进与 directive 顺序（input → output → params → log → conda → shell）；shell 串中 `{}` 转义遵循 `{{}}`。
- Shell 脚本：`set -u` 起步、路径加引号、`SCRIPT_DIR` 自解析；优先可移植写法（GNU coreutils 为主）。
- Python：标准库优先（测试生成器/校验器不得引入第三方依赖）；R：保持现有 `getopt` 风格。
