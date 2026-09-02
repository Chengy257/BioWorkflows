# 优化、完善与规范化计划（IMPROVEMENT PLAN）

> 依据 `docs/REVIEW.md` 的审查结论制定。原则：先让流程"跑通且正确"（P0），再"换台机器也能跑"（P1），最后"像标准开源项目"（P2/P3）。
> 每个阶段有明确验收标准，完成一项勾一项。

---

## Phase 0 — 规范化基建（✅ 已完成，2026-09-03）

| 事项 | 状态 |
|---|---|
| 建立 git 仓库，原始代码作为快照入库（保留历史可回溯） | ✅ commit `03f5c0c` |
| `README.md` 使用说明文档 | ✅ |
| `docs/REVIEW.md` 全面审查报告 | ✅ |
| `docs/IMPROVEMENT_PLAN.md` 本计划 | ✅ |
| `.gitignore`（排除 0.index~6.motif 等分析产出与大数据文件） | ✅ |
| `.gitattributes`（强制 LF，防止 Windows CRLF 破坏 Linux 脚本） | ✅ |
| `.editorconfig` / `CHANGELOG.md` | ✅ |

---

## Phase 1 — P0 修复：让流程真正跑通（预计 1~2 天）

目标：在不改变分析逻辑的前提下修复所有阻断性问题，4 种 assay 全链路可达。

1. **统一入口结构**
   - 4 个入口 smk 合并为单一入口 `workflow.smk`（保留 4 个旧文件暂不删，标注 deprecated），按样本表 `seqtype` 列自动路由 assay 分支。
   - 所有 `include` / `configfile` / 脚本引用改为相对 workflow 源码目录的路径（`workflow.basedir`），消除 `/home/chengyu/...`。
2. **接入峰调用链**：`rule all` 恢复 peak → 注释 → QC 目标；include 对应 rules。
3. **修复 `rules/rmDup.smk`**：conda 块还原为 shell 块；`${id}` → `{sample}`；threads `int(config["threads"])`；统一 `3.align/bowtie2/` 路径。
4. **修复 `callpeak_*.smk`**：
   - 删除空的 `callpeak_ATAC` 规则；
   - 修 `for j in ...` 双重调用 bug；
   - `${id}` → `{sample}`；
   - `peakAnno` 改调 `scripts/annoPeak_batch.R`（经 `workflow.basedir` 定位）；
   - bdgcmp/peakAnno 改为按 group 去重后调用，或拆成按组并行规则。
5. **统一样本表 schema**（关键设计决策）：
   建议列：`sample_id, role(treat|control), group, seqtype(chip|cuttag|atac|faire), layout(PE|SE), peak_type(narrow|broad)`；
   同步修改 `get_samples()`（按 `sample_id` 去重）与所有 shell 循环；提供 `sample_info.example.csv` 并在 README 文档化。
6. **conda 环境拆分**：从 `chip_environment.yaml` 拆出 rule 级环境到 `envs/`：`trim_fastqc.yaml`、`bowtie2_samtools.yaml`、`picard.yaml`、`macs2.yaml`、`deeptools.yaml`、`multiqc.yaml`、`r-chipseeker.yaml`；`conda: "chip"` 全部改为对应文件路径。原文件移入 `legacy/` 或加 deprecation 头注释。
7. **threads/参数类型**：config 数值统一 `int()`；MAPQ、基因组大小、q 值等入 config。

**验收标准**
- [ ] `snakemake -n --use-conda`（dry-run）对 4 种 assay 均能生成完整 DAG，无语法/路径报错；
- [ ] `grep -rn "/home/chengyu\|/opt/\|/share/" rules/ workflow.smk` 无结果；
- [ ] 用一套 2 处理 + 2 对照的小样本实测跑通 QC → 比对 → 峰调用 → 注释全链路。

---

## Phase 2 — P1 可移植性与参数化（预计 2~3 天）

目标：换机器只需改 config；每步产物命名与 QC 链路对齐。

1. **config 重构**：`config/config.yaml` 扩展为完整 schema（参考路径四件套 + genome_size + peak 参数 + 去重策略 + 工具路径），提供 `config/config.template.yaml`；本地私有配置用 `config.local.yaml` 覆盖（不入库）。
2. **`main_run.sh` 重写**：getopts 参数化（`-w workflow -d workdir -c config -j jobs`）、自动创建 logs、PBS/SLURM profile 可选、启动前先 `--dry-run` 预检。
3. **CUT&Tag 去重策略落地**：按文献实现"不去重 + 保留 `"--keep-dup all`"（对照 Kaya-Okur et al. 2019），落实 `call_peak.sh:143` 的 TODO；去重开关进 config 按 assay 配置。
4. **QC 链路打通**：
   - bigWig 生成规则声明 `*_FE_bdgcmp.bw` 输出，deeptools 规则改用该命名；
   - `run_deeptools_QC.sh` 片段转成真正的 Snakemake 规则；
   - 修复 `run_ChIPQC.R`（`commandArgs(T)`、路径参数化、reportName 去教程化）；
   - FRiP、peak 数、NSC/RSC 汇总表接入 MultiQC 自定义内容。
5. **脚本去硬编码**：`call_peak.sh`/`bdgcmp_macs2.sh` 的 chromsize、ucsc-tools、HMMRATAC、run_spp 路径全部改为 config 传参；R 脚本移除 `.libPaths` 强制覆盖；DROMPAplus docker 镜像名参数化。
6. **规则粒度优化**：fastQC 按 sample 拆分；callpeak/bdgcmp 按 group 拆分为可并行规则。

**验收标准**
- [ ] 换一台干净 Linux 机器，仅需安装 conda + snakemake、修改 config 即可复跑（grep 检查无机器特定路径）；
- [ ] 全流程日志在 `logs/` 结构化留存，MultiQC 报告含比对率/去重率/FRiP/NSC/RSC 汇总；
- [ ] CUT&Tag 与 ChIP 的去重/峰调用参数差异由 config 显式表达。

---

## Phase 3 — P2/P3 工程化与长期维护（持续）

1. **清理**：删除根目录重复的 `call_peak.sh`（`scripts/` 为唯一真身）；删除或补完 `diffpeak_DiffBind.sh` / `run_DiffBind.R`（补完方案：DiffBind samplesheet 由样本表生成 + contrast 列设计）。
2. **测试**：`tests/` 下放入小型模拟 fastq（可用 public rice 数据截取或 `seqkit` 模拟），`make test` 一键端到端冒烟测试；`snakemake --lint` 通过。
3. **CI**：GitHub Actions：`snakemake --lint` + `shellcheck` + `yamllint`（每次 push）。
4. **文档**：每个脚本头部 usage 说明；README 增加 FAQ、参数依据文献（Kaya-Okur 2019 CUT&Tag、MACS2 手册、ENCODE ChIP QC 阈值 NSC≥1.05/RSC≥0.8 等）。
5. **版本管理**：从 `v0.2.0` 起语义化 tag，维护 `CHANGELOG.md`；主分支保护 + feature 分支工作流。
6. **LICENSE**：作者定夺后添加（个人学术项目建议 MIT 或 Apache-2.0）。
7. **可选进阶**：对标 nf-core/chipseq、nf-core/atacseq、nf-core/cutandrun 的输出结构与 QC 指标，评估是否直接迁移或保持自研对齐其报告规范。

---

## 里程碑建议

| 里程碑 | 内容 | 完成标志 |
|---|---|---|
| M1（Phase 1） | P0 修复 + 小样本全链路实测 | 4 种 assay 产出 peaks + 注释 |
| M2（Phase 2） | 可移植 + QC 汇总 | 干净机器 config 即跑 |
| M3（Phase 3） | 测试/CI/版本化 | CI 绿灯、v0.2.0 tag |
