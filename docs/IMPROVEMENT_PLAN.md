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

## Phase 1 — P0 修复：让流程真正跑通（✅ 已完成，2026-09-03，commit `e29c428`）

目标：在不改变分析逻辑的前提下修复所有阻断性问题，4 种 assay 全链路可达。

1. ✅ **统一入口结构**：4 入口合并为 `workflow.smk`，按样本表 `seqtype` 自动路由；旧入口归档 `legacy/`；全部 include/脚本经 `workflow.basedir` 相对定位，无任何 `/home/chengyu` 路径。
2. ✅ **接入峰调用链**：rule all 覆盖 QC → 比对 → 去重 → 峰 → bigWig → 注释 → FRiP/deeptools。
3. ✅ **修复去重规则**：重写为 `rules/dedup.smk`（shell 块归位、`{sample}` 通配符、threads int、统一路径）。
4. ✅ **修复峰调用规则**：空规则删除、双重调用删除、输出契约与 MACS2 `-n {group}` 命名一致、注释规则调用 `scripts/annoPeak_batch.R`（workflow.basedir 定位）、按组并行。
5. ✅ **统一样本表 schema**：`sample_id,role,group,seqtype,layout,peak_type`，入口逐行校验（错误含行号），`get_samples()` 重写为 `load_sample_table()`，`sample_info.example.csv` 文档化。
6. ✅ **conda 环境拆分**：`envs/` 下 11 个 per-rule 环境；原整环境导出归档 `legacy/chip_environment.yaml`。
7. ✅ **参数类型与可配置**：threads int；MAPQ/genome_size/qvalue/broad_cutoff/keepdup 全部入 config。

**验收标准**
- [x] 样本表解析/校验/路由逻辑：31 项单元测试全部通过（`python tests/run_tests.py`，测试直接抽取 workflow.smk 真实源码执行）
- [x] `grep -rn "/home/chengyu\|/opt/\|/share/" rules/ workflow.smk` 无结果（仅 config 中参考路径默认值，属用户配置）
- [ ] 服务器小样本端到端实测（**遗留**：本机为 Windows 无法运行 snakemake；CI 已配置 `--lint` + `--list-rules`，首次推送后运行）
- [x] 独立代码审查（fix-first 裁决）→ 3 项 P1 全部修复并复验（r-chipseeker R 版本、ATAC mode 死旋钮、bigwig 排序顺序）

---

## Phase 2 — P1 可移植性与参数化（✅ 已完成，2026-09-03）

1. ✅ **config 重构**：完整 schema + `config/config.template.yaml` + `config.local.yaml` 自动叠加（main_run.sh -l 或自动检测）。
2. ✅ **main_run.sh 重写**：getopts 参数化、logs 自动创建、PBS 可选、dry-run 预检（-n）、集群/本机模式的 -j/--cores 正确拆分、不再删除 .snakemake。
3. ✅ **CUT&Tag 去重策略落地**：`dedup.cuttag: false`（Kaya-Okur et al. 2019），按 assay 配置化，`sample_bam()` 全流程一致路由。
4. ✅ **QC 链路打通**：bigWig 规则声明 `*_FE.bw` 输出；deeptools 全套规则化（`rules/qc_deeptools.smk`）；FRiP 规则 + 汇总表；multiqc 聚合 fastqc + bowtie2 日志 + picard 指标；SPP NSC/RSC 可选模块；`run_ChIPQC.R` 修复（commandArgs bug/私人路径/教程残留）。
5. ✅ **脚本去硬编码**：chromsize/ucsc-tools/docker 镜像参数化；R 脚本去 `.libPaths`；DROMPAplus 参数化。
6. ✅ **规则粒度优化**：fastqc 按 sample 拆分；callpeak/bigwig/frip 按 group 并行。

**验收标准**
- [x] 代码目录无机器特定路径（换机器只改 config）
- [x] 日志结构化留存 `logs/`；multiqc 含比对与去重指标（FRiP 为独立汇总表）
- [x] CUT&Tag 与 ChIP 的去重/峰调用参数差异由 config 显式表达
- [ ] 干净 Linux 机器全流程复跑（**遗留**，与 Phase 1 端到端实测合并）

---

## Phase 3 — P2/P3 工程化与长期维护（✅ 主体完成，2026-09-03）

1. ✅ **清理**：重复 `call_peak.sh` 删除；legacy 归档 + `legacy/README.md` 映射表。
2. ⬜ **测试数据**：`tests/run_tests.py`（41 项零依赖单元测试）+ CI 已就绪；端到端小样本实测待服务器（遗留）。
3. ✅ **CI**：`.github/workflows/ci.yaml`（单元测试 + shellcheck + snakemake --lint + --list-rules）。
4. ✅ **文档**：README 重写（新 schema/新入口/QC 阈值/已知限制）；各脚本头部 usage。
5. ✅ **版本管理**：`v0.2.0` tag；CHANGELOG 维护。
6. ✅ **LICENSE**：MIT（作者已确认）。
7. ⬜ **可选进阶**：对标 nf-core 输出结构（远期）；DiffBind 差异分析补完（**需用户决策 contrast/设计公式，明确排除在本轮外**）。

---

## 遗留事项清单（需用户/服务器条件）

| 事项 | 阻塞原因 |
|---|---|
| 服务器最小样本端到端实跑（myc/IgG 两组） | 本机 Windows 无法运行 snakemake；推送 GitHub 后 CI 亦可先行验证 lint |
| CI 首次运行确认 | 需推送到 GitHub |
| conda 环境真实求解（重点 r-chipseeker / phantompeakqualtools） | 需 Linux + conda |
| MACS2 无对照时 `control_lambda.bdg` 产出确认 | 外部证据（MACS issues #275/#291）表明会产出，需实测 |
| DiffBind 差异分析实现 | 需 contrast 与设计公式决策 |
| bowtie2 `.bt2l` 大基因组支持 | 低频需求，README 已声明 workaround |
