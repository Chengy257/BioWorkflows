# Changelog

所有对项目的显著变更将记录在本文件。格式参考 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/)。

## [Unreleased]

### 待办（见 docs/IMPROVEMENT_PLAN.md 遗留清单）

- 服务器最小样本端到端实跑（CI 首跑 + conda 环境求解 + MACS2 无对照 control_lambda 确认）
- DiffBind 差异分析补完（需 contrast/设计公式决策）
- bowtie2 `.bt2l` 大基因组索引支持

## [0.2.0] - 2026-09-03

> v0.1.0 原始实现整体归档至 `legacy/`。重构经独立代码审查（fix-first），审查发现的 3 项 P1 与全部 P2/P3 已在同版本内修复。

### Added

- 统一入口 `workflow.smk`：新样本表 schema（sample_id/role/group/seqtype/layout/peak_type）、逐行校验（行号报错、名字字符集校验）、按 assay 自动路由；旧 4 入口归档
- 全新规则集 `rules/`：upstream（trim/fastqc/multiqc/bowtie2，fastqc 按样本并行）、dedup（picard，按 assay 开关，CUT&Tag 默认不去重）、callpeak（narrow/broad/atac 三规则按组并行 + bigwig）、annotation（ChIPseeker 批量注释）、frip（FRiP + 汇总表）、qc_deeptools（相关性/PCA/指纹/片段长/基因区信号全套）、spp_qc（可选 NSC/RSC）
- per-rule conda 环境 `envs/`（11 个，仅 conda-forge+bioconda），替代非法的 `conda: "chip"` 写法
- 配置体系：完整 config schema（genome_size/min_mapq/dedup 按 assay/peak 阈值/qc 开关）+ `config.template.yaml` + `config.local.yaml` 叠加机制
- 阈值常规默认：narrow q=0.05、broad_cutoff=0.05、min_mapq=30（ENCODE 常规值），全部可配置
- ATAC/FAIRE 峰调用双模式：`bampe`（默认，ENCODE ATAC v2 做法）| `shifted`（经典 Tn5 偏移配方）
- `main_run.sh` 重写：getopts 参数化（含 dry-run 预检、PBS 可选、config.local 自动叠加、集群/本机 -j/--cores 正确拆分）
- 测试与 CI：`tests/run_tests.py`（31 项零依赖单元测试，直接抽取 workflow.smk 真实源码执行）、Makefile（check/lint/dryrun）、GitHub Actions（测试 + shellcheck + snakemake --lint）
- MIT LICENSE；`docs/REVIEW.md` 修复状态表；`legacy/README.md` 归档映射

### Fixed（对 v0.1.0 的全部 P0/P1，详见 docs/REVIEW.md §七）

- 峰调用链接入 DAG（原 include 断裂、rule all 目标全被注释）
- 去重规则 shell 命令误写进 conda 块、`${id}` 非法通配符、threads 字符串类型
- 峰调用循环变量未用导致每组重复调用两次、空规则输出冲突
- 样本表 schema 与消费代码矛盾（三种互不兼容的定义并存）
- `get_samples()` 用 seqtype 列判断样本去重的逻辑错误
- 注释规则引用不存在的脚本路径
- 全部 16 处机器特定绝对路径（/home/chengyu、/opt、/share）
- `run_ChIPQC.R`：args 未定义即使用、`basename()` 缺参、教程残留名、私人路径
- bigwig 染色体字典序与 chrom.sizes 顺序不匹配（改 `bedtools sort -g`）
- r-chipseeker conda 环境 R 4.2 与 Bioc 3.18 包版本冲突（r-base=4.3）
- ATAC shift/extsize 在 BAMPE 模式下被静默忽略（改双模式开关）
- MACS2 参数 `-q 0.5` 等宽松阈值（改常规默认并配置化）
- 根目录与 scripts/ 重复的 call_peak.sh

### Changed

- README 全面重写（新用法/新 schema/QC 阈值/已知限制）；换行符策略 `.gitattributes` 强制 LF

## [0.1.0] - 2024-03

### Added（原始快照，commit 03f5c0c，已归档至 legacy/）

- ChIP-seq / CUT&Tag / ATAC-seq / FAIRE-seq 四入口 Snakemake 流程
- 上游规则：trim_galore、FastQC、MultiQC、bowtie2 比对
- 峰调用脚本：SPP 片段长估计 + MACS2（narrow/broad/ATAC 模式）、HMMRATAC、FSeq2
- MACS2 bdgcmp → bigWig 转换脚本
- ChIPseeker 批量/单样本峰注释 R 脚本
- deeptools / ChIPQC / DROMPAplus QC 脚本（部分为片段或空壳）
- 旧版整环境导出 `chip_environment.yaml` 与 PBS 启动脚本 `main_run.sh`
