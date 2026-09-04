# Changelog

所有对项目的显著变更将记录在本文件。格式参考 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/)。

## [Unreleased]

### 待办

- 服务器最小样本端到端实跑（CI 首跑 + conda 环境求解 + MACS2 无对照 control_lambda 确认；验证入口 `bash tests/run_test.sh --real-run`）
- DiffBind 差异分析补完（需 contrast/设计公式决策；v0.4.0 起空壳脚本归档于 `legacy/diffbind/`）
- bowtie2 `.bt2l` 大基因组索引支持
- snakemake 8.x 的 executor-plugin 风格 profile（`snakemake-executor-plugin-cluster-generic`）

## [0.4.0] - 2026-09-04

向姊妹项目 rna-seq（v0.8.0）工程体系全面对齐：目录布局、环境管理、启动体验、资源模型、测试与文档。
设计文档见 `docs/superpowers/specs/2026-09-04-rna-seq-alignment-design.md`，实施计划见 `docs/superpowers/plans/2026-09-04-rna-seq-alignment.md`。

### Added（对齐 rna-seq v0.8.0）

- **统一环境体系三件套**：`workflow/environment.yaml`（钉版一体化主环境，合并原 11 个 per-rule envs 的全部版本约束）、`config/software.yaml`（environment.type=system/conda_prefix/conda_name + R runtime + 工具路径覆盖）、`workflow/scripts/runtime_config.py`（`export` 注入 `CHIP_*` 环境变量 / `check` preflight 双子命令）
- **software_versions 规则 + collect_versions.py**：运行期把实际工具版本、git commit、运行时模式记录到 `5.QC/software_versions.yaml`（复现审计与方法节依据）
- **run.sh 完整运维 CLI**（535 行，替代 main_run.sh）：`--profile auto|default|pbs|sge|slurm`（auto 探测：sbatch→slurm；qsub 按 SGE_ROOT 消歧 PBS/SGE）、per-rule 资源占位符集群提交串、`--memory/--runtime/--queue/--partition` 覆盖、`--check-software/--check-r` preflight、`--validate-only`、`--unlock`、`--retries`、`--log` tee + trap 计时、`--` 透传、`-r` 原始数据重命名、config.local.yaml 自动叠加
- **per-rule 资源模型**：22 条规则全部声明 `mem_mb/runtime_min/runtime_sec`（4 条规则补 threads:1）；`config.yaml` 新增 `resources:` 覆盖段（按规则名覆盖）；`res()` helper 支持项目级调参
- **四套 profile**：`workflow/profile/{default,pbs,sge,slurm}/config.yaml`（资源占位符统一；pbs walltime 用秒避免格式歧义）
- **MultiQC 定制**：`workflow/multiqc_config.yaml`（标题/流程标识/说明），multiqc 规则 `-c` 接线
- **测试与 CI**：`tests/lint.sh`（六段静态检查，缺工具自动跳过）；`tests/make_testdata.py`（确定性合成数据生成器：2×100kb 参考基因组 + chr1 三峰区富集的 chip/atac PE reads + 样本表 + 测试 config，seed 固定逐字节可复现）；`tests/run_test.sh`（合成数据 dry-run 回归默认 + `--real-run` 服务器实跑开关 + 产物断言）；CI 重写为双 job（lint + dry-run 回归，后者无需 conda 环境）
- **文档四件**：`docs/使用说明.md`（8 章操作手册）、`CONTRIBUTING.md`、`README.md` 全面重写（保留 mermaid/QC 阈值表/结果速查亮点）、`example/` 真实项目模板（2 样本表 + 项目 config + 一条命令指引）

### Changed（破坏性）

- **目录布局迁移至 Snakemake 标准**：`workflow.smk` → `workflow/Snakefile`（纯编排）；约 200 行共享定义（样本表解析/config 校验/查询函数/目标汇总）拆至 `workflow/rules/common.smk`；rules/ 与 R 脚本迁入 `workflow/{rules,scripts}/`；样本表模板 → `config/samples.csv`；profile → `workflow/profile/`（全部 git mv 保留历史）
- **启动入口更替**：`main_run.sh` 移除，统一 `bash run.sh`（位置参数=工作目录，`-P/-w` 均可）；旧 `-b/-e/-E` conda 部署选项随环境路线切换移除
- **样本表解析失败行为**：`grouplist` 三级解析（绝对 > 工作目录 > 仓库）都找不到文件时直接报错（不再静默回落仓库根示例文件）；默认值改指 `config/samples.csv` 模板
- 5 个未接入 DAG 的独立 QC 空壳脚本（DiffBind/ChIPQC/DROMPAplus 等）归档 `legacy/diffbind/`（映射表见 `legacy/README.md`）
- 单元测试 45 → **55 项**（删 envs 完整性 1 项；增资源声明 5 项 + 合成数据生成器 6 项）

### Removed

- **per-rule conda 体系**（!）：`envs/` 11 个环境文件与规则内 19 处 `conda:` 指令全部移除，改统一环境路线；服务器需一次性 `mamba env create -f workflow/environment.yaml` 或经 `config/software.yaml` 复用已有环境（三条落地路径见使用说明 §1）

### Fixed

- 工具名映射保真：preflight/版本记录的 deeptools → `bamCoverage`（代表二进制）、spp → `run_spp.R`（与 spp_qc.smk 实际调用一致，消除保证性误报）
- meta.smk 的 python 解释器经 `CHIP_PYTHON` 注入（对齐 rna-seq）；PBS `.o` 日志回收加 profile 守卫（不再误吞 slurm 输出）

## [0.3.0] - 2026-09-03

### Added（运维审查 P3 批次）

- **FRiP / NSC-RSC 注入 MultiQC 报告**：frip_summary 与新增 spp_summary 规则产出 `_mqc.tsv` 自定义表（custom content），multiqc 汇总规则自动纳入——QC 指标（fastqc + bowtie2 + picard + FRiP + NSC/RSC）集中单报告；qc 开关闭合与条件 include 联动
- **分析窗口参数统一**：新 `region_flank`（默认 3000）一键控制 ChIPseeker flank/TSS 窗口与 deeptools computeMatrix 上下游长度（原先三处硬编码）；`annoPeak_batch.R` 新增第四参数
- **`profiles/pbs/`**：snakemake 7.x PBS profile（集群参数固化入库：`{rule}` 任务名、`{threads}` 核数、latency-wait/rerun-incomplete 默认值）；8.x 用户继续用 main_run.sh（版本自动适配），README 注明迁移方向
- 测试新增 4 项（region_flank 校验 ×2 + **mqc shell 规则体实测** ×2——按 snakemake 同款解析链 `ast.literal_eval` 解码 + format 渲染 + bash 实际执行并断言产物格式），共 **45 项**全部通过

### Changed

- SPP 规则移除 `-savp`（pdf 副产物文件名依输入 BAM 派生且落 cwd，不可声明管理；全部指标已含于 `-out` 文本）
- `envs/bigwig.yaml`：ucsc-bedclip/bedgraphtobigwig 锁 bioconda 构建号 482（消除无语义版本的漂移风险）
- `envs/trim-galore.yaml`：移除冗余的显式 `cutadapt=4.4` 钉（由 trim-galore 依赖自行拉动，双钉增加求解冲突面）

## [0.2.2] - 2026-09-03

### Added（运维审查 P2 批次）

- **`main_run.sh` 集群健壮性**：默认启用 `--rerun-incomplete` 与 `--latency-wait`（`-t` 可调，默认 90s），覆盖 PBS 断点重跑与共享文件系统输出可见性延迟两类常见假失败
- **共享 conda 环境目录**：新 `-e DIR`（`--conda-prefix`），多项目复用同一套环境；新 `-E` 预建模式（`--conda-create-envs-only`），PBS 计算节点无外网时先在登录节点建环境
- **trim_galore 参数配置化**：新 `trim.quality/stringency/error_rate/extra` 四键（原硬编码 `-q 25 --stringency 3 -e 0.1`），template 同步
- **config 集中校验**（`workflow.smk` `validate_config`）：必需键/子键/类型/取值范围一次汇总报出（含 threads 非整数的友好报错）；参考文件缺失仅警告不中断（保证 --lint/dry-run 在无参考文件的机器可解析）
- 测试新增 10 项（validate_config 真实源码提取执行），共 **45 项**全部通过

### Changed

- **callpeak 三条规则 threads 降为 1**：MACS2 为单线程程序，原先按 config threads=12 申请集群资源造成超订/浪费；配合 v0.2.1 的 `ncpus={threads}` 后 PBS 申请与实际占用精确一致

## [0.2.1] - 2026-09-03

### Fixed（运维审查 2 项 P1）

- **PBS 资源参数与流程线程数脱钩**：`main_run.sh`/README 的集群提交示例改为 `-l ncpus={threads}`——snakemake 按每个任务实际线程数格式化 cluster 串，与 config `threads` 自动对齐（原示例固定 `ncpus=6` 与默认 `threads: 12` 矛盾，导致超订）
- **snakemake 版本矩阵不明确**：`main_run.sh` 启动时探测 snakemake 主版本——≥8 自动使用 `--software-deployment-method conda`（`--use-conda` 在 8.x 已弃用），7.x 继续用 `--use-conda`；未安装/版本不可解析时给出明确报错。README 新增版本支持矩阵（7.32.4 参考版本 / 8.x 自动适配 / <7 与 ≥9 未验证），方式二手动命令同步加注

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
