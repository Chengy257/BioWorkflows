# Changelog

本项目的全部显著变更记录于此。格式参考 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/)。

## [0.8.0] - 2026-09-03

统一软件与 R runtime：工作流默认复用用户已有服务器环境，不再为每个 rule 自动创建独立 Conda 环境。

### Added（新增）
- 新增 `config/software.yaml`：集中配置主 runtime（`system` / Conda prefix / Conda name）、Rscript、R 版本检查、R library path、工具覆盖、CNCI Python2 与 Pfam/NR 数据库。
- 新增 `workflow/environment.yaml`，作为新服务器可选的一体化 Conda 环境模板；仅供用户显式创建，不由 Snakemake 自动部署。
- 新增 `workflow/scripts/runtime_config.py`，统一解析软件配置、构建 PATH / `R_LIBS_USER`、输出 `RNASEQ_TOOL_*` runtime 变量，并按 pipeline 执行 executable / R / R package / database preflight。
- `run.sh` 新增 `--software`、`--check-software`、`--check-r`、`--skip-software-check`；项目目录存在 `software.yaml` 时自动优先使用。
- 新增 `tests/test_runtime_config.sh` 与 `tests/test_r_runtime_propagation.sh`：分别验证 Conda prefix / Rscript / R library 解析，以及 `enrich.sh`、`DEGgroupCompare.sh` 嵌套 R 子任务继承同一 Rscript/R_LIBS。

### Changed（变更）
- 所有 rule 改为复用统一 runtime；STAR、samtools、Trim Galore、MultiQC、StringTie、DIAMOND、bedtools、Rscript/Python 等均支持显式 executable override，并在未覆盖时从主环境/PATH 解析。
- R runtime 成为独立一级配置：允许 Rscript 与主 Conda 环境分离，`r.lib_paths` 通过 `R_LIBS_USER` 统一传递到直接 rule 与 `enrich.sh` / `DEGgroupCompare.sh` 的嵌套 R 调用。
- lncRNA 的 CPC2、CNCI、pfam_scan、Pfam DB、NR DIAMOND DB 从分析 `config.yaml` 迁至 `software.yaml`；CNCI 可独立指定 Python2 解释器。
- `software_versions.yaml` 改为记录实际 runtime：Conda prefix、Rscript/R 版本/R library、解析后的工具路径及大型数据库位置，而不是环境 recipe 列表。
- CI 与回归测试改为先准备一个统一环境再执行 workflow；Snakemake lint 仍检查其它问题，但有意忽略“每 rule 应指定 conda/container”这一与当前架构冲突的提示。

### Removed（移除）
- 删除正式工作流中的 `workflow/envs/*.yaml` 多环境体系、全部 rule `conda:` directive，以及 profile 中的 `use-conda: true`。
- 移除 R 脚本运行期间自动 `install.packages()` 的副作用；所需 R 包必须在运行前已安装于配置的 R library。
- `orgdb_tarball` 不再作为分析参数/运行参数传递；本地 OrgDb tarball 如需保留，可放入 `r.package_sources` 仅用于 preflight 安装提示。

### Migration（迁移）
- 已有服务器环境：复制 `config/software.yaml` 到项目，设置 `environment.conda_prefix`（推荐）或 `conda_name`；若当前 PATH 已完整配置，可直接保持 `environment.type: system`。
- 没有现成环境：可显式执行 `mamba env create -f workflow/environment.yaml`，随后在 `software.yaml` 指向该环境。
- 正式运行前推荐先执行 `run.sh ... --check-software`；R 配置可单独用 `--check-r` 验证。

## [0.7.0] - 2026-09-03

引入 per-rule 资源模型并完成 Snakemake 工程规范化，使 SGE/SLURM 能按任务动态申请 CPU、内存与 walltime，同时让仓库级 lint 真正通过。

### Added（新增）
- 24 个可执行 rule 全部声明 `threads` 与 `resources.mem_mb/runtime_min/runtime_sec`；默认资源覆盖 trim、STAR、featureCounts、DESeq2、StringTie、Pfam、DIAMOND 等不同负载。
- `config/config.yaml` 新增 `resources:` 段，可按 rule 覆盖 `mem_mb`、`runtime_min`，也可显式加入 `threads`。旧项目不含该段时自动使用内置默认值。
- `run.sh` 新增 `--runtime MIN` / `RNASEQ_RUNTIME_MIN`；SGE/SLURM 默认提交模板分别读取 `{resources.mem_mb}` 与 rule walltime，`--memory` / `--runtime` 作为全局强制覆盖。
- 新增 `workflow/envs/utils.yaml` 与 `workflow/rules/meta.smk`，用于 Python 聚合与软件版本记录。

### Changed（变更）
- rule 内实际工具线程统一使用 Snakemake `{threads}`，不再直接读取同一个 `config[threads]`；保留旧顶层 `threads` 与 `lncrna.threads` 的兼容回退。
- helper 函数集中到 `common.smk`，metadata rule 独立到 `meta.smk`；shell 所需脚本/路径/配置通过 `params` 或 input/output 动态推导，减少隐藏全局依赖。
- Conda 环境由 7 个增至 8 个；纯 Python 聚合与 metadata 使用轻量 `utils` 环境。
- `collect_versions.py` 可接收当前运行中的 Snakemake 版本，避免 metadata 环境改变记录值。

### Fixed（修复）
- 四种 pipeline 的 `snakemake --lint` 全部通过；`tests/lint.sh` 也通过 bash、ShellCheck、Python 与 Snakemake lint（当前机器未安装 Rscript，因此 R parse 按设计跳过）。
- 修复 `enrich.sh`、`lncRNA_functions.sh` 与 `tests/lint.sh` 中 `cd` 失败未处理的 ShellCheck SC2164 问题。

## [0.6.0] - 2026-09-03

增强运行入口与集群调度，并修复 Snakemake 7.32.4 的规则解析兼容性问题。

### Added（新增）
- `run.sh` 新增完整英文 CLI：`--help` / `--version` / `--dry-run` / `--validate-only` / `--skip-validation` / `--unlock` / `--log`，并保留原四位置参数调用方式。
- 调度自动选择升级为 SGE（`qsub`）→ SLURM（`sbatch`）→ local；新增 `--profile`、`--queue`、`--partition`、`--memory`、`--sge-mem-resource`、`--scheduler-extra`。
- 新增失败重试与调度限流：`--retries`、`--latency-wait`、`--max-jobs-per-sec`、`--max-status-per-sec`；`--` 后参数可透传 Snakemake；常用长参数支持 `--option=value`。
- `batch_correction: "T"` 时，启动阶段自动要求样本表包含非空 `batch` 列。

### Changed（变更）
- `validate_samples.py` 改为英文 `argparse` CLI，增加结构化 summary、`--require-batch` 和 `--strict-warnings`，同时兼容旧的 `samples.csv control` 调用。
- `workflow/Snakefile`、`workflow/rules/common.smk` 与三套 profile 的运行时提示/核心注释统一为英文。
- README 与使用说明同步新的调度逻辑与启动参数。

### Fixed（修复）
- 修复 Snakemake 7.32.4 将命名 input/output 字段 `count` 视为保留名而导致 workflow parse 失败的问题；相关字段改名后四种 pipeline 均可正常解析。
- 修复 PE 样本 R1 文件（如 `ctrl_1_1.fastq.gz`）被 SE trim 规则误解释为独立样本而触发 `AmbiguousRuleException` 的问题；现在所有 `{sample}` wildcard 严格限制为样本表中的真实 ID。

## [0.4.1] - 2026-09-03

审查复盘修复：补齐配置化最后缺口 + 若干工程细节。三处行为修复均向后兼容（新 config 键缺省时行为与 v0.4.0 完全一致）。

### Fixed（修复）
- **DEGgroupCompare 物种参数断链**：`rules/deg.smk` 此前未把 `species`/`orgdb_tarball` 传入 `DEGgroupCompare.sh`，导致 `run_deg_compare.R` 恒按默认 osa 分支运行——`species=hsa` 时组间比较的 GO 富集会静默失败或错用 OrgDb；osa 则依赖上游富集步骤先安装 OrgDb 的副作用。现两条参数经规则 → 脚本 → 任务列表完整贯通（每任务显式引用传递）。
- `DEGgroupCompare.sh`：交集注释的临时文件由运行目录 `tmp.gene` 改为 `mktemp` + `trap` 清理（fgrep 无命中时不再残留）；并行任务部分失败时输出显式 WARN（不再完全静默，容错语义保持不变）。

### Added（新增）
- `tests/test_deggroupcompare.sh`：DEGgroupCompare 参数贯通回归测试（Rscript shim 记录参数，无需 R/snakemake，bash + python 即可运行）。
- config 新键 `batch_correction`（默认 `"F"`，与旧行为一致；设 `"T"` 时 DESeq2 按 `~ batch + group` 设计，样本表需含第三列 `batch`）。
- config 新键 `lncrna.pfam_scan`（默认 `"pfam_scan.pl"` 走 PATH 查找；pfam_scan.pl 成为最后一个可配置路径的外部工具）。

### Changed（变更）
- 注释修正：`lncRNA_functions.sh`、`envs/{lncrna,qc,assembly}.yaml` 中残留的 v0.2 时代旧文件名引用（`config_lncRNA.yaml`、`RNA-seq_lncRNA_DenovoIdenti.smk`、`multiQC_cleaned` 等）更新为现行 `config/config.yaml` lncrna 段与 `rules/*.smk` 命名。

## [0.4.0] - 2026-09-03

阶段 3（可持续性：测试、CI、文档随代码走）完成，详见 `docs/优化路线图.md`。三阶段改造计划至此全部落地。

### Added（新增）
- `tests/make_testdata.py`：确定性微型测试数据生成器（2 条 100kb 染色体、60 个模拟基因的参考/GTF/BED12/注释表，2 组 × 2 样本模拟 PE reads——含测序错误、10% 接头、5% 噪声 reads，以及预期差异基因 truth 表）；纯标准库实现，数据不入库、测试机即时生成（3.1）。
- `tests/check_outputs.py`：端到端输出断言器——关键文件存在性、mapping_stat 对齐/数值校验（ID 完整、Total_Reads=输入读对数、比对率可解析）、count 矩阵形状、DEG 结果方向校验（对 truth 上/下调基因）、图表与 sessionInfo 产出（3.1）。
- `tests/run_test.sh`：一键回归（生成数据 → dry-run → 端到端 → 断言 → DAG 再生成，`--pipeline`/`--reads`/`--keep` 可配）（3.1）。
- `tests/lint.sh`：静态检查聚合——bash -n、shellcheck、Python 编译、R 解析、snakemake --lint（4 种 pipeline），缺工具自动跳过（3.2）。
- `.github/workflows/ci.yml`：GitHub Actions CI（lint job + miniconda 端到端回归 job + DAG artifact）（3.2）。
- `workflow/profile/slurm/`：SLURM 调度 profile（3.3）。
- `CONTRIBUTING.md`：贡献指南——分支/提交信息规范、CHANGELOG 要求、文档同步 checklist、测试要求（3.4）。

### Changed（变更）
- `DEGgroupCompare.sh`：xargs 增加 `-r`，处理组不足 2 个（任务列表为空）时不再误执行。
- `envs/enrich.yaml`：新增 `bioconductor-org.hs.eg.db`，species=hsa 的富集开箱即用。

## [0.3.0] - 2026-09-03

阶段 2（结构标准化与质量提升）完成，详见 `docs/优化路线图.md`。**注意：旧入口 `rna-seq-workflow/RNA-seq_*.smk` 自本版本起废弃**（保留弃用提示 stub），请改用 `run.sh` / `workflow/Snakefile`。

### Added（新增）
- `workflow/Snakefile` 唯一入口：`pipeline=upstream|deg|as|lncrna` 经 `RNASEQ_PIPELINE` 环境变量注入（解析期确定），项目配置经 `RNASEQ_CONFIG` 注入；未设置时使用仓库默认配置（供 --lint 与示例）。
- `rules/` 模块化：`common`（辅助函数与 STAR 参数）/ `align`（trim/QC/STAR/链型）/ `quant`（定量合并）/ `deg` / `as` / `lncrna`；输出根目录 `results_dir` 可配置（默认 `results/`）。
- 全流程 MultiQC 汇总（FastQC/Trim Galore/STAR/featureCounts）+ `workflow/multiqc_config.yaml`（P2-9）。
- `software_versions.yaml`：每次运行自动记录环境定义、snakemake 版本与 git 提交（`collect_versions.py`）；DESeq2 输出附带 `sessionInfo.txt`（P2-9）。
- `config/` 配置体系（P2-4/P1-11）：`config.yaml` 主配置模板（全键注释，lncRNA 嵌套段）、`species.yaml`（osa/hsa 物种资源预设映射）、`samples.csv` 模板。
- `example/` 示例项目：配置模板 + 真实项目样本表（231107XTL 27 样本，自 `rna-seq-workflow/sample_info.csv` 迁出，P2-3）。
- `run.sh` 移至仓库根目录。

### Changed（变更）
- 统一两套 trim 规则（P2-2）：原始 fastq 双命名约定（`.fastq.gz`/`.fq.gz`）精确探测在全管线通用；trim 统一运行 FastQC（修复 P1-1 MultiQC 空报告）；trim 报告重命名为以样本 id 为键的规范名并声明为规则输出（P1-6）。
- STAR（P2-2）：直接输出 `BAM SortedByCoordinate`（P1-4）；索引改为文件级追踪（P1-5）；组装管线沿用组装优化参数、其余管线标准参数，均可经 `star_extra_args` 追加。
- 图表修复（P1-8）：火山图动态坐标范围（不再截断）、`fontface`、显式颜色映射、空类别容错；PCA `ntop` 可配置（`pca_ntop`，默认 20000）。
- 差异结果目录拼写修正：`Diff_Expr_Analysis_Reults` → `Diff_Expr_Analysis_Results`（`enrich.sh` 同步引用，P2-7/P2-8）。
- 脚本稳定化命名（去日期/连字符后缀）：`run_deseq2.R`、`run_enrichment.R`、`run_gsea.R`、`run_deg_compare.R`、`run_featurecounts.R`。
- KEGG 物种代码经 `species.yaml`（`kegg_organism`）传入富集脚本；`run_deg_compare.R` 输出位置改为随 results_dir。
- lncRNA 管线线程、外部工具路径全部经 config（`lncrna` 段）读取（P1-7）；CPC2/CNCI/Pfam/NR 关键产物声明为规则输出（P1-6）。
- `run_deseq2.R` 对照判定沿用 v0.2 的精确匹配；样本表 `group` 列按表头取值。

### Removed（移除）
- 旧四入口 Snakefile、旧规则文件（`rules/RNA-seq_*.smk`）、旧 `config_*.yaml`（被 `config/` 体系替代）、过期 DAG 图（阶段 3 CI 重新生成）；旧入口位置保留弃用提示 stub。

## [0.2.0] - 2026-09-03

阶段 1（P0 修复）完成，详见 `docs/优化路线图.md`。

### Added（新增）
- `envs/*.yaml`：7 个 conda 环境定义（qc / align / quant / deseq2 / enrich / assembly / lncrna），并为各规则挂接 `conda:` 指令（P0-2）。
- `profile/sge`、`profile/default`：集群与本地调度 profile。
- `run.sh`：统一启动脚本——自动选择调度器、配置优先级解析、运行前样本表校验；替代已删除的 `main_runRNA-seq.sh`。
- `scripts/merge_featurecounts.py`：自包含的定量结果合并脚本（替代依赖仓库外 `njoin.sh`/`transposition.sh` 的 `featureCount.R_result_merge.sh` 与空壳 `.py`，P0-3/P0-4）。
- `scripts/validate_samples.py`：样本表校验器（重复 id、`-` 字符、前缀冲突、对照组存在性、layout 取值）。
- `sample_info.csv` 增加 `layout` 列（PE/SE/auto）。

### Changed（变更）
- 路径参数化（P0-1/P0-8）：入口 Snakefile 的 `configfile:`/`include:` 改为基于 `workflow.basedir`；规则内脚本调用统一走 `SCRIPTS` 变量；R 脚本移除 `.libPaths` 硬编码；配置中机器相关路径改为 `/path/to/...` 占位符。
- `Mapping_stat` 重写（P0-6）：逐样本解析 `*_Log.final.out` 与 trim 报告，样本 id 不再截断、双端样本不再错位。
- 原始 fastq 探测改为精确文件名匹配（不再 `ls {sample}*` 通配，消除前缀互配风险）；`runSTAR --readFilesIn` 使用显式文件；featureCounts 的链型/双端判定不再依赖 `ls` 计数。
- 对照组配置化（P0-7）：新增 `control_group` 配置键，贯通 runDESeq2（`-r`，对照判定由 `grepl` 子串匹配改为精确匹配）/ enrich / getGroups.py / DEGgroupCompare。
- 物种配置化（P0-8）：新增 `species`/`annotation_tsv`/`orgdb_tarball` 配置键；富集 R 脚本支持 `osa`/`hsa` 双物种；KEGG 查询失败时跳过而不中断流程。
- 组间比较并行化：`ParaFly` 替换为 `xargs -P`；`getGroups.py` 改按表头取 `group` 列。
- lncRNA 管线线程数配置化（`lncrna_threads`，原硬编码 30）；CPC2/CNCI/Pfam 路径经 config 注入 `lncRNA_functions.sh`。
- `count_merge`/`count_merge2` 显式声明 `GeneExpression_FPKM.xls`、`GeneCount_Assigned_logs.xls` 与 `{sample}.log` 输入输出；`check_strandedness` 声明 `_infer_experiment.out` 输出。
- `STAR_index` 修复损坏的 stderr 重定向（`; 2>{log}` → `>> {log} 2>&1`）；`runSTAR` 移除从不清理的 `--outReadsUnmapped Fastx` 产物。

### Removed（移除）
- 死/坏脚本：`enrich_KEGG_clusterProfiler.R`（语法错误）、`check_strandness.sh`（未展开模板）、`getExpr_featureCounts.sh`（截断）、`featureCount.R_result_merge.{sh,py}`（被新合并脚本替代）、`enrich_GO_KEGG_clusterProfiler_gProfilerGO_hsa.R`（并入主脚本物种分支）、`test.py`/`test.groups`、空的 `SampleListFile`。

### Known issues（遗留问题）
- upstream/DEG 管线 trim 规则未跑 FastQC（MultiQC 清洗后报告为空）；KEGG 依赖外网；火山图坐标硬编码等 P1/P2 问题留待阶段 2（见 `docs/审查报告.md`）。

## [0.1.0] - 2026-09-03

### Added（新增）
- 建立 git 仓库，原样收录 `rna-seq-workflow/` 流程代码作为基线（Snakemake 四入口 + rules + scripts + 配置）。
- 标准项目脚手架：
  - `README.md`（项目说明与文档索引）
  - `docs/使用说明.md`（操作手册：数据准备/配置/运行/结果解读/FAQ/断点清单）
  - `docs/审查报告.md`（全面代码审查：P0×8 / P1×14 / P2×9，含 file:line 定位与方法学评估）
  - `docs/优化路线图.md`（三阶段改造计划、验收标准、方案对比）
  - `CHANGELOG.md`、`LICENSE`
- `.gitattributes`（脚本统一 LF 行尾，面向 Linux 集群运行）、`.gitignore`（忽略 Snakemake 运行产物）。

### Unchanged（未改动）
- 流程代码本体未做任何修改；全部修复工作按 `docs/优化路线图.md` 阶段 1 起推进。

### Known issues（已知问题）
- 详见 `docs/审查报告.md`：硬编码绝对路径、conda 环境文件缺失、count 合并外部依赖缺失、lncRNA 管线终点空壳脚本、`mapping_stat.xls` 统计错位、"control" 组名硬编码、富集环节物种硬编码。
