# Changelog

本项目的全部显著变更记录于此。格式参考 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/)。

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
