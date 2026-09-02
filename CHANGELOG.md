# Changelog

本项目的全部显著变更记录于此。格式参考 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/)。

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
