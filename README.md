# RNA-seq 转录组分析工作流

基于 **Snakemake** 的 bulk RNA-seq（有参考基因组）端到端分析流程，覆盖从原始 fastq 到差异表达与功能富集的完整链路，并支持可变剪接分析与 de novo lncRNA 鉴定。

> **状态提示**：v0.2.0 已完成阶段 1 修复——路径参数化、conda 环境入库（`envs/*.yaml`）、比对统计修正、对照组/物种配置化。首次运行前需将配置中的 `/path/to/...` 占位符替换为真实路径。遗留 P1/P2 问题与后续计划见 [docs/审查报告.md](docs/审查报告.md) 与 [docs/优化路线图.md](docs/优化路线图.md)。

## 功能总览

| 管线入口（`rna-seq-workflow/`） | 内容 |
|---|---|
| `RNA-seq_upstream_only.smk` | FastQC → Trim Galore → 链特异性推断 → STAR 比对 → featureCounts 定量 → TPM/count 矩阵 |
| `RNA-seq_up_DEG.smk` | 上游 + DESeq2 差异表达 → GO/KEGG 富集（clusterProfiler + aPEAR）→ GSEA → 组间比较 |
| `RNA-seq_up_AS.smk` | 上游 + StringTie 转录本组装/合并/gffcompare → isoform 定量 |
| `RNA-seq_up_lncRNA.smk` | 上游 + StringTie 组装 → 编码势过滤（CPC2/CNCI/长度/Pfam/NR）→ lncRNA 表达矩阵 |

运行环境：Linux + SGE 集群（qsub）+ conda + Snakemake ≥7；方法主线：STAR two-pass → featureCounts → DESeq2 → clusterProfiler。

流程 DAG：[dag.svg](rna-seq-workflow/dag.svg) ｜ [dag_lnc.svg](rna-seq-workflow/dag_lnc.svg)

## 目录结构

```
rna-seq/
├── README.md                  # 本文件
├── CHANGELOG.md               # 变更记录
├── LICENSE
├── docs/
│   ├── 使用说明.md            # 详细操作手册（数据准备/配置/运行/结果解读/FAQ）
│   ├── 审查报告.md            # 全量代码审查：P0/P1/P2 分级问题清单
│   └── 优化路线图.md          # 三阶段改造计划与验收标准
└── rna-seq-workflow/          # 流程代码
    ├── run.sh                 # 统一启动脚本（自动选择 SGE/本地 profile）
    ├── RNA-seq_*.smk          # 四个管线入口
    ├── rules/                 # Snakemake 规则
    ├── scripts/               # R/Python/Shell 脚本
    ├── envs/                  # conda 环境定义（7 个，--use-conda 自动构建）
    ├── profile/               # 调度 profile（sge / default）
    ├── config_*.yaml          # 配置文件（默认配置 + 模板，路径需替换）
    └── sample_info.csv        # 样本分组表（真实项目示例，含 layout 列）
```

## 快速开始

完整步骤见 [docs/使用说明.md](docs/使用说明.md)，概要：

1. **准备数据**：工作目录下放置 `1.rawdata/{sample}_1.fastq.gz`、`{sample}_2.fastq.gz`（双端）及 `sample_info.csv`（含对照组，命名规范见使用说明）。
2. **写配置**：复制 `config_user_defined.yaml` 模板（或直接修改 `config_basic_defaulted.yaml`），把所有 `/path/to/...` 占位符换成集群真实路径。
3. **运行**：`bash rna-seq-workflow/run.sh deg /path/to/project my_config.yaml 10`（建议先 `snakemake -n` dry-run）。

> 首次运行会自动构建 7 个 conda 环境（耗时较长）；CPC2/CNCI/pfam_scan.pl 等外部工具路径在 `config_lncRNA.yaml` 配置。

## 文档索引

| 文档 | 内容 |
|---|---|
| [docs/使用说明.md](docs/使用说明.md) | 数据准备、样本表约束、配置说明、运行监控、结果目录解读、FAQ |
| [docs/审查报告.md](docs/审查报告.md) | 2026-09 全面审查：8 个阻断性问题（P0）+ 14 个重要缺陷（P1）+ 9 个规范问题（P2），含 file:line 定位 |
| [docs/优化路线图.md](docs/优化路线图.md) | 阶段 1 修复阻断问题 → 阶段 2 标准布局重构 → 阶段 3 测试与 CI，附验收标准与方案对比 |
| [CHANGELOG.md](CHANGELOG.md) | 版本变更记录 |

## 已知问题（v0.2.0 修复情况）

阶段 1 已修复全部 8 个 P0：路径参数化（~~P0-1~~）、conda 环境入库（~~P0-2~~）、合并脚本自包含（~~P0-3/P0-4~~）、死脚本清理（~~P0-5~~）、mapping_stat 修正（~~P0-6~~）、对照组配置化（~~P0-7~~）、物种配置化 + 去 ParaFly（~~P0-8~~）。

仍遗留的主要 P1/P2（详见 [审查报告](docs/审查报告.md)）：

1. upstream/DEG 管线 trim 规则未跑 FastQC，MultiQC 清洗后报告为空（P1-1，阶段 2 统一 trim 规则时修复）。
2. KEGG 富集依赖外网，离线节点自动跳过（阶段 2 提供离线方案）。
3. 火山图坐标硬编码、`frontface` 无效参数等图表问题（P1-8，阶段 2）。
4. 批次校正参数 `-b` 未配置化（阶段 2）。

## 贡献与变更

- 修改流程代码请先阅读 [docs/优化路线图.md](docs/优化路线图.md)，避免与改造方向冲突；所有变更记入 `CHANGELOG.md`。
- 行尾统一 LF（见 `.gitattributes`），运行产物不入库（见 `.gitignore`）。

## License

[MIT](LICENSE)
