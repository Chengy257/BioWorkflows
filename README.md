# RNA-seq 转录组分析工作流

基于 **Snakemake** 的 bulk RNA-seq（有参考基因组）端到端分析流程，覆盖从原始 fastq 到差异表达与功能富集的完整链路，并支持可变剪接分析与 de novo lncRNA 鉴定。

> **状态提示**：项目正处于标准化改造初期。流程方法学主线可用，但存在多项待修复的工程问题（硬编码路径、conda 环境缺失等），详见 [docs/审查报告.md](docs/审查报告.md) 与 [docs/优化路线图.md](docs/优化路线图.md)。

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
└── rna-seq-workflow/          # 流程代码（原样基线，待阶段 1-2 改造）
    ├── main_runRNA-seq.sh     # 主控脚本（SGE 调度）
    ├── RNA-seq_*.smk          # 四个管线入口
    ├── rules/                 # Snakemake 规则
    ├── scripts/               # R/Python/Shell 脚本
    ├── config_*.yaml          # 配置文件
    └── sample_info.csv        # 样本分组表（示例）
```

## 快速开始

完整步骤见 [docs/使用说明.md](docs/使用说明.md)，概要：

1. **准备数据**：工作目录下放置 `1.rawdata/{sample}_1.fastq.gz`、`{sample}_2.fastq.gz`（双端）及 `sample_info.csv`（须含 `control` 组）。
2. **写配置**：复制 `config_user_defined.yaml` 模板，填入 `gtf`/`genome`/`threads` 等。
3. **改主控脚本**：核对 `main_runRNA-seq.sh` 中的 `smk`/`DIR`/`config`/`conda_path` 路径。
4. **运行**：`bash main_runRNA-seq.sh`（建议先 `snakemake -n` dry-run）。

⚠️ **当前版本运行前请先通读 [使用说明 §6 断点清单](docs/使用说明.md#6-已知限制与断点清单阶段-1-修复前必须知晓)**——conda 环境文件与部分依赖脚本尚不在仓库中。

## 文档索引

| 文档 | 内容 |
|---|---|
| [docs/使用说明.md](docs/使用说明.md) | 数据准备、样本表约束、配置说明、运行监控、结果目录解读、FAQ |
| [docs/审查报告.md](docs/审查报告.md) | 2026-09 全面审查：8 个阻断性问题（P0）+ 14 个重要缺陷（P1）+ 9 个规范问题（P2），含 file:line 定位 |
| [docs/优化路线图.md](docs/优化路线图.md) | 阶段 1 修复阻断问题 → 阶段 2 标准布局重构 → 阶段 3 测试与 CI，附验收标准与方案对比 |
| [CHANGELOG.md](CHANGELOG.md) | 版本变更记录 |

## 已知问题 Top 5

1. **路径硬编码**：入口/规则/脚本中大量 `/home/chengyu/...` 与 `/share/...` 绝对路径，换环境必挂（P0-1）。
2. **conda 环境缺失**：规则引用 `rna-seq`/`common`/`lncrna` 环境，但 `envs/*.yaml` 不在仓库，`--use-conda` 无法工作（P0-2）。
3. **count 合并依赖缺失脚本**：`njoin.sh`/`transposition.sh` 未入库，DEG 链路断（P0-3）；lncRNA 管线终点调用空壳脚本（P0-4）。
4. **比对统计表对实际样本数据即出错**：`mapping_stat.xls` 截断含下划线样本 ID、双端样本行错位（P0-6）。
5. **"control" 组名硬编码**：无名为 control 的分组时差异分析静默无结果（P0-7）。

## 贡献与变更

- 修改流程代码请先阅读 [docs/优化路线图.md](docs/优化路线图.md)，避免与改造方向冲突；所有变更记入 `CHANGELOG.md`。
- 行尾统一 LF（见 `.gitattributes`），运行产物不入库（见 `.gitignore`）。

## License

[MIT](LICENSE)
