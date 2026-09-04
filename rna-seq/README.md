# RNA-seq 转录组分析工作流

基于 **Snakemake** 的 bulk RNA-seq（有参考基因组）端到端分析流程，覆盖从原始 fastq 到差异表达与功能富集的完整链路，并支持可变剪接分析与 de novo lncRNA 鉴定。

> **状态**：v0.4.0——阶段 1（P0 修复）、阶段 2（标准布局重构）、阶段 3（测试与 CI）已完成：单入口 `workflow/Snakefile`（`pipeline=upstream|deg|as|lncrna`）、conda 环境入库、路径全参数化、全流程 MultiQC、`software_versions.yaml`、一键回归测试（`tests/run_test.sh`）与静态检查（`tests/lint.sh`）。遗留事项见 [docs/审查报告.md](docs/审查报告.md)。

## 功能总览

| pipeline（经 `run.sh` 第一个参数选择） | 内容 |
|---|---|
| `upstream` | FastQC → Trim Galore（含 FastQC）→ 链特异性推断 → STAR 比对 → featureCounts 定量 → TPM/count 矩阵 |
| `deg` | upstream + DESeq2 差异表达 → GO/KEGG 富集（clusterProfiler + aPEAR）→ GSEA → 组间比较 |
| `as` | upstream（STAR 组装优化参数）+ StringTie 转录本组装/合并/gffcompare → isoform 定量 |
| `lncrna` | upstream + StringTie 组装 → 编码势过滤（CPC2/CNCI/长度/Pfam/NR）→ lncRNA 表达矩阵 |

运行环境：Linux + SGE 集群（qsub，可选）+ conda + Snakemake ≥7；方法主线：STAR two-pass → featureCounts → DESeq2 → clusterProfiler。

## 目录结构

```
rna-seq/
├── run.sh                    # 统一启动脚本（SGE/SLURM/本地 profile 自动或手动选择）
├── workflow/
│   ├── Snakefile             # 唯一入口（pipeline 参数选择模式）
│   ├── rules/                # common / align / quant / deg / as / lncrna
│   ├── scripts/              # R/Python/Shell 脚本（稳定命名，无日期后缀）
│   ├── envs/                 # conda 环境定义（7 个，--use-conda 自动构建）
│   ├── profile/              # 调度 profile（sge / slurm / default）
│   └── multiqc_config.yaml
├── config/
│   ├── config.yaml           # 主配置模板（全键注释，含 lncRNA 段）
│   ├── species.yaml          # osa/hsa 物种资源预设
│   └── samples.csv           # 样本表模板
├── example/                  # 示例项目（真实样本表参考 + 一条命令上手）
├── tests/                    # 回归测试：make_testdata.py（合成数据）/ run_test.sh / lint.sh
├── docs/                     # 使用说明 / 审查报告 / 优化路线图
├── CONTRIBUTING.md           # 贡献指南（文档同步 checklist / CHANGELOG 要求）
└── results/                  # 运行产物（gitignore；根目录可配 results_dir）
```

## 快速开始

完整步骤见 [docs/使用说明.md](docs/使用说明.md)，概要：

1. **建项目**：项目目录放 `1.rawdata/`（`{sample}_1/_2.fastq.gz` 或单端 `{sample}.fastq.gz`）与样本表（模板 `config/samples.csv`，须含对照组）。
2. **写配置**：复制 `config/config.yaml` 到项目目录，替换 `/path/to/...` 占位符（未给出的资源键回落到 `config/species.yaml` 物种预设）。
3. **运行**：`bash run.sh deg /path/to/myproject myproject/config.yaml 10`（建议先 `snakemake -n` dry-run）。

> 首次运行自动构建 7 个 conda 环境；每次运行产出 `results/software_versions.yaml` 与全流程 `multiqc_report.html`。CPC2/CNCI/pfam_scan.pl 等外部工具在 config 的 `lncrna` 段配置。

## 文档索引

| 文档 | 内容 |
|---|---|
| [docs/使用说明.md](docs/使用说明.md) | 数据准备、样本表约束、配置说明、运行监控、结果目录解读、FAQ |
| [docs/审查报告.md](docs/审查报告.md) | 2026-09 全面审查：8 个阻断性问题（P0）+ 14 个重要缺陷（P1）+ 9 个规范问题（P2），含 file:line 定位 |
| [docs/优化路线图.md](docs/优化路线图.md) | 三阶段改造计划与验收标准（阶段 0-3 全部完成） |
| [CHANGELOG.md](CHANGELOG.md) | 版本变更记录 |

## 已知问题（v0.4.0 修复情况）

阶段 1 已修复全部 8 个 P0；阶段 2 完成结构标准化、统一 trim/FastQC（MultiQC 空报告已修复）、火山图修复、`results_dir` 输出根目录参数化、物种预设映射、版本记录；阶段 3 提供一键回归测试、静态检查、CI 配置与 slurm profile。

仍遗留的主要事项：

1. KEGG 富集依赖外网，离线节点自动跳过（附录 A 提供离线 gmt 备选方案）。
2. 批次校正参数 `-b` 未配置化。
3. `lncrna` 管线未纳入自动回归测试（依赖 CPC2/CNCI/pfam_scan.pl 外部工具）。
4. CI 需托管到 GitHub/Gitee 后启用（`.github/workflows/ci.yml` 已就绪）；端到端回归与 lint 也可本地运行（`tests/`）。

## 测试

```bash
bash tests/lint.sh                        # 静态检查：bash/shellcheck/Python/R/snakemake --lint
bash tests/run_test.sh --pipeline deg     # 端到端回归：生成合成数据 → dry-run → 运行 → 35 项断言
```

测试数据由 `tests/make_testdata.py` 确定性生成（2 条 100kb 染色体、60 个模拟基因、2 组 × 2 样本模拟 reads，含接头/测序错误/预期差异基因 truth），不入库。断言覆盖关键文件存在性、mapping_stat 对齐与数值校验、count 矩阵形状、DEG 方向校验（对 truth 上/下调基因）。

## 贡献与变更

- 修改流程代码请先阅读 [docs/优化路线图.md](docs/优化路线图.md)，避免与改造方向冲突；所有变更记入 `CHANGELOG.md`。
- 行尾统一 LF（见 `.gitattributes`），运行产物不入库（见 `.gitignore`）。

## License

[MIT](LICENSE)
