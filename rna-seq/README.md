# RNA-seq 转录组分析工作流

基于 **Snakemake** 的 bulk RNA-seq（有参考基因组）端到端分析流程，覆盖从原始 fastq 到差异表达与功能富集的完整链路，并支持可变剪接分析与 de novo lncRNA 鉴定。

> **状态**：v0.8.0——运行环境改为“单一主环境 + `software.yaml` 覆盖”模式：可直接复用服务器已有 Conda 环境或系统 PATH，并统一配置 Rscript、R 版本、R library path、外部工具与数据库；Snakemake 不再为每个 rule 自动创建独立 Conda 环境。v0.7 的 per-rule CPU/内存/walltime 调度模型继续保留。

## 功能总览

| pipeline（经 `run.sh` 第一个参数选择） | 内容 |
|---|---|
| `upstream` | FastQC → Trim Galore（含 FastQC）→ 链特异性推断 → STAR 比对 → featureCounts 定量 → TPM/count 矩阵 |
| `deg` | upstream + DESeq2 差异表达 → GO/KEGG 富集（clusterProfiler + aPEAR）→ GSEA → 组间比较 |
| `as` | upstream（STAR 组装优化参数）+ StringTie 转录本组装/合并/gffcompare → isoform 定量 |
| `lncrna` | upstream + StringTie 组装 → 编码势过滤（CPC2/CNCI/长度/Pfam/NR）→ lncRNA 表达矩阵 |

运行环境：Linux + local / SGE（qsub）/ SLURM（sbatch）+ Snakemake ≥7。推荐在服务器维护一个统一的 `rna-seq` Conda 环境，由 `config/software.yaml` 指定其 prefix/name；也可直接使用系统 PATH。`run.sh` 负责解析 runtime、校验软件/R 包并按 SGE → SLURM → local 选择调度 profile。

调度资源按 rule 独立声明；例如 STAR index 默认 12 threads / 48 GB / 480 min，STAR alignment 默认 12 threads / 32 GB / 360 min，DESeq2 默认 4 threads / 16 GB / 240 min。项目可在 `config.yaml` 的 `resources:` 段覆盖，`run.sh --memory` / `--runtime` 则用于全局强制覆盖集群请求。

## 目录结构

```
rna-seq/
├── run.sh                    # 统一启动脚本（SGE/SLURM/本地 profile 自动或手动选择）
├── workflow/
│   ├── Snakefile             # 唯一入口（pipeline 参数选择模式）
│   ├── rules/                # common / meta / align / quant / deg / as / lncrna
│   ├── scripts/              # R/Python/Shell 脚本（稳定命名，无日期后缀）
│   ├── environment.yaml      # 可选的一体化 Conda 环境模板（不自动创建）
│   ├── profile/              # 调度 profile（sge / slurm / default）
│   └── multiqc_config.yaml
├── config/
│   ├── config.yaml           # 分析参数与任务资源
│   ├── software.yaml         # 主环境、R runtime、工具/数据库覆盖
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
2. **写配置**：复制 `config/config.yaml` 与 `config/software.yaml` 到项目目录。前者配置分析/参考资源，后者配置服务器软件环境；已有统一 Conda 环境时通常只需填写 `environment.conda_prefix`。
3. **检查环境**：`bash run.sh -p deg -P /path/to/myproject --software /path/to/myproject/software.yaml --check-software`。
4. **运行**：`bash run.sh deg /path/to/myproject myproject/config.yaml 10`；建议先用 `bash run.sh -p deg -P /path/to/myproject -c /path/to/myproject/config.yaml --dry-run` 预览 DAG。

> Snakemake 不再自动创建 Conda 环境。没有现成环境时，可用 `mamba env create -f workflow/environment.yaml` 创建推荐的一体化环境；已有环境则直接复用。每次运行的实际环境、Rscript/R library 和工具解析结果会记录到 `results/software_versions.yaml`。

## 文档索引

| 文档 | 内容 |
|---|---|
| [docs/使用说明.md](docs/使用说明.md) | 数据准备、样本表约束、配置说明、运行监控、结果目录解读、FAQ |
| [docs/审查报告.md](docs/审查报告.md) | 2026-09 全面审查：8 个阻断性问题（P0）+ 14 个重要缺陷（P1）+ 9 个规范问题（P2），含 file:line 定位 |
| [docs/优化路线图.md](docs/优化路线图.md) | 三阶段改造计划与验收标准（阶段 0-3 全部完成） |
| [CHANGELOG.md](CHANGELOG.md) | 版本变更记录 |

## 已知问题（v0.8.0）

阶段 1 已修复全部 8 个 P0；阶段 2 完成结构标准化、统一 trim/FastQC（MultiQC 空报告已修复）、火山图修复、`results_dir` 输出根目录参数化、物种预设映射、版本记录；阶段 3 提供一键回归测试、静态检查、CI 配置与 slurm profile；v0.8.0 进一步将软件/R 环境从分析参数中分离为统一 `software.yaml` runtime。

仍遗留的主要事项：

1. KEGG 富集依赖外网，离线节点自动跳过（附录 A 提供离线 gmt 备选方案）。
2. `lncrna` 管线未纳入完整端到端自动回归测试（依赖 CPC2/CNCI/pfam_scan.pl 与大型 Pfam/NR 数据库）；这些依赖通过 `software.yaml` 配置，CNCI 的 Python2 解释器可单独指定。
3. CI 需托管到 GitHub/Gitee 后启用（`.github/workflows/ci.yml` 已就绪）；端到端回归与 lint 也可本地运行（`tests/`）。

## 测试

```bash
bash tests/lint.sh                        # 静态检查：bash/shellcheck/Python/R/snakemake --lint
bash tests/test_deggroupcompare.sh        # DEGgroupCompare 参数贯通回归（无需 R/snakemake）
bash tests/test_runtime_config.sh         # Conda prefix + Rscript/Rlib 解析回归
bash tests/test_r_runtime_propagation.sh # 嵌套 shell/R 子任务继承 Rscript/R_LIBS
bash tests/run_test.sh --pipeline deg     # 端到端回归：生成合成数据 → dry-run → 运行 → 35 项断言
```

测试数据由 `tests/make_testdata.py` 确定性生成（2 条 100kb 染色体、60 个模拟基因、2 组 × 2 样本模拟 reads，含接头/测序错误/预期差异基因 truth），不入库。断言覆盖关键文件存在性、mapping_stat 对齐与数值校验、count 矩阵形状、DEG 方向校验（对 truth 上/下调基因）。

## 贡献与变更

- 修改流程代码请先阅读 [docs/优化路线图.md](docs/优化路线图.md)，避免与改造方向冲突；所有变更记入 `CHANGELOG.md`。
- 行尾统一 LF（见 `.gitattributes`），运行产物不入库（见 `.gitignore`）。

## License

[MIT](LICENSE)
