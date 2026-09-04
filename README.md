# BioWorkflows

个人生物信息学 Snakemake 工作流集合，目前包含两个独立的工作流项目：

| 目录 | 简介 |
|------|------|
| [`rna-seq/`](rna-seq/) | bulk RNA-seq（有参考基因组）端到端分析流程：从原始 fastq 到差异表达与功能富集，支持可变剪接分析与 de novo lncRNA 鉴定 |
| [`chip_cuttag_atac_faire/`](chip_cuttag_atac_faire/) | 基于 Snakemake 的植物表观组学一站式分析流程，单入口 `workflow/Snakefile` 同时支持 ChIP-seq / CUT&Tag / ATAC-seq / FAIRE-seq 四种数据类型（可混型项目） |

每个子目录均为独立工作流项目，自带配置模板、文档（`docs/`）、示例数据（`example/`）、
一键回归测试（`tests/`）与运行脚本，详见各自目录下的 `README.md` 与 `docs/`。

## 许可证

本仓库整体以 [Apache-2.0](LICENSE) 发布；`rna-seq/` 与 `chip_cuttag_atac_faire/`
子目录的代码保留其原 MIT 许可声明（见各子目录 `LICENSE` 文件）。
