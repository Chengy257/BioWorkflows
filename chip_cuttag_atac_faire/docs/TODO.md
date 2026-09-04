# chip_cuttag_atac_faire 后续待做项（TODO）

> 产生于 enhancer_lncRNA_2026 项目实测准备阶段（2026-09-05），均为已识别、暂缓实施的事项。

## 1. FASTQ 命名约定兼容

- **现状**：`workflow/rules/upstream.smk` 的 fastq 输入仅识别 `1.rawdata/{sample}_1.fq.gz + {sample}_2.fq.gz`；测序交付普遍使用的 `{sample}_R1/_R2.fq.gz`、`.fastq.gz` 后缀变体不被识别。
- **当前临时方案**：在项目 `1.rawdata/` 内补一套 `{sample}_1.fq.gz`/`{sample}_2.fq.gz` 软链接指向同一批原始文件（纯新增，不动原始数据）。
- **待做**：在样本表 → fastq 路径解析处增加 `_R1/_R2` 与 `.fastq.gz` 变体支持（与 rna-seq 工作流的同类待办保持一致，见 `rna-seq/docs/TODO.md`），或在用户文档中显著说明命名约定。

## 2. `--validate-only` 与锁定的 Snakemake 7 不兼容

- **现状**：`run.sh` 的 `--validate-only` 调用 `snakemake --list-rules`，该参数是 Snakemake 8 的改名（7.x 为 `--list`）；在锁定的 snakemake-minimal=7.32.4 下报 `unrecognized arguments: --list-rules`。
- **待做**：改回 7.x 兼容的 `--list`（或做版本分支）；日常验证用 `-n` dry-run 可完全替代。
