# chip_cuttag_atac_faire 后续待做项（TODO）

> 产生于 enhancer_lncRNA_2026 项目实测准备阶段（2026-09-05），均为已识别、暂缓实施的事项。

## 1. FASTQ 命名约定兼容

- **现状**：`workflow/rules/upstream.smk` 的 fastq 输入仅识别 `1.rawdata/{sample}_1.fq.gz + {sample}_2.fq.gz`；测序交付普遍使用的 `{sample}_R1/_R2.fq.gz`、`.fastq.gz` 后缀变体不被识别。
- **当前临时方案**：在项目 `1.rawdata/` 内补一套 `{sample}_1.fq.gz`/`{sample}_2.fq.gz` 软链接指向同一批原始文件（纯新增，不动原始数据）。
- **待做**：在样本表 → fastq 路径解析处增加 `_R1/_R2` 与 `.fastq.gz` 变体支持（与 rna-seq 工作流的同类待办保持一致，见 `rna-seq/docs/TODO.md`），或在用户文档中显著说明命名约定。

## 2. `--validate-only` 与锁定的 Snakemake 7 不兼容

- **现状**：`run.sh` 的 `--validate-only` 调用 `snakemake --list-rules`，该参数是 Snakemake 8 的改名（7.x 为 `--list`）；在锁定的 snakemake-minimal=7.32.4 下报 `unrecognized arguments: --list-rules`。
- **待做**：改回 7.x 兼容的 `--list`（或做版本分支）；日常验证用 `-n` dry-run 可完全替代。

## 3. 规则内裸调 `Rscript` 绕过 `r.rscript` 配置

- **现状**：`annoPeak_batch.R` 等规则直接写 `Rscript ...`（见 `callpeak/annotation` 规则），不经过
  `software.yaml` 的 `r.rscript`/`CHIP_RSCRIPT` 通道。当项目配置指向独立的 R（如 R4.2.3 包装脚本）
  时，作业里裸调的 `Rscript` 仍会解析到主环境自带的 R（enhancer_lncRNA_2026 部署中为 chip env 的
  R 4.1.3 + 高版本 R 库 → `rlang.so: undefined symbol: EXTPTR_PROT`，peak_annotation 失败）。
- **当前绕过**：项目 `bin/env.sh`（BASH_ENV 注入）把含 Rscript 包装脚本的 bin 目录放到 PATH 最前。
- **待做**：规则改用 `tool("rscript")` 或 `${CHIP_RSCRIPT:-Rscript}`，与 rna-seq 的做法对齐。
