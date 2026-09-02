# legacy/

历史版本归档（v0.1.0 原始实现，2024-03）。这些文件已被 `workflow.smk` + `rules/` 的新实现取代，**不再维护、不可直接运行**，仅作历史参考保留。

| 归档文件 | 原位置 | 取代者 |
|---|---|---|
| `chip-seq.smk` / `atac.smk` / `cuttag.smk` / `faire.smk` | 根目录 | `workflow.smk`（按 seqtype 自动路由） |
| `upsteam.smk` | `rules/` | `rules/upstream.smk`（修正拼写） |
| `rmDup.smk` | `rules/` | `rules/dedup.smk` |
| `callpeak_chip/atac/cuttag/faire.smk` | `rules/` | `rules/callpeak.smk`（按组并行规则） |
| `chip_environment.yaml` | 根目录 | `envs/*.yaml`（per-rule 环境） |
| `call_peak.sh`（根目录与 scripts/ 各一份） | — | `rules/callpeak.smk` |
| `bdgcmp_macs2.sh` | `scripts/` | `rules/callpeak.smk` 中的 bigwig 规则（排序改用 bedtools sort -g） |
| `run_deeptools_QC.sh` | `scripts/` | `rules/qc_deeptools.smk`（全量规则化 + FRiP） |

仍在 `scripts/` 维持的文件（不属于 legacy）：`annoPeak_batch.R`（已参数化重写）、`annoPeak_single.R`（已重写）、`run_ChIPQC.R`（已修复参数化）、`run_chipqc_DROMPAplus.sh`（已参数化）；`diffpeak_DiffBind.sh` / `run_DiffBind.R` 为待实现的空壳（需 contrast 设计决策）。
