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

## legacy/diffbind/（v0.4.0 归档）

未接入主流水线 DAG 的独立 QC / 差异分析空壳脚本，v0.4.0 起从 `scripts/` 归档至此，**不可直接运行**：

| 归档文件 | 原位置 | 取代者 / 恢复路线 |
|---|---|---|
| `annoPeak_single.R` | `scripts/` | `workflow/scripts/annoPeak_batch.R`（已参数化重写，按组批量注释） |
| `run_ChIPQC.R` | `scripts/` | `rules/qc_deeptools.smk` + FRiP（QC 已规则化接入 DAG）；ChIPQC 恢复需先修复 docs/REVIEW.md 所列参数化问题 |
| `run_chipqc_DROMPAplus.sh` | `scripts/` | `rules/qc_deeptools.smk` + FRiP（同上；DROMPAplus 依赖 docker 镜像未参数化，见 docs/REVIEW.md） |
| `diffpeak_DiffBind.sh` | `scripts/` | 暂无——差异分析待实现（contrast 与设计公式待定，见 docs/REVIEW.md） |
| `run_DiffBind.R` | `scripts/` | 同上（空壳文件，仅有文件头） |

仍在 `workflow/scripts/` 维持的文件（不属于 legacy）：`annoPeak_batch.R`（已参数化重写）、`runtime_config.py`、`collect_versions.py`。
