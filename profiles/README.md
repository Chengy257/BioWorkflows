# profiles/

Snakemake 运行 profile（把集群参数固化进版本库）。

| profile | 适用 | 说明 |
|---|---|---|
| `pbs/` | snakemake **7.x** | PBS/qsub；`--profile profiles/pbs` 使用；任务名取 `{rule}`、核数取 `{threads}`（与各规则声明自动对齐）；日志直接落入 `logs/` |

snakemake 8.x 改用 executor 插件体系（`--cluster` 语义移除），8.x 用户请使用仓库根目录的 `main_run.sh -p "qsub ..."`（已做 7/8 版本自动适配）。如后续需要在 8.x 上使用 profile，可安装 `snakemake-executor-plugin-cluster-generic` 后新增对应 profile。
