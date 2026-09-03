# workflow/profile/

Snakemake 运行 profile（把集群参数固化进版本库）。调度方式通过 `--profile workflow/profile/<name>` 选择；auto 探测规则见 `run.sh --help`（v0.4.0）。

| profile | 调度器 | cluster 串要点 |
|---|---|---|
| `default/` | 本机/独立服务器 | 无 cluster 串，`--cores 8` 本地直接执行；`latency-wait: 90` |
| `pbs/` | PBS (Torque)，snakemake **7.x** 经典接口 | `qsub -V -N {rule} -l select=1:ncpus={threads}:mem={resources.mem_mb}mb -l walltime={resources.runtime_sec} -j oe`；walltime 用秒避免 `[[HH:]MM:]SS` 歧义 |
| `sge/` | SGE，snakemake **7.x** 经典接口 | `qsub -V -N {rule} -l ncpus={threads} -l h_vmem={resources.mem_mb}M -l h_rt={resources.runtime_sec} -j oe` |
| `slurm/` | SLURM，snakemake **7.x** 经典接口 | `sbatch --parsable -J {rule} -c {threads} --mem={resources.mem_mb}M --time={resources.runtime_min} -o slurm-{rule}-%j.out` |

cluster 串中的资源占位符（`{threads}` / `{resources.mem_mb}` / `{resources.runtime_sec}` / `{resources.runtime_min}`）由 snakemake 运行期从各规则的 resources 声明自动填充。

snakemake 8.x 改用 executor 插件体系（`--cluster` 语义移除）；版本自动探测与适配规则见 `run.sh --help`（v0.4.0）。
