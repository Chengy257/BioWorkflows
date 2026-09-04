# workflow/profile/

Snakemake run profiles (cluster parameters pinned into the repository). The
scheduling mode is selected via `--profile workflow/profile/<name>`; for the
auto-detection rules see `run.sh --help` (v0.4.0).

| profile | scheduler | cluster command highlights |
|---|---|---|
| `default/` | local / standalone server | no cluster command, runs directly with `--cores 8`; `latency-wait: 90` |
| `pbs/` | PBS (Torque), snakemake **7.x** classic interface | `qsub -V -N {rule} -l select=1:ncpus={threads}:mem={resources.mem_mb}mb -l walltime={resources.runtime_sec} -j oe`; walltime in seconds avoids `[HH:]MM:]SS` ambiguity |
| `sge/` | SGE, snakemake **7.x** classic interface | `qsub -V -N {rule} -l ncpus={threads} -l h_vmem={resources.mem_mb}M -l h_rt={resources.runtime_sec} -j oe` |
| `slurm/` | SLURM, snakemake **7.x** classic interface | `sbatch --parsable -J {rule} -c {threads} --mem={resources.mem_mb}M --time={resources.runtime_min} -o slurm-{rule}-%j.out` |

Resource placeholders in the cluster command (`{threads}` /
`{resources.mem_mb}` / `{resources.runtime_sec}` / `{resources.runtime_min}`)
are filled at runtime by snakemake from each rule's resources declaration
(centralized in config/resources.yaml).

snakemake 8.x moved to the executor plugin system (`--cluster` semantics
removed); for version auto-detection and adaptation rules see `run.sh --help`
(v0.4.0).
