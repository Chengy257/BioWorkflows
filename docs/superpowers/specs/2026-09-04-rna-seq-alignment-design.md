# v0.4.0 对齐 rna-seq v0.8.0 工程体系 — 设计文档

- 日期：2026-09-04
- 状态：已获用户批准（含三项关键决策拍板）
- 参考标准：`D:\BioWorkflows\rna-seq`（v0.8.0 工作区状态，只读，禁止修改）
- 改造对象：本项目 chip_cuttag_atac_faire（v0.3.0 → v0.4.0）

## 1. 背景与目标

BioWorkflows 家族下的各流程子项目需要在架构、启动使用、环境软件管理等方面保持一致：
用户在任何一个子项目学到的用法可直接迁移到另一个。rna-seq 完成度最高，作为对齐标杆。

**一致性原则**：用户体验层完全一致（启动命令、配置结构、环境管理、测试入口、文档组织）；
保留 chip 合理的数据模型差异——按样本表 `seqtype` 列自动路由（chip/cuttag/atac/faire 混合项目）
优于 rna-seq 的 pipeline 选模式，属于流程逻辑差异而非工程完成度差异。

## 2. 已确认的关键决策（用户拍板）

| 决策点 | 结论 |
|---|---|
| 环境/软件管理路线 | **完全跟随 rna-seq v0.8.0**：删除 per-rule conda（envs/），改单一主环境 + `config/software.yaml` 运行时覆盖 + preflight |
| 调度器范围 | **机制完全对齐 rna-seq**：`--profile auto` 自动探测 + profile 目录结构；profile 集合 = `{default, pbs, sge, slurm}` 四套（pbs 为现实需求，sge/slurm 从 rna-seq 移植） |
| 测试深度 | **dry-run 级折中**：做合成数据生成器，CI 跑到 dry-run（不实跑工具）；实跑断言留服务器手工验证（run_test.sh 保留 `--real-run` 开关） |
| DiffBind 等死脚本 | 5 个未接入 DAG 的脚本移入 `legacy/`，README 待办保留 DiffBind 决策 |
| 版本号 | v0.4.0 |
| git 策略 | master 直接按 phase 提交，最后打 tag v0.4.0（项目惯例） |
| 工具版本 | 合并 envs 时沿用现有钉死版本（macs2=2.2.7.1、bowtie2=2.5.1、samtools=1.17、picard=3.0.0 等），不重新选版 |

**保留的 chip 反超项**（对齐不削足适履）：git tag 策略、Makefile、.editorconfig、
集中式 `validate_config`、MultiQC 注入 FRiP/NSC（`_mqc.tsv`）、`config.local.yaml` 自动叠加、
snakemake 7/8 flag 自适应（system runtime 兼容场景）、样本表 6 列 schema 与集中校验。

**规则命名**：保留 snake_case（`trim_adapter` 等）。rna-seq 的 PascalCase 是历史遗留，不学。

## 3. 目录布局迁移（对齐 rna-seq 标准布局）

| 旧（仓库根） | 新 |
|---|---|
| `workflow.smk`（303 行，含约 200 行 helper） | `workflow/Snakefile`（纯编排瘦身）；样本解析/校验/路由 helper 拆到 `workflow/rules/common.smk` |
| `rules/*.smk`（7 个模块） | `workflow/rules/*.smk`（模块划分不变：upstream/dedup/callpeak/annotation/frip/qc_deeptools/spp_qc，新增 common） |
| `scripts/*`（6 个，5 个死脚本除外） | `workflow/scripts/*` |
| `scripts/{diffpeak_DiffBind.sh, run_DiffBind.R, run_ChIPQC.R, run_chipqc_DROMPAplus.sh, annoPeak_single.R}` | `legacy/diffbind/`（归档，legacy/README.md 增加映射说明） |
| `envs/*.yaml`（11 个） | 删除（git 历史保留；版本清单并入 `workflow/environment.yaml`） |
| `profiles/pbs/` | `workflow/profile/{default,pbs,sge,slurm}/config.yaml` |
| `sample_info.example.csv` | `config/samples.csv`（模板） |
| `sample_info.csv`（真实数据） | `example/samples.csv`（随 example/ 示例项目） |
| `main_run.sh` | `run.sh`（重写，见第 5 节） |
| —（新增） | `workflow/environment.yaml`、`config/software.yaml`、`workflow/scripts/runtime_config.py`、`workflow/scripts/collect_versions.py`、`workflow/multiqc_config.yaml` |
| —（新增） | `docs/使用说明.md`、`CONTRIBUTING.md`、`example/{config.yaml, samples.csv, README.md}`、`tests/{make_testdata.py, run_test.sh, lint.sh}` |

`Makefile`、`.editorconfig`、`.gitattributes`、`.github/`、`docs/REVIEW.md`、`docs/IMPROVEMENT_PLAN.md`、
`legacy/`、`LICENSE`、`CHANGELOG.md`、`README.md`、`tests/run_tests.py` 原地保留（Makefile/CI/tests 适配新路径）。

## 4. 环境与软件管理（核心切换）

参照 rna-seq v0.8.0 的三件套 + 版本记录：

1. **`workflow/environment.yaml`**：一体化主环境模板，不自动创建，仅供
   `mamba env create -f`。钉 `snakemake-minimal=7.32.4`；合并现有 11 个 envs 的工具与版本：
   bowtie2=2.5.1、samtools=1.17、fastqc、trim-galore、picard=3.0.0、macs2=2.2.7.1、bedtools、
   deeptools、phantompeakqualtools、multiqc 等（仅 conda-forge + bioconda）。
2. **`config/software.yaml`**：`environment.type`（system / conda_prefix / conda_name）、
   R runtime（R 路径 + R_LIBS，覆盖 spp、ChIPseeker、DiffBind 等 R 包）、工具路径覆盖键。
   对齐 rna-seq 的 schema 并适配 chip 工具集。
3. **`workflow/scripts/runtime_config.py`**：自 rna-seq 移植适配（解析 / export 导出
   PATH/R_LIBS / preflight 检查三模式）。`run.sh --check-software` / `--check-r` 调用。
4. **`software_versions` 规则 + `collect_versions.py`**：运行期把实际工具版本记录到
   `results/software_versions.yaml`（对齐 rna-seq meta.smk 机制）。
5. 规则内 19 处 `conda:` 指令全部移除；`main_run.sh` 的 `-b/-e/-E` conda 部署选项移除，
   由 preflight 取代。

服务器落地（文档写清，不代做）：`mamba env create -f workflow/environment.yaml` 一次性建环境，
或 software.yaml 走 `system`/`conda_name` 复用已有环境；旧 11 个 per-rule env 用户自行清理。

## 5. 启动脚本 run.sh（完整运维 CLI，对齐 rna-seq 517 行量级）

接口能力清单（移植 rna-seq 并适配）：

- 长短选项 + 位置参数向后兼容（老 main_run.sh 调用方式可继续用）
- `--profile auto|default|pbs|sge|slurm`；auto 探测：`sbatch`→slurm；`qsub` 存在时用
  `SGE_ROOT` 环境变量区分 SGE/PBS（qsub 二义性），默认 pbs
- 集群提交串带 per-rule 资源占位符：PBS `select=1:ncpus={threads}:mem={resources.mem_mb}mb,
  walltime={resources.runtime_min}`；SGE/SLURM 对齐 rna-seq 模板（h_vmem/mem_mb、h_rt/runtime）
- `--memory`/`--runtime` 全局覆盖；`--queue`/`--partition`；`--retries`；`--latency-wait`；
  `--unlock`；`--validate-only`；`--check-software`/`--check-r`（preflight 独立模式）；
  `--skip-validation`/`--skip-software-check`
- `--log FILE` + tee 全量日志 + trap 计时/状态汇报 + INT/TERM 捕获
- `--` 后参数透传 snakemake；`--version`；`-q` 静默；`CHIP_*` 环境变量默认值机制
- 保留 chip 特有：`-r` 原始数据 R1/R2 重命名、`config.local.yaml` 自动叠加、多层 `-l` 配置、
  snakemake 7/8 探测显示（system runtime 场景的 flag 自适应保留）

## 6. per-rule 资源模型

- 21 条规则全部补 `resources.mem_mb` + `resources.runtime_min` 声明（对齐 rna-seq align.smk 风格）
- `config/config.yaml` 新增 `resources:` 覆盖段：按规则名覆盖 mem_mb/runtime_min
- 三个集群 profile 的提交串全部走占位符（见第 5 节）
- callpeak 等单线程规则保持 threads=1，mem/runtime 仍声明

## 7. 测试与 CI（dry-run 折中）

- `tests/make_testdata.py`：确定性合成数据生成器——小参考基因组（2×~100kb，含 BED/GTF/chromsize
  等注释配套）、合成 chip 型与 atac 型 reads（含接头）、6 列样本表（多组 treat/control）
- `tests/run_test.sh`：生成 → 组测试项目目录 → dry-run（CI 模式）；`--real-run` 开关保留，
  供将来服务器手工实跑；实跑后断言以内置输出文件存在性检查为准（不单独移植
  rna-seq 的 check_outputs.py 35 项断言——那是实跑折中的自然结果，服务器首跑后按需扩充）
- `tests/lint.sh`：bash -n / shellcheck / py_compile / R parse / snakemake --lint（四 profile 全过）
- `tests/run_tests.py`（45 项单测）保留并适配新路径（helper 拆分后按提取逻辑调整）
- CI 两 job：lint + dry-run（合成数据 → snakemake -n → 断言 DAG 完整：期望规则数/文件列表）；
  dry-run job 只需 pip 装 snakemake==7.32.4 + python 标准库，不需完整工具环境
- 本地验证边界（本机无 snakemake/shellcheck）：单测 + bash -n + python -S 语法检查；
  lint 全量与 dry-run 在 CI 首跑验证

## 8. 文档与清理

- `docs/使用说明.md`：对齐 rna-seq 8 章结构（安装、快速开始、样本表、配置、运行模式、
  集群提交、输出解读、FAQ）
- `CONTRIBUTING.md`：对齐 rna-seq（开发流程/CHANGELOG 要求/文档同步 checklist/测试要求/代码风格）
- `example/`：真实项目模板（现有 sample_info.csv 27 样本级数据 + 项目 config + 一条命令启动指引）
- `README.md`：重写环境/启动/目录章节对齐新体系，保留 mermaid 图、QC 阈值表、版本矩阵、
  结果路径速查等 chip 亮点；更新 DiffBind 待办说明
- `CHANGELOG.md`：v0.4.0 条目（Keep a Changelog 格式）
- git tag `v0.4.0`

## 9. 交付节奏（5 个 phase，各自独立 commit + 可验证）

1. **布局迁移**：纯移动 + 路径引用修正 + helper 拆分，行为不变；验证 = 单测全绿 + CI lint
2. **环境切换**：envs/ 删除、三件套移植（environment.yaml / software.yaml / runtime_config.py）、
   software_versions 规则、conda: 指令移除；验证 = 单测 + lint + dry-run（CI）
3. **run.sh 重写 + 资源模型**：CLI 移植、四 profile、资源占位符、21 条规则 resources、
   config resources 段；验证 = bash -n + shellcheck + CI dry-run
4. **测试与 CI**：make_testdata.py、run_test.sh、lint.sh、CI dry-run job、单测适配；
   验证 = 本地单测 + CI 两 job 全绿
5. **文档与清理**：使用说明、CONTRIBUTING、example/、README 重写、DiffBind 归档、
   CHANGELOG、tag v0.4.0；验证 = 文档交叉引用检查 + 全量 CI

## 10. 风险与缓解

| 风险 | 缓解 |
|---|---|
| 目录大迁移的路径引用遗漏（rules 内 os.path、run.sh、CI、Makefile、tests） | Phase 1 后跑单测 + lint + dry-run 三层拦截；迁移映射表写入 commit message |
| 统一环境与 PBS 服务器现实不匹配（无外网/已装环境） | software.yaml 的 system/conda_prefix 模式即为此设计；文档写清三条落地路径 |
| runtime_config.py 移植适配不全（R 包 spp/ChIPseeker、phantompeakqualtools 的 R 依赖） | preflight 独立模式可在服务器单跑检查；software_versions 规则记录实际值 |
| dry-run CI 仍发现不了 shell 运行期错误 | 已知边界，用户已确认；服务器 --real-run 为最终验证 |
| qsub 二义性导致 auto 探测错误 | SGE_ROOT 区分 + 探测不到时显式报错提示 --profile 指定 |

## 11. 验证策略

- 每个 phase 完成即 commit，commit 前：本机可跑的验证全绿（单测 + bash -n + py 语法）
- CI（lint + dry-run）为每 phase 的自动化关卡；CI 首跑问题在下一 phase 开头修复
- 最终验收 = 5 个 phase 全部合并 + CI 全绿 + tag v0.4.0；服务器实跑属后续独立任务
