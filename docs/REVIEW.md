# 项目全面审查报告（REVIEW）

> 审查日期：2026-09-03
> 审查范围：仓库内全部 24 个代码/配置文件（4 个入口 Snakefile、6 个 rules、8 个脚本、conda 环境、config、示例 samplesheet）
> 严重程度定义：
> - **P0 阻断**：必然导致流程无法运行、产出错误结果或规则永不生效
> - **P1 严重**：硬编码/环境耦合，换机器即失败；或明显的语义错误
> - **P2 一般**：代码质量、一致性、可维护性问题
> - **P3 规范**：工程化缺失（文档、测试、版本管理等）

---

## 一、项目现状概述

本项目是基于 Snakemake 的植物（水稻 IRGSP-1.0 示例）表观组学分析流程集合，覆盖 ChIP-seq / CUT&Tag / ATAC-seq / FAIRE-seq 四种数据类型，思路为：trim_galore 质控修剪 → FastQC/MultiQC → bowtie2 比对 → (picard 去重) → MACS2/SPP/HMMRATAC/FSeq2 峰调用 → bdgcmp 生成 bigWig → ChIPseeker 注释 → deeptools/DROMPAplus/ChIPQC 质控 → DiffBind 差异分析（未完成）。入口为 4 个几乎相同的 Snakefile（`chip-seq.smk` / `atac.smk` / `cuttag.smk` / `faire.smk`），配套 PBS 集群启动脚本 `main_run.sh`。

整体设计方向合理（一 套规则复用到多种 assay），但目前处于"能跑通比对、峰调用环节未接通"的状态，且大量依赖原作者机器的绝对路径。

---

## 二、P0 阻断性问题

### 2.1 峰调用规则链从未接入 DAG
- `chip-seq.smk:29`、`atac.smk:29`、`cuttag.smk:29`、`faire.smk:29` 只 `include` 了 `rules/upsteam.smk`（QC + 比对）；`rules/rmDup.smk` 与 `rules/callpeak_*.smk` 没有被任何入口 include。
- 4 个入口的 `rule all` 中，`4.peak`、注释、deeptools、motif 相关目标全部处于注释状态（如 `atac.smk:24-27`）。
- **后果**：当前任何入口都只能跑完 multiQC → bowtie2 比对，流程核心产出（peaks）不可达。

### 2.2 `rules/rmDup.smk` 整体不可用
- `rmDup.smk:11-13`：picard/samtools 的 **shell 命令被误写进 `conda:` 指令块**，`shell:` 块不存在，命令永远不会被执行；conda 指令期望的是环境 yaml 路径。
- `rmDup.smk:5`：output 使用 bash 变量语法 `${id}`，Snakemake 通配符应为 `{sample}`。
- `rmDup.smk:6-7`：`threads: {config["threads"]}`，而 config 中 threads 是字符串 `"12"`，未 `int()` 转换，会直接报错。
- `rmDup.smk:12`：路径写成 `3.align/...`，与全流程的 `3.align/bowtie2/...` 不一致。

### 2.3 `callpeak_*.smk` 系列的共性问题（4 份文件同源复制）
- output 同样使用非法的 `${id}` 语法（`callpeak_chip.smk:5` 等）。
- `callpeak_macs2` 的循环 `for j in $ID_control $ID_treat` 循环变量未使用，循环体固定 `-c control -t treat` → **每组峰调用重复执行两次**（如 `callpeak_chip.smk:32-37`）。
- `bdgcmp_macs2` 与 `peakAnno` 按 CSV 第 4 列逐行调用脚本，同一 group 的处理/对照两行会重复调用同一脚本；且隐式要求 group 列恰好写成 `treat_control` 格式。
- `peakAnno` 调用的 `Rscript /home/chengyu/workflows/annoPeak_chipseeker.R` 在仓库中不存在（实际文件是 `scripts/annoPeak_batch.R`），规则必然失败（`callpeak_chip.smk:71` 等）。
- `callpeak_atac.smk:63-78`：存在一个 **shell 为空的 `callpeak_ATAC` 规则**，其 output 与 `callpeak_macs2` 完全相同（`4.peak/all_peakfiles.txt`），一旦同时 include 会造成 DAG 冲突。

### 2.4 样本表 schema 与消费代码互相矛盾
仓库内并存两套不兼容的 CSV schema：
- `rules/callpeak_*.smk` 期望：第 1 列 id、第 4 列 group（`treat_control`）。
- `scripts/call_peak.sh:105-106` 期望：`ID_treat, ID_control, Seq_type, Peak_type`（第 2 列是对照名、第 4 列是 narrow/broad）。
- 实际 `sample_info.csv`：`ids,seqtype,layout,group` → `IgG,chip,PAIRED,myc_IgG`。
- **后果**：`call_peak.sh` 读到的 `ID_control="chip"`、`Peak_type="myc_IgG"`，case 分支必然落入报错路径。另外 `seqtype`/`layout` 列没有任何代码使用，无法自动路由到对应 assay 入口。

### 2.5 `get_samples()` 样本解析逻辑错误
`chip-seq.smk:11`（4 个入口相同）：`if line[1] not in ids: ids.append(line[0])` —— 用 **seqtype 列**判断是否加入 **样本 id**。当前 2 行样本碰巧能跑通，但属于侥幸；且不去重（同一 id 出现两次会重复追加）、不区分处理组与对照组。

### 2.6 `conda: "chip"` 指令非法
`upsteam.smk` 多处、`callpeak_*.smk` 使用 `conda: "chip"`。Snakemake 的 conda 指令只接受环境 yaml 文件路径，不接受环境名。配合 `--use-conda` 时该写法直接报错；而根目录 `chip_environment.yaml` 是整环境 export（400+ 个 pin 死的依赖、带 `prefix: /home/chengyu/...` 行、清华镜像），不是可用的 rule 级环境定义。

### 2.7 启动脚本无法直接使用
`main_run.sh:4`：`smk=` 为空变量，`-s $smk` 展开为空 → snakemake 报错，必须手改脚本。
`main_run.sh:6`：config 指向旧项目目录 `/home/chengyu/workflows/snakemake/cuttag_chip-workflow/config.yaml`。
`main_run.sh:31`：`mv ${job_title}.o* ./logs/` 假定 logs 目录存在。

---

## 三、P1 严重问题（硬编码与可移植性）

| 位置 | 硬编码内容 |
|---|---|
| `chip-seq.smk:1` | configfile 指向 `/home/chengyu/workflows/snakemake/chip_cuttag_atac_faire-workflows/config/config.yaml`（且与其余 3 个入口的相对路径写法不一致、目录名是旧的 `-workflows` 后缀） |
| `chip-seq.smk:29` 等 4 处 | include 绝对路径 `/home/chengyu/.../rules/upsteam.smk` |
| `callpeak_*.smk` bdgcmp/peakAnno | 脚本绝对路径 `/home/chengyu/workflows/snakemakem/...`、`/home/chengyu/workflows/annoPeak_chipseeker.R` |
| `upsteam.smk:80-82` | samtools 写死 `/opt/anaconda3/bin/samtools`，绕开 conda env |
| `scripts/call_peak.sh:14-15` | `/opt/anaconda3/envs/phantompeakqualtools/bin/...`（run_spp.R） |
| `scripts/call_peak.sh:77` | `/opt/biosoft/HMMRATAC_V1.2.10_exe.jar` |
| `scripts/call_peak.sh:103`、`scripts/bdgcmp_macs2.sh:11` | `/share/data/reference/osa/chrom.sizes`（参考基因组路径，config 中已有 chromsize 却不用） |
| `scripts/bdgcmp_macs2.sh:15,20` | `/share/ucsc-tools/bedClip`、`/share/ucsc-tools/bedGraphToBigWig` |
| `scripts/run_chipqc_DROMPAplus.sh:7,12` | `/share/data/reference/osa/chrom.sizes`、docker 镜像 `rnakato/ssp_drompa` 未参数化 |
| `scripts/annoPeak_batch.R:3`、`annoPeak_single.R:3` | R 包路径 `.libPaths(c("/home/chengyu/R/Rlib_4.2.3", ...))` |
| `scripts/run_ChIPQC.R:14,37,40` | `~/chipseq/results/...` 私人路径 |
| `callpeak_*.smk` params、`call_peak.sh:38,51,65` | 基因组大小 `-g 3.7e8`（水稻）写死 |
| `upsteam.smk:73` | MAPQ 过滤阈值 `-q 5` 写死在 shell 串里 |

其他脚本级 bug（P1）：
- `scripts/run_ChIPQC.R:29`：使用 `args[1]` 但从未定义 `args <- commandArgs(T)`，运行必崩；`:27` `Dir <- basename()` 缺参数；`:40` `reportName="Nanog_and_Pou5f1"` 是 ChIPQC 教程里的鼠 TF 例子名，属复制残留。
- `scripts/run_deeptools_QC.sh`：含 `{input.bed}`、`{output.matrix}`、`{log}`、`{config[threads]}` 等 Snakemake 占位符，是从 Snakefile 拷出的片段，独立执行必失败；`ids` 从 `*rmdup.bw` 提取，但流程中没有任何规则生成 `*_rmdup.bw`（bdgcmp 产出的是 `*_FE_bdgcmp.bw`，命名对不上）。

---

## 四、P2 一般问题（质量与一致性）

1. **文件重复**：根目录 `call_peak.sh` 与 `scripts/call_peak.sh` 内容完全相同；且 `scripts/call_peak.sh:143` 自带注释 `To improve: use _sorted.bam not _rmdup.bam`（CUT&Tag 去重策略 TODO 未落地）。
2. **入口文件 90% 重复**：4 个入口 smk 仅在 `rmdup.bam` 目标有无上不同，其余逐行相同，应参数化合并。
3. **命名**：`rules/upsteam.smk` 拼写错误（应为 upstream）；入口文件名 `chip-seq.smk` 与其他 `*.smk` 命名风格不统一。
4. **空壳文件**：`scripts/diffpeak_DiffBind.sh` 只有 PBS 头，`scripts/run_DiffBind.R` 完全为空 —— 差异分析有目标无实现。
5. **参数语义需确认**（不一定是错误，但需作者确认是否刻意）：
   - `callpeak_macs2` 的 `-q 0.5`（narrow peak 常规用 0.05/0.01）与 `--broad-cutoff 0.5`；
   - `--keep-dup all` 与 picard `REMOVE_DUPLICATES=true` 并存，不同 assay 的去重策略（尤其 CUT&Tag 一般不去重）没有按 assay 区分；
   - `MACS2_narrow_peak` 使用 `--nomodel --extsize ${fragment_len}`（来自 SPP 估计），而 `callpeak_macs2` 规则里 `--nomodel` 后没有 extsize，两条路径参数不一致。
6. **规则粒度**：`fastQC`/`multiQC` 是"全样本单规则"，损失并行度；`callpeak_macs2`/`bdgcmp` 把整个 group 循环塞进单条 shell，无法按组并行、失败无法断点重跑。
7. **产物命名断层**：deeptools QC 脚本期望 `*_rmdup.bw`，bigWig 实际命名为 `*_FE_bdgcmp.bw`，且生成 bigWig 的规则没有产出声明（`touch` flag 文件代替），`5.QC_deeptools` 相关规则整体缺失。
8. **config 覆盖不全**：`config.yaml` 只有 7 个键；gtf/bed/GENOME/chromsize 都是 rice 专用，没有 genome size、peak 参数、工具路径、样本表校验等。

---

## 五、P3 工程化缺失（本次已补齐部分）

| 缺失项 | 状态 |
|---|---|
| git 仓库 / 版本控制 | ✅ 本次建立（原始快照 commit `03f5c0c`） |
| README / 使用文档 | ✅ 本次补齐（见根目录 `README.md`） |
| 优化计划 | ✅ 本次补齐（见 `docs/IMPROVEMENT_PLAN.md`） |
| `.gitignore`（防分析产出误入库） | ✅ 本次补齐 |
| `.gitattributes`（强制 LF，避免 Windows CRLF 破坏 Linux 脚本） | ✅ 本次补齐（首次提交时的 CRLF 警告即风险实证） |
| `.editorconfig` / `CHANGELOG.md` | ✅ 本次补齐 |
| LICENSE | ⬜ 待作者选择（MIT / Apache-2.0 / GPL-3.0） |
| 测试数据 + 端到端小样本验证 | ⬜ 计划 Phase 3 |
| CI（snakemake --lint / shellcheck / yamllint） | ⬜ 计划 Phase 3 |
| per-rule conda 环境 | ⬜ 计划 Phase 1（拆分 `chip_environment.yaml`） |

---

## 六、结论

流程的**分析设计是合理的**（QC → 比对 → 去重策略按 assay 区分 → 峰调用按 narrow/broad 与 assay 区分 → bigWig → 注释 → QC 汇总），但重构前**只能跑通到比对一步**；峰调用链因 include 断裂、`${id}` 语法、shell/conda 块错位、CSV schema 矛盾、双重调用 bug 等问题处于不可用状态；可移植性方面几乎所有脚本都与 `/home/chengyu`、`/opt`、`/share` 路径耦合。

---

## 七、修复状态（2026-09-03，v0.2.0 重构完成）

> 重构实施见 commit `e29c428` 及后续修复 commit；v0.1.0 原文件整体归档至 `legacy/`。修复经独立代码审查（fix-first 裁决）复核，审查提出的 3 项 P1 与全部 P2/P3 均已处理。

| 问题 | 状态 | 修复方式 |
|---|---|---|
| §2.1 峰调用链未接入 DAG | ✅ 已修 | 统一入口 `workflow.smk`，rule all 全目标可达，按组并行规则 |
| §2.2 rmDup.smk 整体不可用 | ✅ 已修 | 重写为 `rules/dedup.smk`（conda/shell 块归位、`{sample}` 通配符、threads int） |
| §2.3 callpeak_*.smk 共性问题 | ✅ 已修 | 重写为 `rules/callpeak.smk`：删除空规则与双重调用，按 assay/peak_type 用互斥 wildcard_constraints 路由三规则 |
| §2.4 样本表 schema 矛盾 | ✅ 已修 | 新 schema（sample_id/role/group/seqtype/layout/peak_type）+ 入口逐行校验（含行号报错） |
| §2.5 get_samples() 解析错误 | ✅ 已修 | 重写 `load_sample_table()`（按 sample_id 去重、role 分组、混型校验），26 项单元测试覆盖 |
| §2.6 `conda: "chip"` 非法 | ✅ 已修 | 拆分为 `envs/` 下 11 个 per-rule 环境文件 |
| §2.7 main_run.sh 不可用 | ✅ 已修 | getopts 重写（-w/-s/-c/-j/-C/-p/-b/-l/-r/-n），dry-run 预检，不删 .snakemake |
| §三 硬编码路径（16 处） | ✅ 已修 | 全部经 `workflow.basedir`/config 参数化；脚本去 `.libPaths`/私人路径 |
| §四.1 文件重复 | ✅ 已修 | 根目录 `call_peak.sh` 删除，`scripts/` 为唯一真身并归档至 legacy |
| §四.2 入口 90% 重复 | ✅ 已修 | 4 入口合并为 `workflow.smk` |
| §四.3 命名不规范 | ✅ 已修 | `rules/upstream.smk`（修正拼写）等统一命名 |
| §四.5 峰参数（-q 0.5 等） | ✅ 已修 | 阈值全部入 config，默认常规值 q=0.05/broad_cutoff=0.05；CUT&Tag 不去重按 Kaya-Okur 2019 落地 |
| §四.7 产物命名断层 | ✅ 已修 | deeptools QC 规则化并使用 `{group}_FE.bw` 命名 |
| §四.8 config 覆盖不全 | ✅ 已修 | 完整 schema + `config.template.yaml` + `config.local.yaml` 叠加机制 |
| §五 LICENSE 缺失 | ✅ 已修 | MIT |
| §五 测试/CI | ✅ 已修 | `tests/run_tests.py`（26 项）+ Makefile + GitHub Actions CI |

**审查修复补充项**（独立审查 fix-first 发现）：r-chipseeker 环境 R 版本冲突（升级 r-base=4.3 + Bioc 3.18 对齐）；ATAC shift/extsize 在 BAMPE 下被静默忽略（改为 `peak.atac.mode` 双模式开关，默认 ENCODE ATAC v2 的 bampe）；bigwig 字典序排序与 chrom.sizes 顺序冲突（改 `bedtools sort -g`，chromsize 提为 input）；main_run.sh `-j/--cores` 同参数覆盖（按集群/本机模式拆分）；样本/分组名非法字符校验；envs 移除 defaults 频道。

**遗留（不阻塞交付）**：服务器真实数据端到端实跑（含 MACS2 无对照 control_lambda.bdg 产出确认、conda 环境求解、CI 首跑）；DiffBind 差异分析（需 contrast 设计决策）；bowtie2 `.bt2l` 大基因组限制（已在 README 声明）。
