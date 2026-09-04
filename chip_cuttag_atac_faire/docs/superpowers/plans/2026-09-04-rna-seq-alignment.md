# v0.4.0 对齐 rna-seq v0.8.0 工程体系 — 实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 把 chip_cuttag_atac_faire 的目录布局、环境管理、启动脚本、资源模型、测试 CI、文档全面对齐 rna-seq v0.8.0，交付 v0.4.0。

**Architecture:** 五阶段线性改造：布局迁移（行为不变）→ 环境路线切换（per-rule conda → 统一环境 + software.yaml + preflight）→ run.sh 重写 + per-rule 资源模型 → 测试/CI 对齐（dry-run 级）→ 文档与清理。每阶段独立 commit、独立可验证。

**Tech Stack:** Snakemake 7.32.4（钉版）、bash、Python 3（PyYAML）、conda/mamba（仅 environment.yaml 手工创建）、GitHub Actions。

**设计文档:** `docs/superpowers/specs/2026-09-04-rna-seq-alignment-design.md`（所有决策依据）

## Global Constraints

- **绝对禁止修改 `D:\BioWorkflows\rna-seq` 的任何文件**——它是只读移植源。移植 = 读取原件 + 在本仓库创建适配版本。
- 每个 task 结束必须 `make check PYTHON=python` 全绿（本机无 snakemake/shellcheck；本机全局 Python 损坏，本地跑单测用 venv 中的 python，CI 用 `python`）。
- commit 信息风格沿用现有历史：`feat(scope): ...` / `fix: ...` / `docs: ...` / `refactor: ...`。
- 不引入新依赖到 tests（run_tests.py 保持零依赖自研 check 框架）。
- 环境变量前缀统一 `CHIP_`；脚本版本号对齐目标 `0.4.0`。
- 规则命名保持 snake_case 不变；rules 模块文件名与划分不变（仅新增 common/meta 两块）。

---

## Phase 1 — 目录布局迁移（行为不变）

**Phase 验收：** 所有文件迁入 `workflow/` 标准布局；`make check` 全绿；`git mv` 保留历史；`--list-rules` 在 CI 可解析。

### Task 1.1: 文件树迁移 + 入口改写

**Files:**
- Create: `workflow/Snakefile`（由 `workflow.smk` `git mv` 而来后修改）
- Create: `workflow/rules/common.smk`（新建，承接 helper）
- Modify: `workflow.smk` → `workflow/Snakefile`
- Delete: 根目录 `workflow.smk`（即 mv 结果）
- Move: `rules/*.smk` → `workflow/rules/*.smk`（git mv）
- Move: `scripts/annoPeak_batch.R` → `workflow/scripts/annoPeak_batch.R`（git mv；其余 5 个死脚本 Phase 5 归档，本 task 不动）
- Modify: `Makefile`、`.github/workflows/ci.yaml`（路径引用）

**Interfaces:**
- Produces: `workflow/Snakefile` 中定义的全局名（供各 rules 与测试引用）：`BASE_DIR`（仓库根）、`WORKFLOW_DIR`（workflow/ 目录）、`SAMPLES, GROUPS, SEQTYPE_OF, load_sample_table, validate_config, assay_needs_dedup, sample_bam, group_bams, group_control_arg, group_peak_file, _groups_of, _group_regex, _resolve_sample_table, _NAME_RE, REQUIRED_COLUMNS, ASSAYS`——全部迁入 `workflow/rules/common.smk`，Snakefile 通过 `include` 获得。

- [ ] **Step 1: git mv 移动文件**

```bash
mkdir -p workflow
git mv workflow.smk workflow/Snakefile
git mv rules workflow/rules
git mv scripts/annoPeak_batch.R workflow/scripts/annoPeak_batch.R   # 先 mkdir -p workflow/scripts
git status   # 确认全部为 renamed
```

- [ ] **Step 2: 改写 workflow/Snakefile 为纯编排（照 rna-seq Snakefile 模式）**

保留的头部（替换原 `configfile`/`REPO_DIR` 行）：

```python
import os

BASE_DIR = os.path.dirname(workflow.basedir)   # Repository root.
WORKFLOW_DIR = workflow.basedir                # workflow/ directory.

configfile: os.environ.get("CHIP_CONFIG", os.path.join(BASE_DIR, "config", "config.yaml"))
```

`CHIP_CONFIG` 环境变量机制为 run.sh 预留（Task 3.3 使用）。紧随其后保留原 `shell.executable`/`shell.prefix` 两行，然后：

```python
# 样本表解析、config 校验、常用查询函数、目标汇总常量的全部实现
# 迁移至 rules/common.smk（本文件按出现顺序使用它们）。
include: os.path.join(WORKFLOW_DIR, "rules", "common.smk")
```

随后保留原 `rule all:` 与末尾 7 个 `include:`（路径改 `os.path.join(WORKFLOW_DIR, "rules", ...)`）。原 293-302 行的 include 路径 `REPO_DIR, "rules"` 全部改为 `WORKFLOW_DIR, "rules"`。

- [ ] **Step 3: 新建 workflow/rules/common.smk**

把 `workflow/Snakefile` 原第 13-278 行的内容（`import csv/os/re`、`from snakemake.exceptions import WorkflowError`、`shell.executable(...)`、`shell.prefix(...)`、常量 `_NAME_RE/ASSAYS/REQUIRED_COLUMNS`、`_resolve_sample_table/load_sample_table`、threads/mode 白名单校验、`validate_config`、`assay_needs_dedup` 到 `_group_regex`、`BAM_TARGETS/QC_TARGETS/PEAK_TARGETS/BW_TARGETS`）整体搬入 `common.smk`，文件头注释：

```python
# ---------------------------------------------------------------------
# 共享定义：路径常量、样本表解析、config 校验、规则查询函数、目标汇总
# 由 Snakefile 第一个 include；各 rules/*.smk 直接使用这里定义的名字。
# ---------------------------------------------------------------------
```

关键适配（逐条）：
- `ENVS = os.path.join(REPO_DIR, "envs")` 一行**删除**（Phase 2 前规则仍引用 envs 的相对写法？否——Phase 1 只动布局：rules 里的 `conda:` 指令原样保留，它们引用 `ENVS` 全局变量。所以 Phase 1 **保留** `ENVS = os.path.join(BASE_DIR, "envs")`，Phase 2 Task 2.4 再删）。注意 `REPO_DIR` 全局名在 common.smk 中改为 `BASE_DIR`（与 rna-seq 命名对齐），规则文件中所有 `REPO_DIR` 引用随本 task 全部替换为 `BASE_DIR`（`grep -rn REPO_DIR workflow/` 归零）。
- `_resolve_sample_table` 内 `os.path.join(REPO_DIR, p)` → `os.path.join(BASE_DIR, p)`；`validate_config` 内同。
- configfile 加载移到 Snakefile（common.smk 不含 `configfile:` 指令，但含读取 `config[...]` 的语句——`SAMPLES, GROUPS, SEQTYPE_OF = load_sample_table(...)` 保持在 common.smk，因为它在 include 时求值，Snakefile 的 configfile 已先执行）。

- [ ] **Step 4: 修正 workflow/rules/*.smk 内路径引用**

```bash
grep -rn "REPO_DIR\|scripts/\|\"envs\"" workflow/rules/   # 逐个替换
```
- `REPO_DIR` → `BASE_DIR`（约 19 处 conda 引用 + scripts 引用）
- `os.path.join(REPO_DIR, "scripts", ...)`（annotation.smk 的 R 脚本引用）→ `os.path.join(WORKFLOW_DIR, "scripts", ...)`
- `os.path.join(ENVS, ...)` 保持（Phase 2 删）

- [ ] **Step 5: 适配 Makefile 与 CI 路径**

Makefile：`-s workflow.smk` → `-s workflow/Snakefile`；`SHELL_SCRIPTS := main_run.sh $(wildcard workflow/scripts/*.sh)`。ci.yaml：`snakemake -s workflow.smk` → `-s workflow/Snakefile`，shellcheck 行 `scripts/*.sh` → `workflow/scripts/*.sh`。

- [ ] **Step 6: 适配 tests/run_tests.py 的提取源**

run_tests.py 三处从 `workflow.smk` 文本提取（`extract_workflow_functions` 读 `os.path.join(REPO, "workflow.smk")`；`mode_block` 提取；`vc_block` 提取）。改为读 `workflow/rules/common.smk`：

```python
COMMON = os.path.join(REPO, "workflow", "rules", "common.smk")
```
提取标记不变（函数 def 行照样存在），`segment()`/`vc_block` 的起止锚点核对一遍（`BAM_TARGETS` 锚点仍在 common.smk）。`render_rule_body` 的 smk 路径参数 `rules/xxx.smk` → `workflow/rules/xxx.smk`（两处调用点）。

- [ ] **Step 7: 验证 + commit**

```bash
make check PYTHON=<venv-python>     # 45 项全过
bash -n main_run.sh
git add -A && git commit -m "refactor: migrate to standard workflow/ layout (Snakefile + rules/common.smk split)"
```

### Task 1.2: 样本表模板迁移 + profile 目录骨架

**Files:**
- Move: `sample_info.example.csv` → `config/samples.csv`（git mv）
- Move: `profiles/pbs/` → `workflow/profile/pbs/`（git mv，仅迁移位置；内容 Phase 3 重写）
- Create: `workflow/profile/README.md`（从 `profiles/README.md` git mv）
- Modify: `config/config.yaml`（grouplist 注释指向新模板名）、`README.md` 中 sample_info.example.csv 引用

**Interfaces:**
- Produces: `config/samples.csv` 为官方模板（列头不变：`sample_id,role,group,seqtype,layout,peak_type`）；根目录 `sample_info.csv`（真实数据）保持原地不动，Phase 5 移入 example/。

- [ ] **Step 1: 移动与引用修正**

```bash
git mv sample_info.example.csv config/samples.csv
git mv profiles workflow/profile
```
`config/config.yaml` 的 `grouplist` 注释改为：`# 样本表；相对路径相对工作目录解析；模板见 config/samples.csv`。`README.md` 里所有 `sample_info.example.csv` → `config/samples.csv`。`grep -rn "profiles/" Makefile .github/ tests/` 核对无旧引用。

- [ ] **Step 2: 验证 + commit**

```bash
make check PYTHON=<venv-python>
git add -A && git commit -m "refactor: move sample template to config/samples.csv, profiles into workflow/profile"
```

---

## Phase 2 — 环境路线切换（per-rule conda → 统一环境 + software.yaml）

**Phase 验收：** `envs/` 删除；规则无 `conda:` 指令；三件套（environment.yaml / software.yaml / runtime_config.py）可用；`software_versions` 规则进 DAG；CI 全绿。

### Task 2.1: workflow/environment.yaml（统一主环境）

**Files:**
- Create: `workflow/environment.yaml`

- [ ] **Step 1: 写入完整内容（版本 = 现有 11 个 envs 的钉版合并，不新选版）**

```yaml
# 全新服务器的一体化环境模板。
# 此文件不会被 Snakemake 自动创建；示例：mamba env create -f workflow/environment.yaml
name: chip-cuttag-atac-faire
channels:
  - conda-forge
  - bioconda
  - nodefaults
dependencies:
  - python=3.11
  - git
  - pyyaml
  - snakemake-minimal=7.32.4
  # 上游 QC 与修剪
  - fastqc=0.12.1
  - trim-galore=0.6.10
  - multiqc=1.14
  # 比对
  - bowtie2=2.5.1
  - samtools=1.17
  # 去重 / 峰调用 / 信号
  - picard=3.0.0
  - macs2=2.2.7.1
  - bedtools=2.31.0
  - deeptools=3.5.1
  - ucsc-bedclip=482
  - ucsc-bedgraphtobigwig=482
  # SPP NSC/RSC
  - phantompeakqualtools=1.2.2
  - r-base=4.3
  - bioconductor-chipseeker=1.38.0
  - bioconductor-genomicfeatures=1.54.0

# R 包注：ChIPseeker/GenomicFeatures 已含 conda 版；若服务器复用已有 R 库，
# 走 config/software.yaml 的 r.lib_paths 指定，不必安装本段 bioconductor-*。
```

- [ ] **Step 2: 验证 + commit**

```bash
python -c "import yaml;yaml.safe_load(open('workflow/environment.yaml'))"   # 用 venv python
git add workflow/environment.yaml && git commit -m "feat(env): unified all-in-one environment template (merged pinned versions from 11 per-rule envs)"
```

### Task 2.2: config/software.yaml + workflow/scripts/runtime_config.py

**Files:**
- Create: `config/software.yaml`
- Create: `workflow/scripts/runtime_config.py`（移植源：`D:\BioWorkflows\rna-seq\workflow\scripts\runtime_config.py`，291 行，**只读**）

**Interfaces:**
- Produces: CLI 子命令 `export --config FILE`（stdout 输出可 eval 的 `export KEY=VALUE` 行）与 `check --config FILE --scope {all,software,r}`（preflight，非零退出码 + 人读错误清单）；输出环境变量名 `CHIP_SOFTWARE_TYPE / CHIP_ENV_PREFIX / CHIP_RSCRIPT / R_LIBS / CHIP_TOOL_*`。

- [ ] **Step 1: 写 config/software.yaml（对齐 rna-seq schema，工具集换成 chip 的）**

```yaml
# 统一软件/运行时配置。
# 复制到项目目录为 software.yaml 并仅覆盖差异项。

environment:
  type: system                  # system | conda
  conda_prefix: ""              # HPC 推荐，如 /share/conda/envs/chip-cuttag-atac-faire
  conda_name: ""                # conda_prefix 的替代
  strict: true

r:
  rscript: "Rscript"            # 命令名或绝对路径
  version: ""                   # 可选，如 "4.3"
  version_check: major_minor    # major_minor | exact | warn | off
  lib_paths: []                 # 可选 R 库目录（复用服务器已有 ChIPseeker 等时填）
  lib_mode: prepend             # prepend | append | replace
  package_sources: {}           # 可选 包名 -> 本地 tarball 路径（仅用于安装提示）

# 常规工具自动从主环境/PATH 解析；仅覆盖差异项。
tools:
  trim_galore: "trim_galore"
  bowtie2: "bowtie2"
  samtools: "samtools"
  picard: "picard"              # picard.jar 完整路径亦可
  macs2: "macs2"
  bedtools: "bedtools"
  deeptools: "deeptools"
  spp: "spp.R"                  # phantompeakqualtools 的 R 入口

paths: {}

databases: {}
```

- [ ] **Step 2: 移植 runtime_config.py**

读取 rna-seq 原件，逐段适配（其余逻辑保持原样，含错误信息双语风格对齐）：
1. 全部 `RNASEQ_` 前缀 → `CHIP_`（环境变量导出名与读取处）。
2. `PIPELINE_TOOLS` 分 pipeline 的工具映射表 → 单一 `WORKFLOW_TOOLS`（chip 无 pipeline 概念）：`["python", "snakemake", "rscript", "trim_galore", "bowtie2", "bowtie2-build", "samtools", "picard", "macs2", "bedtools", "deeptools", "spp", "multiqc", "fastqc"]`。
3. argparse：删 `--pipeline` 参数；`check` 保留 `--config`/`--scope`/`--analysis-config`（analysis-config 可留作扩展，默认不启用——chip 校验集中在 Snakefile）。
4. R 检查段（r_eval/包检查）保留，R 包清单改为 chip 的：`GenomicFeatures, ChIPseeker`（spp 检查归入工具级 `spp.R --version`，不做包级）。
5. `resolve_runtime`/`shell_exports`/`check_runtime` 函数签名与返回结构不变（run.sh 依赖）。

- [ ] **Step 3: 本地验证（venv python，无外部工具）**

```bash
<venv-python> -m py_compile workflow/scripts/runtime_config.py
printf 'environment:\n  type: system\n  strict: false\nr: {rscript: Rscript, version_check: "off"}\ntools: {}\npaths: {}\ndatabases: {}\n' > /tmp/soft.yaml
<venv-python> workflow/scripts/runtime_config.py export --config /tmp/soft.yaml   # 输出 export CHIP_... 行
<venv-python> workflow/scripts/runtime_config.py check --config /tmp/soft.yaml --scope software; echo "exit=$?"   # 本机缺工具 → 非零退出 + 清单（预期行为）
```

- [ ] **Step 4: commit**

```bash
git add config/software.yaml workflow/scripts/runtime_config.py
git commit -m "feat(env): software.yaml runtime config + runtime_config.py (export/preflight), ported from rna-seq v0.8.0"
```

### Task 2.3: software_versions 规则 + collect_versions.py

**Files:**
- Create: `workflow/rules/meta.smk`
- Create: `workflow/scripts/collect_versions.py`（移植源：`D:\BioWorkflows\rna-seq\workflow\scripts\collect_versions.py`，只读）
- Modify: `workflow/Snakefile`（rule all 目标 + include meta.smk；common.smk 的 QC/BW 目标汇总处加版本文件）

**Interfaces:**
- Consumes: Task 2.2 的 `CHIP_TOOL_*` 导出约定；common.smk 的 `config`。
- Produces: 规则 `software_versions`，输出 `5.QC/software_versions.yaml`（chip 无 rna-seq 的 `R()` 结果目录函数——直接写相对路径常量，见下）。

- [ ] **Step 1: 移植 collect_versions.py**

照 rna-seq 原件（结构：command_output/git_commit/resolved_tool + argparse 主流程），适配：
- `TOOLS` 清单替换为 Task 2.2 Step 2 的 `WORKFLOW_TOOLS`；
- 环境变量前缀 `RNASEQ_TOOL_` → `CHIP_TOOL_`；
- defaults 映射保留 `{"python": "python3", "rscript": "Rscript"}`，追加 `{"spp": "spp.R"}`。

- [ ] **Step 2: 写 workflow/rules/meta.smk**

```python
# 运行元数据：实际软件版本记录。
rule software_versions:
    output:
        "5.QC/software_versions.yaml",
    log:
        "5.QC/logs/software_versions.log.txt",
    params:
        script=os.path.join(WORKFLOW_DIR, "scripts", "collect_versions.py"),
        software_config=os.environ.get("CHIP_SOFTWARE_CONFIG",
                                       os.path.join(BASE_DIR, "config", "software.yaml")),
        workflow_dir=WORKFLOW_DIR,
    threads: 1
    resources:
        mem_mb=1024,
        runtime_min=10,
        runtime_sec=600,
    shell:
        "python3 {params.script} --software-config {params.software_config} "
        "--workflow {params.workflow_dir} --out {output} >> {log} 2>&1"
```
（snakemake 版本参数先不带——rna-seq 用 importlib.metadata 取，chip 版 collect_versions.py 内部自行 `importlib.metadata.version("snakemake")`，更简单。）

- [ ] **Step 3: 接线**

common.smk 的目标汇总追加 `VERSION_TARGETS = ["5.QC/software_versions.yaml"]`；Snakefile `rule all` 的 input 加 `VERSION_TARGETS`；末尾 `include: os.path.join(WORKFLOW_DIR, "rules", "meta.smk")`。

- [ ] **Step 4: 验证 + commit**

```bash
<venv-python> -m py_compile workflow/scripts/collect_versions.py
make check PYTHON=<venv-python>   # --list-rules 本机不可用，CI 验证
git add -A && git commit -m "feat(meta): software_versions rule + collect_versions.py capture runtime provenance"
```

### Task 2.4: 删除 per-rule conda 体系

**Files:**
- Delete: `envs/`（git rm -r，11 个 yaml）
- Modify: `workflow/rules/*.smk`（删全部 `conda:` 指令块，约 19 处，位置见 grep 清单：annotation.smk:15、callpeak.smk:28/61/102/129、dedup.smk:14、upstream.smk:19/42/67/84/107、frip.smk:21、spp_qc.smk:18、qc_deeptools.smk:17/36/54/73/94/115）
- Modify: `workflow/rules/common.smk`（删 `ENVS` 常量）
- Modify: `Makefile`（dryrun 目标删 `--use-conda`）
- Modify: `main_run.sh`（最小修补，Phase 3 重写前保持可用：删 `-b/-e/-E` 选项与 `--use-conda/--software-deployment-method conda/--conda-*` 拼接段，即 47-65 行的 b/e/E case、100-120 的版本分支、133-142 的 conda 参数段）
- Modify: `README.md` 环境章节加临时说明（Phase 5 重写）：环境按 `workflow/environment.yaml` 统一创建

- [ ] **Step 1: 执行删除与修补**

```bash
git rm -r envs/
grep -rn "conda" workflow/rules/ Makefile main_run.sh   # 逐处清理至无 per-rule conda 引用
```

- [ ] **Step 2: 验证 + commit**

```bash
make check PYTHON=<venv-python>
bash -n main_run.sh
git add -A && git commit -m "feat(env)!: drop per-rule conda envs in favor of unified environment + software.yaml runtime"
```
（`!` 表示 breaking change，CHANGELOG Phase 5 记录。）

---

## Phase 3 — run.sh 重写 + per-rule 资源模型

**Phase 验收：** `run.sh` 具备 rna-seq 全部运维能力 + pbs；四 profile 就位；21+1 条规则有 mem/runtime 声明；config `resources:` 覆盖段生效；`main_run.sh` 删除；CI dry-run 绿。

### Task 3.1: 四套 profile

**Files:**
- Create: `workflow/profile/default/config.yaml`、`workflow/profile/sge/config.yaml`、`workflow/profile/slurm/config.yaml`
- Modify: `workflow/profile/pbs/config.yaml`（重写内容）

- [ ] **Step 1: 写入四份 profile（default 对齐 rna-seq；sge/slurm 抄 rna-seq 原文；pbs 新写）**

default（同 rna-seq，仅 latency-wait 90 保留 chip 对 PBS 共享盘的经验值）：

```yaml
# 本机/独立服务器执行 profile。
cores: 8
keep-going: true
rerun-incomplete: true
latency-wait: 90
printshellcmds: true
```

pbs（资源占位符与 SGE/SLURM 同源，walltime 用秒避免 [[HH:]MM:]SS 歧义）：

```yaml
# PBS (Torque) 执行 profile，Snakemake >= 7 经典集群接口。
# 资源占位符由每条规则的 resources 声明填充。
cluster: "qsub -V -N {rule} -l select=1:ncpus={threads}:mem={resources.mem_mb}mb -l walltime={resources.runtime_sec} -j oe"
jobs: 20
keep-going: true
rerun-incomplete: true
latency-wait: 90
printshellcmds: true
```

sge / slurm：逐字复制 `D:\BioWorkflows\rna-seq\workflow\profile\sge\config.yaml` 与 `...\slurm\config.yaml`（只读源；内容见设计文档引用，两行 cluster 串不变）。

- [ ] **Step 2: 更新 workflow/profile/README.md**：四 profile 一览表 + auto 探测规则说明（指向 run.sh --help）。

- [ ] **Step 3: 验证 + commit**

```bash
python -c "import yaml,glob;[yaml.safe_load(open(f)) for f in glob.glob('workflow/profile/*/config.yaml')]"   # venv python
git add -A && git commit -m "feat(profile): default/pbs/sge/slurm profiles with per-rule resource placeholders"
```

### Task 3.2: 21+1 条规则补 resources + config 覆盖段

**Files:**
- Modify: `workflow/rules/*.smk` 全部规则 + `workflow/rules/meta.smk`（已带）
- Modify: `workflow/rules/common.smk`（加资源 helper）
- Modify: `config/config.yaml`（加 `resources:` 段）
- Modify: `tests/run_tests.py`（加资源声明覆盖单测）

**Interfaces:**
- Produces（common.smk）:

```python
def res(rule, key, default):
    """按规则名读取 config resources 覆盖段，缺省回落到规则默认值。"""
    override = (config.get("resources") or {}).get(rule) or {}
    return override.get(key, default)
```

- [ ] **Step 1: 每条规则在 threads 后追加 resources**

统一模板（数值表如下；threads 保持现状不动，无 threads 的 4 条规则补 `threads: 1`：`frip_summary`、`spp_summary`、`deeptools_correlation`、`deeptools_pca`）：

```python
    resources:
        mem_mb=res("<rule>", <MEM>),
        runtime_min=res("<rule>", <RUNTIME_MIN>),
        runtime_sec=res("<rule>", <RUNTIME_MIN>) * 60,
```

数值表（mem_mb / runtime_min）：

| 规则 | mem_mb | runtime_min |
|---|---|---|
| trim_adapter | 4096 | 60 |
| fastqc | 2048 | 30 |
| multiqc | 4096 | 30 |
| bowtie2_index | 8192 | 120 |
| bowtie2_mapping | 16384 | 240 |
| dedup | 8192 | 120 |
| callpeak_narrow | 8192 | 180 |
| callpeak_broad | 8192 | 180 |
| callpeak_atac | 8192 | 180 |
| bigwig | 4096 | 60 |
| peak_annotation | 8192 | 120 |
| frip | 4096 | 60 |
| frip_summary | 1024 | 10 |
| spp_crosscorr | 8192 | 180 |
| spp_summary | 1024 | 10 |
| deeptools_multibamsummary | 8192 | 120 |
| deeptools_correlation | 4096 | 30 |
| deeptools_pca | 4096 | 30 |
| deeptools_fingerprint | 8192 | 60 |
| deeptools_fragmentsize | 8192 | 60 |
| deeptools_profile | 8192 | 120 |
| software_versions | 1024 | 10 |

注意：`res` 返回 runtime_sec 时若覆盖段只给 runtime_min，则 `* 60` 联动正确；覆盖段直接给 runtime_sec 不支持（YAGNI，文档注明）。

- [ ] **Step 2: config/config.yaml 追加覆盖段（示例注释即文档）**

```yaml
# ---------- 集群资源覆盖（可选；按规则名覆盖，未列出的用规则默认值） ----------
# 只支持 mem_mb / runtime_min 两个键；runtime_sec 由 runtime_min 自动换算。
resources:
#  bowtie2_mapping:
#    mem_mb: 32768
#    runtime_min: 480
  {}
```
（YAML 空段写 `resources: {}` 或全部注释掉 `resources:` + 注释示例——采用后者时 Snakefile 端 `config.get("resources")` 已兼容 None。）

- [ ] **Step 3: 单测覆盖**

run_tests.py 新增 check 组（从 common.smk 提取 `res` 后 exec 测）：

```python
ns = {"config": {}}   # 提取 res 的 exec 命名空间，config 可注入
check("res 无覆盖段回落默认值", ns["res"]("x", "mem_mb", 4096) == 4096)
ns2 = {"config": {"resources": {"bowtie2_mapping": {"mem_mb": 32768}}}}
check("res 覆盖 mem_mb 生效", ns2["res"]("bowtie2_mapping", "mem_mb", 16384) == 32768)
check("res 覆盖段缺 key 回落", ns2["res"]("bowtie2_mapping", "runtime_min", 240) == 240)
check("res 其他规则不受影响", ns2["res"]("frip", "mem_mb", 4096) == 4096)
```
（提取方式复用现有 `extract_workflow_functions` 的 segment 机制，锚点 `def res(`。）
另加静态扫描 check：`workflow/rules/*.smk` 中每个 `rule ` 块含 `runtime_sec`（渲染级断言，防漏声明）。

- [ ] **Step 4: 验证 + commit**

```bash
make check PYTHON=<venv-python>
git add -A && git commit -m "feat(resources): per-rule mem_mb/runtime declarations + config resources override"
```

### Task 3.3: run.sh 完整重写（替换 main_run.sh）

**Files:**
- Create: `run.sh`（底稿：`D:\BioWorkflows\rna-seq\run.sh` 517 行，只读）
- Delete: `main_run.sh`（git rm）
- Modify: `Makefile`（SHELL_SCRIPTS 与文档引用）、`.github/workflows/ci.yaml`（shellcheck 目标）、`README.md` 快速开始临时改指 run.sh（Phase 5 重写）

**Interfaces:**
- Consumes: Task 2.2 runtime_config.py（`export`/`check`）、Task 3.1 四 profile、`CHIP_CONFIG` 环境变量（Task 1.1）。
- Produces: `bash run.sh [选项] [工作目录]`，向后兼容旧位置参数形态 `bash main_run.sh -w DIR ...`。

- [ ] **Step 1: 以 rna-seq run.sh 为底稿逐节适配**

结构保留（usage/选项解析/路径解析/runtime 解析/select_profile/校验/CLUSTER_CMD/执行/日志），适配清单：

1. **头部常量**：`SCRIPT_VERSION="0.4.0"`；`SNAKEFILE="$WORKFLOW_DIR/Snakefile"`；`DEFAULT_SOFTWARE=.../config/software.yaml`；`RUNTIME_HELPER="$WORKFLOW_DIR/scripts/runtime_config.py"`。
2. **pipeline 概念删除**：删 `-p/--pipeline`、PIPELINE case 校验、`RNASEQ_PIPELINE` 导出；usage 的 Pipelines 节换成 Assays 说明（"assay 由样本表 seqtype 列自动路由，无需选择"）。
3. **位置参数向后兼容**（对齐旧 main_run.sh 习惯而非 rna-seq 形态）：`PROJECT_DIR="${POSITIONAL[0]:-}"`（工作目录必填），其余位置参数不接（-c/-j 由选项给）。
4. **环境变量前缀**：`RNASEQ_*` → `CHIP_*`（JOBS/LOG/QUEUE/PARTITION/MEMORY/RUNTIME_MIN/RETRIES/LATENCY_WAIT/SCHEDULER_EXTRA 等，默认值 JOBS=3、LOG=snakemake.logs.txt、LATENCY_WAIT 默认空走 profile 的 90）。
5. **config 解析**：`export CHIP_CONFIG` + `--configfile` 组合——chip 的 configfile 由 Snakefile 的 `CHIP_CONFIG` 读取（Task 1.1 已接），run.sh 不再单独传 `--configfile`（保留 config.local.yaml 叠加逻辑：`CHIP_CONFIG` 设为多值时 Snakefile 端 `os.environ.get` 只取一个——**修正方案**：run.sh 把最终 config 链通过 `--configfile a b` 传给 snakemake 命令行（snakemake 支持多 --configfile 依序覆盖），`CHIP_CONFIG` 仅传第一个做解析用。在 run.sh 注释里写明）。
6. **select_profile 增加 pbs**：

```bash
    if [[ "$requested" == "auto" ]]; then
        if command -v sbatch >/dev/null 2>&1; then
            requested="slurm"
        elif command -v qsub >/dev/null 2>&1; then
            # qsub 二义性：SGE 环境必有 SGE_ROOT；否则按 PBS 处理
            if [[ -n "${SGE_ROOT:-}" ]]; then requested="sge"; else requested="pbs"; fi
        else
            requested="default"
        fi
    fi
    case "$requested" in
        default|pbs|sge|slurm) ;;
        *) die "Unknown profile '$PROFILE_REQUEST'. Choose auto, default, pbs, sge, or slurm." ;;
    esac
```
（探测顺序 sbatch 优先——SLURM 集群一般同时无 qsub；qsub 场景用 SGE_ROOT 区分。）pbs 的 command 校验：`command -v qsub`。

7. **CLUSTER_CMD 增加 pbs case**（覆盖 --memory/--runtime 时拼串，否则依赖 profile yaml 内的串——与 rna-seq 相反的方向：rna-seq 总是拼。为一致性照 rna-seq 总是拼）：

```bash
    pbs)
        if [[ -n "$MEMORY" ]]; then
            PBS_MEM_REQUEST="$MEMORY"
        else
            PBS_MEM_REQUEST='{resources.mem_mb}mb'
        fi
        if [[ -n "$RUNTIME_MIN" ]]; then
            PBS_RUNTIME_REQUEST="$((RUNTIME_MIN * 60))"
        else
            PBS_RUNTIME_REQUEST='{resources.runtime_sec}'
        fi
        CLUSTER_CMD="qsub -V -N {rule} -l select=1:ncpus={threads}:mem=${PBS_MEM_REQUEST} -l walltime=${PBS_RUNTIME_REQUEST} -j oe"
        [[ -n "$QUEUE" ]] && CLUSTER_CMD+=" -q $QUEUE"
        ;;
```
（--queue 允许 pbs 与 sge；`--partition` 仅 slurm；--sge-mem-resource 仅 sge——校验分支相应调整。）

8. **runtime_check 调用**：`"$RUNTIME_PYTHON" "$RUNTIME_HELPER" check --config "$SOFTWARE_PATH" --scope "$scope"`（无 --pipeline）。
9. **样本表校验**：chip 无 validate_samples.py——对应段删除，替换为提示"样本表与 config 校验在 Snakefile 解析期集中执行"；`--validate-only` 实现为 `snakemake --list-rules` + dry-run 前退出（解析期校验自然触发）。
10. **保留 chip 特有能力**（从 main_run.sh 移植，插在 cd "$PROJECT_DIR" 后）：
    - `-r/--rename`：原 main_run.sh 79-92 行的重命名块原样保留；
    - config.local.yaml 自动叠加：原 94-98 行逻辑 + `-l/--extra-config FILE` 选项；
    - snakemake 7/8 版本探测仅用于 info 显示与告警（不再拼 conda flag）。
11. **usage/examples 全面改写为 chip 语境**：

```
Examples:
  bash run.sh /data/project -j 3 --profile pbs
  bash run.sh -P /data/project -c config.yaml --check-software
  bash run.sh -P . -n
  bash run.sh -P . --profile pbs --queue workq --memory 16G --runtime 600
  bash run.sh -P . --unlock
  bash run.sh -P . -- --rerun-triggers mtime
```
12. **集群输出回收**：PBS `mv ./[a-zA-Z]*.o* ./logs/ 2>/dev/null || true` 保留在 on_exit 成功分支。

- [ ] **Step 2: shellcheck + bash -n（本机 bash -n；shellcheck 在 CI）**

```bash
bash -n run.sh
```

- [ ] **Step 3: 删除 main_run.sh + 引用清理**

```bash
git rm main_run.sh
grep -rn "main_run" Makefile .github/ README.md tests/ docs/   # 逐处改为 run.sh
```

- [ ] **Step 4: 验证 + commit**

```bash
make check PYTHON=<venv-python>
git add -A && git commit -m "feat(launcher): rewrite run.sh as full ops CLI (profiles/pbs, resource overrides, preflight, unlock, logging); drop main_run.sh"
```

### Task 3.4: MultiQC 定制配置

**Files:**
- Create: `workflow/multiqc_config.yaml`
- Modify: `workflow/rules/upstream.smk` 的 `multiqc` 规则（shell 命令追加 `-c {params.mqc_config}`，params 加 `mqc_config=os.path.join(WORKFLOW_DIR, "multiqc_config.yaml")`）

- [ ] **Step 1: 写入 multiqc_config.yaml（对齐 rna-seq 定制模式，chip 语境）**

```yaml
## MultiQC 定制配置（全流程 QC 汇总：FastQC / Trim Galore / Bowtie2 / Picard / MACS2 / deeptools / FRiP / NSC-RSC）
title: "ChIP/CUT&Tag/ATAC/FAIRE QC Summary"
subtitle: "FastQC / Trim Galore / Bowtie2 / Picard / MACS2 / deeptools"
intro_text: "全流程 QC 汇总。FRiP 与 NSC/RSC 为流程注入的 custom content（_mqc.tsv）。"
report_header_info:
  - Workflow: "chip_cuttag_atac_faire (Snakemake)"
```

- [ ] **Step 2: 验证 + commit**

```bash
python -c "import yaml;yaml.safe_load(open('workflow/multiqc_config.yaml'))"   # venv python
make check PYTHON=<venv-python>   # multiqc 规则体渲染若被单测覆盖则同步断言
git add -A && git commit -m "feat(qc): multiqc_config.yaml customization wired into multiqc rule"
```

---

## Phase 4 — 测试与 CI 对齐（dry-run 级）

**Phase 验收：** `tests/lint.sh` 本地可跑（缺工具自动跳过）+ CI 全量；`make_testdata.py` 生成合成 chip/atac 数据集；`run_test.sh` 默认 dry-run 路径 CI 全绿；CI 两 job。

### Task 4.1: tests/lint.sh

**Files:**
- Create: `tests/lint.sh`（移植源：`D:\BioWorkflows\rna-seq\tests\lint.sh` 69 行，只读）
- Modify: `Makefile`（`make lint` 改调 `bash tests/lint.sh`；保留 `make check`）

- [ ] **Step 1: 移植适配**：结构照原件（bash -n + shellcheck 可选 + py_compile 全部 py + R parse 可选 + snakemake --lint）。适配点：脚本清单 = `run.sh tests/*.sh`；py 清单 = `tests/*.py workflow/scripts/*.py`；R 清单 = `workflow/scripts/*.R example/*.R 2>/dev/null`（可选）；snakemake --lint 用 `-s workflow/Snakefile --configfile config/config.yaml`（chip 无多 pipeline，跑一次；`--lint` 前置 `CHIP_CONFIG` 不需要）。
- [ ] **Step 2: 验证 + commit**

```bash
bash tests/lint.sh   # 本机：snakemake/shellcheck 缺失自动跳过并提示，bash -n/py_compile 必须过
git add -A && git commit -m "test: add tests/lint.sh aligned with rna-seq lint suite"
```

### Task 4.2: tests/make_testdata.py（合成数据生成器）

**Files:**
- Create: `tests/make_testdata.py`（参照源：`D:\BioWorkflows\rna-seq\tests\make_testdata.py` 219 行的确定性生成器模式，只读）

**Interfaces:**
- Produces: `--outdir DIR [--reads N]` 生成：
  - `ref/genome.fa`：两条染色体 `chr1/chr2` 各 100kb，固定随机种子（`random.Random(42)`）的 ACGS 序列；
  - `ref/genes.gtf`、`ref/genes.bed`、`ref/chrom.sizes`：每染色体 30 个"基因"（bed 12 列简化为 6 列 BED：chrom,start,end,name,score,strand；gtf 同源生成 exon 记录）；
  - `1.rawdata/{sample}_1.fq.gz` / `_2.fq.gz`：合成 PE reads（长度 50，片段 150-300 片段长度，`chr1` 上预置 3 个富集峰区各 2kb——treat 从峰区按 70% 概率取样，control 全基因组均匀），gzip 写出，样本集：`chip_treat_rep1/rep2`（narrow）、`chip_control`、`atac_treat_rep1/rep2`（atac 峰区两侧 Tn5 式分布可简化为同峰区）；
  - `samples.csv`：上述 5 样本的 6 列表（chip 组 group=`g1`，atac 组 group=`g2`，均 PE）；
  - `config.yaml`：基于 `config/config.yaml` 的测试覆盖（genome_fa/gtf/bed/chromsize 指向 ref/，genome_size: "180000"，threads: 2，qc 三开关全 true 含 nsc_rsc）。
- 全部输出确定性（固定 seed），无外部依赖（纯标准库 + gzip）。

- [ ] **Step 1: 实现**（TDD 式：先在 run_tests.py 加生成器断言组，再写实现至 check 过）。run_tests.py 断言骨架：

```python
def check_testdata(outdir):
    fq = sorted(os.listdir(os.path.join(outdir, "1.rawdata")))
    check("5 样本 × PE 共 10 个 fq.gz", len([f for f in fq if f.endswith(".fq.gz")]) == 10)
    with gzip.open(os.path.join(outdir, "1.rawdata", "chip_treat_rep1_1.fq.gz"), "rt") as fh:
        lines = [fh.readline().strip() for _ in range(4)]
    check("FASTQ 四行结构", lines[0].startswith("@") and lines[2].startswith("+") and len(lines[1]) == 50)
    rows = open(os.path.join(outdir, "samples.csv")).read().strip().splitlines()
    check("样本表 1 表头 + 5 数据行", len(rows) == 6 and rows[0].startswith("sample_id,role,group"))
    check("参考四件套齐全", all(os.path.exists(os.path.join(outdir, "ref", f))
          for f in ("genome.fa", "genes.gtf", "genes.bed", "chrom.sizes")))
    check("基因组两条染色体各 100kb", all(len(s) == 100000 for s in
          [l.strip() for l in open(os.path.join(outdir, "ref", "genome.fa")) if not l.startswith(">")]))
    # 确定性：同 seed 两次生成 genome.fa 哈希一致
```
（生成器实现要点：`random.Random(42)` 贯穿；峰区 `[(20000,22000),(50000,52000),(80000,82000)]` 固定；reads 写 gzip 用 `gzip.open(..., "wt")`。）
- [ ] **Step 2: 本地验证**

```bash
<venv-python> tests/make_testdata.py --outdir /tmp/chip-test-data --reads 2000
ls /tmp/chip-test-data/1.rawdata | head   # 10 个 fq.gz
zcat /tmp/chip-test-data/1.rawdata/chip_treat_rep1_1.fq.gz | head -4   # FASTQ 四行结构
```
- [ ] **Step 3: commit**：`git add tests/make_testdata.py tests/run_tests.py && git commit -m "test: deterministic synthetic testdata generator (chip + atac, 2x100kb genome)"`

### Task 4.3: tests/run_test.sh（dry-run 默认 + --real-run 开关）

**Files:**
- Create: `tests/run_test.sh`（移植源：`D:\BioWorkflows\rna-seq\tests\run_test.sh` 113 行，只读）

- [ ] **Step 1: 移植适配**（步骤骨架照原件：清理 → make_testdata → 组工作目录与测试 config → 跑 → 断言）：
  - `--pipeline` 参数删除（无 pipeline）；`--reads N` 默认 50000；`--keep` 保留；
  - 执行段：默认 `bash "$REPO_DIR/run.sh" -P "$WORK_DIR" -c "$WORK_DIR/config.yaml" -n -j 4`（dry-run）；`--real-run` 时不加 `-n`，跑完做输出存在性断言（`3.align/bowtie2/*_sorted.bam` 计数 = 样本数、`4.peak/g1_peaks.narrowPeak`、`4.peak/g2_peaks.narrowPeak`、`5.QC/frip/FRiP_summary.tsv`、`multiqc_report.html`、`5.QC/software_versions.yaml`）；
  - dry-run 断言：snakemake 退出码 0 且输出含 `bowtie2_mapping` 与 `callpeak_narrow` 两个规则名；
  - DAG 再生成段照原件（dot 存在才执行）。

- [ ] **Step 2: 本地验证**（无 snakemake → 只验证脚本能被 bash -n 与 usage 输出）

```bash
bash -n tests/run_test.sh && bash tests/run_test.sh --help
```

- [ ] **Step 3: commit**：`git add tests/run_test.sh && git commit -m "test: end-to-end harness (dry-run default, --real-run for server validation)"`

### Task 4.4: CI 重写（lint + dry-run 两 job）

**Files:**
- Modify: `.github/workflows/ci.yaml`（全文重写）
- Modify: `Makefile`（`make test` = check + lint.sh；删 dryrun 目标里过时引用）

- [ ] **Step 1: 新 ci.yaml**

```yaml
name: CI

on:
  push:
    branches: [master, main]
  pull_request:

jobs:
  lint:
    name: 静态检查
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: "3.11"
      - name: 安装工具
        run: |
          sudo apt-get update
          sudo apt-get install -y shellcheck
          pip install "snakemake==7.32.4" pyyaml
      - name: 单元测试
        run: make check PYTHON=python
      - name: 运行 lint
        run: bash tests/lint.sh

  dryrun:
    name: 合成数据 dry-run 回归
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: "3.11"
      - name: 安装 snakemake（dry-run 无需分析工具）
        run: pip install "snakemake==7.32.4" pyyaml
      - name: 生成合成数据并 dry-run
        run: bash tests/run_test.sh --reads 2000
```
（dry-run 不创建 conda 环境——dry-run 只构建 DAG，不执行工具；这是与 rna-seq CI 的有意差异，spec §7 已定。）

- [ ] **Step 2: Makefile 更新**：`lint: bash tests/lint.sh`；`test: check lint`；`dryrun` 目标删除（由 run_test.sh 取代）。
- [ ] **Step 3: commit**：`git add -A && git commit -m "ci: two-job pipeline (lint + synthetic-data dry-run regression)"`

---

## Phase 5 — 文档、example、清理与发布

**Phase 验收:** 文档结构与 rna-seq 对齐；死脚本归档；CHANGELOG/tag v0.4.0；全量 CI 绿。

### Task 5.1: docs/使用说明.md + CONTRIBUTING.md

**Files:**
- Create: `docs/使用说明.md`（章节对齐 rna-seq 8 章：安装与环境准备 / 快速开始 / 样本表 / 配置详解 / 运行模式与 assay 路由 / 集群提交（四 profile + auto 探测 + 资源覆盖）/ 输出解读（编号目录树 + QC 指标 NSC/RSC/FRiP 参考阈值表迁移自 README）/ FAQ）
- Create: `CONTRIBUTING.md`（移植 rna-seq 同名文件结构：开发流程 / CHANGELOG 要求 / 文档同步 checklist / 测试要求（make check + lint + run_test dry-run）/ 代码风格（snake_case、集中校验、注释双语风格说明））
- Modify: `README.md`——保留 mermaid/QC 阈值表/版本矩阵/结果路径速查，重写安装、启动（run.sh 全选项）、环境（environment.yaml + software.yaml 三条落地路径：mamba create / conda_prefix 复用 / system+lib_paths）、目录结构章节；"尚未服务器实跑"章节改为"服务器验证步骤（run_test.sh --real-run）"。

- [ ] **Step 1: 写三份文档**（内容自本计划与 spec 展开；使用说明的每条命令必须真实对应 run.sh 选项）。
- [ ] **Step 2: 交叉引用检查**

```bash
grep -rn "main_run\|workflow.smk\|envs/" README.md docs/ CONTRIBUTING.md Makefile .github/ tests/ | grep -v legacy | grep -v CHANGELOG
# 期望：无输出（历史文档 CHANGELOG/REVIEW/IMPROVEMENT_PLAN 保留原文不改写）
```
- [ ] **Step 3: commit**：`git add -A && git commit -m "docs: user guide, contributing guide, README rewrite for v0.4.0 launcher/env system"`

### Task 5.2: example/ 项目模板 + 死脚本归档

**Files:**
- Create: `example/samples.csv`（git mv 根目录 `sample_info.csv` 真实数据）、`example/config.yaml`（项目级覆盖示例：grouplist 指向 ./samples.csv + results 相关注释）、`example/README.md`（一条命令启动指引）
- Move: `scripts/`（剩余 5 个死脚本）→ `legacy/diffbind/`（git mv；`scripts/` 目录清空后删除）
- Modify: `legacy/README.md`（追加 diffbind 归档映射表：`diffpeak_DiffBind.sh/run_DiffBind.R/run_ChIPQC.R/run_chipqc_DROMPAplus.sh/annoPeak_single.R → 归档原因：未接入 DAG 的空壳，恢复路线见 docs/REVIEW.md`）

- [ ] **Step 1: 执行移动与写入**
- [ ] **Step 2: 验证 + commit**

```bash
make check PYTHON=<venv-python>
git add -A && git commit -m "chore: archive unused diffbind/chipqc scripts to legacy; add example/ project template"
```

### Task 5.3: CHANGELOG + tag

- [ ] **Step 1: CHANGELOG.md 顶部加 v0.4.0 条目**（Keep a Changelog；分 Added/Changed/Removed/Deprecated，核心：!统一环境体系、run.sh 重写、四 profile、per-rule 资源、dry-run CI、文档四件、布局迁移映射表）。
- [ ] **Step 2: 全量验证**

```bash
make check PYTHON=<venv-python> && bash tests/lint.sh && bash -n run.sh tests/*.sh
git status   # 干净
```

- [ ] **Step 3: commit + tag**

```bash
git add CHANGELOG.md && git commit -m "docs: changelog for v0.4.0"
git tag -a v0.4.0 -m "Align engineering system with rna-seq v0.8.0 (layout, unified env, run.sh, resources, tests, docs)"
```

---

## 风险提示（执行时留意）

1. **Phase 1 是全计划的地基**：任何 `REPO_DIR`/路径漏改都会在 `make check` 或 CI `--list-rules` 暴露；不要跳过 grep 核对步骤。
2. **Task 2.4 删 conda 后、Task 3.3 完成前**，`main_run.sh` 处于最小修补态——这段窗口期不要在服务器实跑（CHANGELOG 与 README 已声明 breaking）。
3. **Task 3.3 的 run.sh 移植量最大**：底稿逐节对照，先骨架后分支；`--` 透传/PIPESTATUS/trap 三处细节 rna-seq 原件已处理，勿自行简化。
4. **本机无 snakemake**：Phase 3/4 的 run_test.sh 真实执行只能在 CI 验证；本地仅 bash -n/usage 级验证，CI 首跑问题按惯例在下一 task 开头修复。
5. **提交节奏**：严格每 task 一 commit；phase 边界若 CI 红，先修红再进下一 phase。
