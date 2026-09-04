# rna-seq 后续待做项（TODO）

> 产生于 enhancer_lncRNA_2026 项目实测准备阶段（2026-09-05），均为已识别、暂缓实施的事项。

## 1. UMI 文库支持

- **现状**：trim 规则（`workflow/rules/align.smk`）的 trim_galore 参数硬编码（`-q 30 --stringency 3 -e 0.1`），没有像 chip 工作流那样的 `trim.extra` 配置通道，无法按项目注入 cutadapt 参数；工作流无 UMI 提取/去重逻辑。
- **当前临时方案**：umi-mRNA 项目通过 `config.yaml` 的 `star_extra_args: "--clip5pNbasesRead1 9 --clip5pNbasesRead2 9"` 在比对时剥掉 5' 端 UMI（云序 oligo_dT UMI 文库结构已经数据实证：R1/R2 前 8bp 随机 UMI + 第 9bp 固定 T 锚，插入序列自第 10bp 起），不做去重。
- **待做**：
  1. 给 trim 规则增加 `trim.extra`（或等价）配置键，对齐 chip 工作流的做法；
  2. 评估集成 `umi_tools extract`（fastq 级提取 UMI 到 read name）与 `umi_tools dedup`（基于 STAR 定位去重），作为可选的 `umi:` 配置段（enable、pattern、长度、锚定碱基）。
- **数据依据**：mRNA-ZH11-0H-1 前 50 万 reads 逐位碱基组成：R1/R2 的 pos1–8 随机、pos9 = 99.5% T、pos10 起为基因组组成；caRNA 对照无此特征。

## 2. FASTQ 命名约定兼容

- **现状**：`get_fastq`（`workflow/rules/common.smk`）仅识别 `1.rawdata/{id}_1.fastq.gz + {id}_2.fastq.gz`（及单端 `{id}.fastq.gz`）；测序交付普遍使用的 `{id}_R1/_R2.fastq(.gz)`、`.fq.gz` 后缀变体不被识别。
- **当前临时方案**：在项目 `1.rawdata/` 内补一套 `{id}_1.fq.gz`/`{id}_2.fq.gz` 软链接指向同一批原始文件（纯新增，不动原始数据）。
- **待做**：在 `get_fastq` 的 pattern 列表中增加 `_R1/_R2` 与 `.fq.gz` 变体（匹配优先级需保证不歧义），或在用户文档中显著说明命名约定；相应更新 `validate_samples.py` 的报错提示。

## 3. PBS 集群支持

- **状态（2026-09-05 已解决）**：已移植 chip 的 `profile/pbs` 并更新 `run.sh` 的
  profile 白名单与 auto 检测（qsub 且无 SGE_ROOT → pbs），本服务器实测生效；
  待随仓库一并提交。
