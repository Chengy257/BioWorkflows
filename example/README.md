# 示例项目

本目录提供启动一个新分析项目所需的全部模板。

## 目录内容

| 文件 | 说明 |
|---|---|
| `config.yaml` | 示例项目配置（复制到项目目录后修改路径） |
| `samples_231107XTL.csv` | 真实项目样本表示例（231107XTL 项目，27 样本 / 8 处理组 + 对照），可作为样本表格式参考 |

微型测试数据集（2 组 × 2 样本降采样 reads）计划在路线图阶段 3 提供（`tests/`），当前示例项目需自备真实 fastq。

## 一条命令启动新项目

```bash
# 1. 建项目目录，复制配置与样本表
mkdir -p myproject && cd myproject
cp ../example/config.yaml .
cp ../example/samples_231107XTL.csv samples.csv   # 或换成你自己的样本表

# 2. 放置原始数据（PE 双端示例）
mkdir -p 1.rawdata
#   cp /data/raw/{sample}_1.fastq.gz 1.rawdata/
#   cp /data/raw/{sample}_2.fastq.gz 1.rawdata/

# 3. 修改 config.yaml：把 /path/to/... 替换为集群真实路径，核对样本表文件名

# 4. 校验样本表（可选但推荐）
python3 ../workflow/scripts/validate_samples.py samples.csv control

# 5. dry-run 预览
RNASEQ_PIPELINE=deg RNASEQ_CONFIG=$PWD/config.yaml \
  snakemake -s ../workflow/Snakefile --profile ../workflow/profile/default -n --quiet

# 6. 正式运行（自动选择 SGE / 本地）
bash ../run.sh deg . config.yaml 10
```

## 样本表格式

```csv
id,group,layout
WT_1,control,auto
drugA_1,drugA,auto
```

- `id`：样本名，禁止 `-`、空格、`/`；避免互为前缀；
- `group`：分组名；对照组名与 config 的 `control_group` 一致；
- `layout`：`PE`/`SE`/`auto`（当前版本以实际文件探测为准，`auto` 即可）。
