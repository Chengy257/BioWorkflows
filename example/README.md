# 示例项目

本目录提供启动一个新分析项目所需的全部模板。

## 目录内容

| 文件 | 说明 |
|---|---|
| `samples.csv` | 真实项目样本表示例（myc / IgG 两样本，chip narrow，组 `myc_vs_IgG`），可作为样本表格式参考 |
| `config.yaml` | 项目级覆盖示例，可直接 `-c` 传入；未列出的键沿用仓库 `config/config.yaml` 默认值 |

微型测试数据集（合成 fastq）见 `tests/make_testdata.py`（`bash tests/run_test.sh` 自动组装工作目录并 dry-run），正式分析项目需自备真实 fastq。

## 三条命令启动新项目

```bash
# 1) 建工作目录，复制样本表（config 可直接用仓库内 example/config.yaml，也可复制后修改）
mkdir -p ~/work/demo
cp example/samples.csv ~/work/demo/samples.csv

# 2) 放置原始数据（PE 双端：<sample_id>_1.fastq.gz / <sample_id>_2.fastq.gz；
#    常见 R1/R2 后缀可 bash run.sh -r 批量重命名）
mkdir -p ~/work/demo/1.rawdata
#   cp /data/raw/myc_R1.fastq.gz ~/work/demo/1.rawdata/myc_1.fastq.gz
#   cp /data/raw/myc_R2.fastq.gz ~/work/demo/1.rawdata/myc_2.fastq.gz

# 3) dry-run 预览 → 正式运行（在仓库根目录执行；先在 config 中替换参考基因组路径）
bash run.sh -P ~/work/demo -n
bash run.sh -P ~/work/demo -c example/config.yaml --profile auto -j 10
```

说明：

- `-c example/config.yaml` 显式指定项目配置；若把 `example/config.yaml` 复制为工作目录内的 `config.yaml`，可省略 `-c`（自动探测）。`config.grouplist` 相对**工作目录**解析，故样本表须放在工作目录内（或用绝对路径）。
- 启动前建议先 `bash run.sh -P ~/work/demo --check-software` 预检工具与 R 环境。
- 集群提交示例：`--profile pbs --queue workq --memory 16G --runtime 600`，完整选项见 `run.sh --help`。

## 样本表格式

6 列 schema 模板与逐行校验规则见 [config/samples.csv](../config/samples.csv)：

```csv
sample_id,role,group,seqtype,layout,peak_type
myc,treat,myc_vs_IgG,chip,PE,narrow
IgG,control,myc_vs_IgG,chip,PE,narrow
```

- `sample_id`：样本名，禁止 `-`、空格、`/`、双下划线，避免互为前缀；
- `role`：`treat` / `control`（无对照的组可只写 treat 行）；
- `group`：分组名，每个组内 `seqtype`、`peak_type` 必须一致；
- `seqtype`：`chip` / `cuttag` / `atac` / `faire`，混型项目按行自动路由；
- `layout`：当前版本仅支持 `PE`；
- `peak_type`：仅 chip 必填（`narrow` / `broad`），atac/faire 留 `none`。
