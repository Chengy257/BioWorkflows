# Changelog

本项目的全部显著变更记录于此。格式参考 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/)。

## [0.1.0] - 2026-09-03

### Added（新增）
- 建立 git 仓库，原样收录 `rna-seq-workflow/` 流程代码作为基线（Snakemake 四入口 + rules + scripts + 配置）。
- 标准项目脚手架：
  - `README.md`（项目说明与文档索引）
  - `docs/使用说明.md`（操作手册：数据准备/配置/运行/结果解读/FAQ/断点清单）
  - `docs/审查报告.md`（全面代码审查：P0×8 / P1×14 / P2×9，含 file:line 定位与方法学评估）
  - `docs/优化路线图.md`（三阶段改造计划、验收标准、方案对比）
  - `CHANGELOG.md`、`LICENSE`
- `.gitattributes`（脚本统一 LF 行尾，面向 Linux 集群运行）、`.gitignore`（忽略 Snakemake 运行产物）。

### Unchanged（未改动）
- 流程代码本体未做任何修改；全部修复工作按 `docs/优化路线图.md` 阶段 1 起推进。

### Known issues（已知问题）
- 详见 `docs/审查报告.md`：硬编码绝对路径、conda 环境文件缺失、count 合并外部依赖缺失、lncRNA 管线终点空壳脚本、`mapping_stat.xls` 统计错位、"control" 组名硬编码、富集环节物种硬编码。
