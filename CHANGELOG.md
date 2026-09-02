# Changelog

所有对项目的显著变更将记录在本文件。格式参考 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/)。

## [Unreleased]

### 计划中（见 docs/IMPROVEMENT_PLAN.md）

- Phase 1：P0 修复（入口统一、峰调用链接入、rmDup/callpeak 规则修复、样本表 schema 统一、conda 环境拆分）
- Phase 2：可移植性（config 重构、main_run.sh 参数化、去硬编码、QC 链路打通）
- Phase 3：测试、CI、LICENSE、DiffBind 补完

## [0.2.0] - 2026-09-03

### Added

- 建立 git 仓库，规范化项目基建
- `README.md` 使用说明文档（流程总览、目录结构、快速开始、输出说明、QC 阈值）
- `docs/REVIEW.md` 全面审查报告（P0~P3 分级问题清单，含文件行号）
- `docs/IMPROVEMENT_PLAN.md` 分阶段优化计划（含验收标准）
- `.gitignore`（排除分析产出目录与大数据文件）
- `.gitattributes`（强制 LF 换行，防止 Windows CRLF 破坏 Linux 脚本）
- `.editorconfig`、`CHANGELOG.md`

## [0.1.0] - 2024-03

### Added（原始快照，commit 03f5c0c）

- ChIP-seq / CUT&Tag / ATAC-seq / FAIRE-seq 四入口 Snakemake 流程
- 上游规则：trim_galore、FastQC、MultiQC、bowtie2 比对
- 峰调用脚本：SPP 片段长估计 + MACS2（narrow/broad/ATAC 模式）、HMMRATAC、FSeq2
- MACS2 bdgcmp → bigWig 转换脚本
- ChIPseeker 批量/单样本峰注释 R 脚本
- deeptools / ChIPQC / DROMPAplus QC 脚本（部分为片段或空壳）
- 旧版整环境导出 `chip_environment.yaml` 与 PBS 启动脚本 `main_run.sh`
