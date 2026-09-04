import sys

# rna-seq-workflow 旧入口已于 v0.3.0 迁移至 workflow/Snakefile（pipeline 参数选择）
sys.stderr.write(
    "[DEPRECATED] 该入口文件已废弃（v0.3.0 迁移至标准布局）。\n"
    "请改用唯一入口 workflow/Snakefile：bash run.sh as <project_dir> [config.yaml] [jobs]\n"
    "详见 README.md 与 docs/使用说明.md。\n"
)
sys.exit(1)
