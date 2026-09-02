#!/usr/bin/env bash
# =====================================================================
# chip_cuttag_atac_faire 启动脚本
#
# 用法:
#   bash main_run.sh -w /path/to/workdir [选项]
#
# 选项:
#   -w DIR     数据工作目录（必选，含 1.rawdata/ 与 sample_info.csv）
#   -s FILE    入口 Snakefile（默认仓库内 workflow.smk）
#   -c FILE    配置文件（默认仓库内 config/config.yaml）
#   -j N       集群并发任务数（默认 3）
#   -C N       本机总核数（默认 18）
#   -p CMD     集群提交命令（PBS 示例: "qsub -V -N chipseq -l ncpus=6 -j oe"）
#              不传 -p 则在本机直接运行（无集群）
#   -b PATH    conda base 路径（传给 --conda-base-path，如 /opt/anaconda3）
#   -l FILE    额外配置文件（如 config.local.yaml，后加载者覆盖前者的键；
#              不传时自动检测工作目录下的 config.local.yaml，存在即叠加）
#   -r         启动前把 1.rawdata 中的常见 R1/R2 后缀统一重命名为 _1/_2.fq.gz
#   -n         dry-run，仅打印 DAG 不执行
# =====================================================================
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
workdir=""
smk="${REPO_DIR}/workflow.smk"
config="${REPO_DIR}/config/config.yaml"
jobs=3
cores=18
cluster_cmd=""
conda_base=""
do_rename=0
dryrun=""

usage() { grep '^#' "$0" | cut -c 3-; exit "${1:-0}"; }

while getopts "w:s:c:j:C:p:b:l:rnh" opt; do
    case ${opt} in
        w) workdir=${OPTARG} ;;
        s) smk=${OPTARG} ;;
        c) config=${OPTARG} ;;
        j) jobs=${OPTARG} ;;
        C) cores=${OPTARG} ;;
        p) cluster_cmd=${OPTARG} ;;
        b) conda_base=${OPTARG} ;;
        l) extra_config=${OPTARG} ;;
        r) do_rename=1 ;;
        n) dryrun="--dry-run" ;;
        h) usage 0 ;;
        *) usage 1 ;;
    esac
done

if [[ -z "${workdir}" ]]; then
    echo "[ERROR] 必须通过 -w 指定数据工作目录" >&2
    usage 1
fi
if [[ ! -d "${workdir}" ]]; then
    echo "[ERROR] 工作目录不存在: ${workdir}" >&2
    exit 1
fi

cd "${workdir}"
mkdir -p logs

# 可选：统一原始数据命名（perl rename 语法；系统无 rename 则跳过并提示）
if [[ ${do_rename} -eq 1 && -d "1.rawdata" ]]; then
    if command -v rename >/dev/null 2>&1; then
        (
            cd 1.rawdata
            rename _R1.fastq.gz _1.fq.gz ./*gz 2>/dev/null || true
            rename _R2.fastq.gz _2.fq.gz ./*gz 2>/dev/null || true
            rename _1.fastq.gz _1.fq.gz ./*gz 2>/dev/null || true
            rename _2.fastq.gz _2.fq.gz ./*gz 2>/dev/null || true
        )
    else
        echo "[WARN] 未找到 rename 命令，跳过重命名；请确保 fastq 命名为 {sample}_1.fq.gz / {sample}_2.fq.gz"
    fi
fi

# 额外配置：显式 -l 优先；否则自动检测工作目录下的 config.local.yaml
if [[ -z "${extra_config:-}" && -f "config.local.yaml" ]]; then
    extra_config="config.local.yaml"
    echo "[INFO] 检测到 config.local.yaml，将叠加覆盖默认配置"
fi

cmd=(snakemake -s "${smk}" --configfile "${config}"
     --use-conda --keep-going
     -j "${jobs}" --cores "${cores}")

if [[ -n "${extra_config:-}" ]]; then
    cmd+=(--configfile "${extra_config}")
fi

if [[ -n "${cluster_cmd}" ]]; then
    cmd+=(--cluster "${cluster_cmd}")
fi
if [[ -n "${conda_base}" ]]; then
    cmd+=(--conda-base-path "${conda_base}")
fi

echo "[INFO] 运行: ${cmd[*]}"
"${cmd[@]}" ${dryrun}

# 集群模式下回收 PBS 输出日志（存在才移动）
mv ./[a-zA-Z]*.o* ./logs/ 2>/dev/null || true

# 注意：不再删除 .snakemake/ —— 它保存运行元数据与 conda 环境缓存，
# 删除会导致断点续跑与增量重跑失效。
echo "[INFO] 完成。结果见 ${workdir} 下的编号目录，日志见 ${workdir}/logs/"
