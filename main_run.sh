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
#   -p CMD     集群提交命令（PBS 示例: "qsub -V -N chipseq -l ncpus={threads} -j oe"）
#              {threads} 由 snakemake 按每个任务实际线程数填充，与 config threads 自动对齐；
#              不传 -p 则在本机直接运行（无集群）
#   -b PATH    conda base 路径（传给 --conda-base-path，如 /opt/anaconda3）
#   -e DIR     共享 conda 环境目录（传给 --conda-prefix；强烈推荐集群使用，
#              多个工作目录复用同一套环境，避免每个项目重建约数 GB 的 11 个环境）
#   -E         仅预建 conda 环境后退出（--conda-create-envs-only；PBS 计算节点
#              无外网时，先在登录节点执行本模式再正式投递）
#   -t N       输出可见性等待秒数（--latency-wait，默认 90；PBS + 共享文件系统
#              上输出延迟是常见的假失败原因，不建议调小）
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
conda_envs_dir=""
prebuild=0
latency=90
do_rename=0
dryrun=""

usage() { grep '^#' "$0" | cut -c 3-; exit "${1:-0}"; }

while getopts "w:s:c:j:C:p:b:l:e:t:Ernh" opt; do
    case ${opt} in
        w) workdir=${OPTARG} ;;
        s) smk=${OPTARG} ;;
        c) config=${OPTARG} ;;
        j) jobs=${OPTARG} ;;
        C) cores=${OPTARG} ;;
        p) cluster_cmd=${OPTARG} ;;
        b) conda_base=${OPTARG} ;;
        e) conda_envs_dir=${OPTARG} ;;
        E) prebuild=1 ;;
        t) latency=${OPTARG} ;;
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

# snakemake 版本探测：8+ 的 conda 部署 flag 是 --software-deployment-method conda
# （--use-conda 在 8.x 为弃用别名，9.x 可能移除）；7.x 使用 --use-conda
smk_version="$(snakemake --version 2>/dev/null | head -1 || true)"
if [[ -z "${smk_version}" ]]; then
    echo "[ERROR] 未找到 snakemake 命令，请先安装（参考版本 7.32.4；8.x 亦可，脚本自动适配 flag）" >&2
    exit 1
fi
if ! [[ "${smk_version}" =~ ^[0-9]+ ]]; then
    echo "[ERROR] 无法解析 snakemake 版本号: ${smk_version}" >&2
    exit 1
fi
smk_major="${smk_version%%.*}"

cmd=(snakemake -s "${smk}" --configfile "${config}" --keep-going
     --rerun-incomplete --latency-wait "${latency}")

if (( smk_major >= 8 )); then
    cmd+=(--software-deployment-method conda)
else
    cmd+=(--use-conda)
fi

if [[ -n "${extra_config:-}" ]]; then
    cmd+=(--configfile "${extra_config}")
fi

# 注意：snakemake 中 -j 与 --cores 是同一参数（后者覆盖前者）。
# 集群模式只限并发任务数（-j）；本机模式只限总核数（--cores）。
if [[ -n "${cluster_cmd}" ]]; then
    cmd+=(-j "${jobs}" --cluster "${cluster_cmd}")
else
    cmd+=(--cores "${cores}")
fi
if [[ -n "${conda_base}" ]]; then
    cmd+=(--conda-base-path "${conda_base}")
fi
if [[ -n "${conda_envs_dir}" ]]; then
    cmd+=(--conda-prefix "${conda_envs_dir}")
fi
if [[ ${prebuild} -eq 1 ]]; then
    cmd+=(--conda-create-envs-only)
    echo "[INFO] 预建模式：仅创建 conda 环境后退出"
fi

echo "[INFO] 运行: ${cmd[*]}${dryrun:+ --dry-run}"
if [[ -n "${dryrun}" ]]; then
    cmd+=("${dryrun}")
fi
"${cmd[@]}"

# 集群模式下回收 PBS 输出日志（存在才移动）
mv ./[a-zA-Z]*.o* ./logs/ 2>/dev/null || true

# 注意：不再删除 .snakemake/ —— 它保存运行元数据与 conda 环境缓存，
# 删除会导致断点续跑与增量重跑失效。
if [[ ${prebuild} -eq 1 ]]; then
    echo "[INFO] conda 环境预建完成（目录：${conda_envs_dir:-默认位置}），可正式投递任务。"
else
    echo "[INFO] 完成。结果见 ${workdir} 下的编号目录，日志见 ${workdir}/logs/"
fi
