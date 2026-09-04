#!/usr/bin/env bash
# =====================================================================
# DROMPAplus QC（docker 版，独立工具，不在主流程 DAG 中）
# 用法: run_chipqc_DROMPAplus.sh <workdir> <qc_samples.txt> [chrom_size] [docker_image]
#   qc_samples.txt 每行: <bam相对路径> <peak相对路径>
# 默认 docker 镜像 rnakato/ssp_drompa；需本机已安装 docker 并可拉取镜像。
# =====================================================================
set -euo pipefail

DIR=${1:?用法: run_chipqc_DROMPAplus.sh <workdir> <qc_samples.txt> [chrom_size] [docker_image]}
QC_samples=${2:?缺少 qc_samples.txt}
chrom_size=${3:-chrom.sizes}
docker_image=${4:-rnakato/ssp_drompa}

cd "${DIR}"
mkdir -p tempdir_runDROMPAplus 5.QC/DROMPAplus

# chrom_size 支持绝对路径或相对本仓库 config 的文件名
if [[ ! -f "${chrom_size}" ]]; then
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
    cp "${script_dir}/config/${chrom_size}" tempdir_runDROMPAplus/ 2>/dev/null || {
        echo "[ERROR] 找不到 ${chrom_size}（既非现有路径，也不在 $(pwd) 或 config/ 下）" >&2
        exit 1
    }
else
    cp "${chrom_size}" tempdir_runDROMPAplus/
fi

while IFS=" " read -r bam peak; do
    [ -z "${bam}" ] && continue
    name=$(basename "${bam}" | cut -d_ -f1)
    docker run --rm -v "$(pwd)":/mnt "${docker_image}" parse2wig+ \
        -i "/mnt/${bam}" \
        -o "${name}" --odir /mnt/5.QC/DROMPAplus/ \
        --gt /mnt/tempdir_runDROMPAplus/"$(basename "${chrom_size}")" --pair \
        --bed "/mnt/${peak}"
done < "${QC_samples}"

echo "[INFO] DROMPAplus QC 完成: ${DIR}/5.QC/DROMPAplus/"
