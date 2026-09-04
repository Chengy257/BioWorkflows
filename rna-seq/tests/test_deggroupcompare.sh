#!/bin/bash
#########################################################################
# DEGgroupCompare.sh 参数贯通回归测试（无需 R/snakemake）
#   用记录参数的 Rscript shim 验证 species 贯通到每个 run_deg_compare.R 任务；
#   同时验证临时文件不残留在运行目录。
# 用法: bash tests/test_deggroupcompare.sh
# 依赖: bash、python3（或 python）、mktemp；结果写入临时目录并自动清理。
#########################################################################
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TARGET="$REPO_DIR/workflow/scripts/DEGgroupCompare.sh"

T="$(mktemp -d)"
trap 'rm -rf "$T"' EXIT
mkdir -p "$T/cwd" "$T/bin" "$T/results/5.DEG/DEGs" "$T/results/6.DEGcompare"
: > "$T/args.log"

## ---- 测试夹具：样本表（2 个处理组 + 对照）、注释表、两个交集文件 ----
## 交集 1 的基因存在于注释表（fgrep 成功路径）；交集 2 的基因不存在（fgrep 失败路径，
## 旧实现会残留 tmp.gene 临时文件）
cat > "$T/samples.csv" <<'EOF'
id,group
s1,treatA
s2,treatB
s3,control
EOF
printf 'geneX\tannotation of geneX\ngeneY\tannotation of geneY\n' > "$T/anno.tsv"
printf 'treatA\ttreatB\t2\ttreatA\342\210\251treatB\tgeneX geneY\n' \
    > "$T/results/6.DEGcompare/treatA_treatB_intersetion.xls"
printf 'treatA\ttreatC\t2\ttreatA\342\210\251treatC\tgeneZ\n' \
    > "$T/results/6.DEGcompare/treatA_treatC_intersetion.xls"

## ---- shims：python 规避本机 site 问题；Rscript 逐行记录参数 ----
PY=python3
command -v python3 >/dev/null 2>&1 || PY=python
printf '#!/bin/bash\nexec %s -S "$@"\n' "$PY" > "$T/bin/python"
printf '#!/bin/bash\nprintf "RLIB:%%s\\n" "${R_LIBS_USER:-}" >> %q\nprintf "ARG:%%s\\n" "$@" >> %q\n' "$T/args.log" "$T/args.log" > "$T/bin/Rscript"
chmod +x "$T/bin/python" "$T/bin/Rscript"

fails=0
check() {  ## check <描述> <命令...>
    local desc=$1; shift
    if "$@"; then echo "  PASS  $desc"; else echo "  FAIL  $desc"; fails=$((fails + 1)); fi
}

run_case() {  ## run_case <species>
    local spe=$1
    : > "$T/args.log"
    rm -f "$T/cwd/tmp.gene"
    (cd "$T/cwd" && PATH="$T/bin:$PATH" RNASEQ_RSCRIPT="$T/bin/Rscript" RNASEQ_PYTHON="$T/bin/python" R_LIBS_USER="$T/rlibA:$T/rlibB" "$TARGET" \
        "$T/results" control "$T/anno.tsv" "$T/samples.csv" 2 "$spe" \
        >/dev/null 2>&1)
}

echo "[test] 用例 1：species=hsa 贯通到全部任务"
run_case hsa
## 2 处理组 → 1 对组合 × UP/DOWN/ALL 3×3 = 9 个任务
check "9 个任务全部收到 species=hsa" \
    test "$(grep -c '^ARG:hsa$' "$T/args.log")" = 9
check "9 个 R 子任务继承 R_LIBS_USER" \
    test "$(grep -c "^RLIB:$T/rlibA:$T/rlibB$" "$T/args.log")" = 9
check "flag.log 产出" test -f "$T/results/6.DEGcompare/flag.log"
check "交集基因注释产出（fgrep 成功路径）" \
    test -f "$T/results/6.DEGcompare/intersetionGene_anno/treatA_treatB_intersetionGene_anno.xls"
check "注释内容含交集基因" \
    grep -q geneX "$T/results/6.DEGcompare/intersetionGene_anno/treatA_treatB_intersetionGene_anno.xls"
check "fgrep 失败路径无 tmp.gene 残留" test ! -e "$T/cwd/tmp.gene"

echo "[test] 用例 2：species=osa 参数贯通"
run_case osa
check "9 个任务全部收到 species=osa" \
    test "$(grep -c '^ARG:osa$' "$T/args.log")" = 9
check "flag.log 产出" test -f "$T/results/6.DEGcompare/flag.log"

if [[ "$fails" -gt 0 ]]; then
    echo "[test] $fails 项断言未通过 ✘"
    exit 1
fi
echo "[test] 全部通过 ✔ （DEGgroupCompare 参数贯通 + 临时文件治理）"
