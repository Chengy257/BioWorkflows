#!/bin/bash
# Unified runtime resolver regression test; uses only temporary fake executables.
set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
T="$(mktemp -d)"
trap 'rm -rf "$T"' EXIT
mkdir -p "$T/env/bin" "$T/override/bin" "$T/rlib" "$T/project"

for tool in fastqc trim_galore multiqc STAR samtools infer_experiment.py; do
    printf '#!/bin/bash\nexit 0\n' > "$T/env/bin/$tool"
    chmod +x "$T/env/bin/$tool"
done
printf '#!/bin/bash\nexit 0\n' > "$T/override/bin/STAR-custom"
printf '#!/bin/bash\nexit 99\n' > "$T/override/bin/samtools"
chmod +x "$T/override/bin/STAR-custom" "$T/override/bin/samtools"
cat > "$T/env/bin/Rscript" <<'EOS'
#!/bin/bash
case "$*" in
  *getRversion*) printf '4.3.2' ;;
  *requireNamespace*) : ;;
  *) : ;;
esac
EOS
chmod +x "$T/env/bin/Rscript"

cat > "$T/project/config.yaml" <<'EOF2'
species: "hsa"
EOF2
cat > "$T/project/software.yaml" <<EOF2
environment:
  type: conda
  conda_prefix: "$T/env"
  strict: true
r:
  rscript: "Rscript"
  version: "4.3"
  version_check: major_minor
  lib_paths:
    - "$T/rlib"
  lib_mode: prepend
  package_sources: {}
tools:
  star: "$T/override/bin/STAR-custom"
paths: {}
databases: {}
EOF2

out=$(bash "$REPO_DIR/run.sh" -p upstream -P "$T/project" -c "$T/project/config.yaml" \
    --check-software 2>&1)
grep -F "Environment: conda ($T/env)" <<<"$out" >/dev/null
grep -F "R: 4.3.2 via $T/env/bin/Rscript" <<<"$out" >/dev/null
grep -F "R library: $T/rlib" <<<"$out" >/dev/null
grep -F "star: $T/override/bin/STAR-custom" <<<"$out" >/dev/null
grep -F "samtools: $T/env/bin/samtools" <<<"$out" >/dev/null
grep -F "Software/runtime check completed successfully." <<<"$out" >/dev/null

echo "[test] unified software/R runtime resolver + project software.yaml precedence PASS"
