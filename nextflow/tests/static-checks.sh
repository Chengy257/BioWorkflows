#!/usr/bin/env bash
# Static checks for the Nextflow migration entries. This is the former
# nextflow-static CI job, kept as a local one-command script after GitHub
# Actions was retired from the repository (owner decision, 2026-10-12; design
# doc Section 15 amendment 8).
#
# Checks: Nextflow config syntax (both ext pipelines, test + test_full
# profiles), YAML validation (adopted params examples, versions.yml, local
# module meta), documentation link check, shellcheck on launchers and fixture
# scripts, py_compile on pipeline bin/ helpers, and the pytest unit suites.
#
# Tools missing from PATH are skipped with a warning (same convention as the
# per-project tests/lint.sh). Toolchain note for the HPC: nextflow lives in
# the `nf` conda env; the pytest suites run under the `ms` env python. The
# simplest full run is therefore two invocations:
#   conda run -n nf  bash nextflow/tests/static-checks.sh
#   conda run -n ms  bash nextflow/tests/static-checks.sh --unit-only
# (on a machine where one env has everything, a single plain run covers all).
set -uo pipefail

repo_root=$(cd "$(dirname "${0}")/../.." && pwd)
cd "$repo_root"

mode=${1:-all}
fail=0
note() { printf '[static-checks] %s\n' "$*"; }
have() { command -v "$1" >/dev/null 2>&1; }

# 1. Nextflow config syntax validation (chipseq-ext, rnaseq-ext, atacseq-ext)
if [ "$mode" = "all" ]; then
    if have nextflow; then
        note "nextflow config syntax validation"
        for d in nextflow/pipelines/chipseq-ext nextflow/pipelines/rnaseq-ext nextflow/pipelines/atacseq-ext; do
            (
                cd "$d" || exit 1
                nextflow config -flat -profile test > /dev/null &&
                    nextflow config -flat -profile test_full > /dev/null
            ) || fail=1
        done
    else
        note "SKIP nextflow config validation (nextflow not on PATH)"
    fi
fi

# 2. YAML validation (adopted params example, versions, module meta)
if [ "$mode" = "all" ] && have python3 && python3 -c "import yaml" 2>/dev/null; then
    note "YAML validation"
    python3 -c "import yaml,glob; files=['nextflow/versions.yml','nextflow/adopted/chipseq/params.example.yaml','nextflow/adopted/rnaseq/params.example.yaml','nextflow/adopted/differentialabundance/params.example.yaml']+glob.glob('nextflow/pipelines/chipseq-ext/modules/local/*/meta.yml')+glob.glob('nextflow/pipelines/rnaseq-ext/modules/local/*/meta.yml')+glob.glob('nextflow/pipelines/atacseq-ext/modules/local/*/meta.yml')+['nextflow/adopted/atacseq/params.example.yaml']; [yaml.safe_load(open(f)) for f in files]; print('yaml OK:', len(files), 'files')" || fail=1
fi

# 3. Documentation link check (design doc Section 11 Stage 0 mandate).
#    Cross-references in these docs are backtick-quoted paths: only
#    repo-rooted spans (nextflow/..., plans/...) and *.md spans are
#    resolvable targets; wildcards, placeholders and machine-local paths are
#    skipped, as are upstream-relative references. The two whitelisted
#    plans/ design docs are scanned as well (everything else under plans/
#    is machine-local and untracked).
if [ "$mode" = "all" ] && have python3; then
    note "documentation link check"
    python3 - <<'EOF' || fail=1
import glob, os, re, sys
SKIP_CHARS = set(' *{}<>~')
bad, scanned = [], 0
files = glob.glob('nextflow/**/*.md', recursive=True)
files += ['plans/2026-10-05-nextflow-migration-design.md',
          'plans/2026-10-05-nextflow-migration-implementation-guidance.md']
for f in files:
    scanned += 1
    text = open(f, encoding='utf-8').read()
    spans = re.findall(r'`([^`\n]+)`', text)
    spans += re.findall(r'\]\(([^)\s]+)\)', text)
    for span in spans:
        s = span.split('#', 1)[0].strip()
        if not s or s.startswith(('http://', 'https://', 'mailto:')):
            continue
        if not (s.startswith(('nextflow/', 'plans/')) or s.endswith('.md')):
            continue
        if any(ch in SKIP_CHARS for ch in s) or '...' in s or s.startswith('--'):
            continue
        if not (s.endswith(('.md', '.yml', '.yaml', '.sh', '.py', '.nf', '.csv', '.txt')) or s.endswith('/')):
            continue
        for base in (os.path.dirname(f), '.'):
            if os.path.exists(os.path.normpath(os.path.join(base, s))):
                break
        else:
            bad.append(f'{f}: unresolved path {s}')
if bad:
    print('\n'.join(sorted(set(bad))))
    sys.exit(1)
print(f'doc link check OK ({scanned} files)')
EOF
fi

# 4. Shell lint (launchers + fixture scripts)
if [ "$mode" = "all" ]; then
    if have shellcheck; then
        note "shellcheck (launchers, fixture scripts)"
        shellcheck nextflow/adopted/chipseq/run-chipseq.sh \
            nextflow/adopted/rnaseq/run-rnaseq.sh \
            nextflow/adopted/differentialabundance/run-differentialabundance.sh \
            nextflow/pipelines/chipseq-ext/tests/make_fixtures.sh \
            nextflow/pipelines/chipseq-ext/tests/make_real_bams.sh || fail=1
    else
        note "SKIP shellcheck (not on PATH)"
    fi
fi

# 5. Extension scripts compile.
#    Helper scripts live in the pipeline-level bin/ (module-level bin/ is not
#    auto-staged in this layout - see validation notes).
if [ "$mode" = "all" ] && have python3; then
    note "py_compile (pipeline bin/ helpers)"
    ok=1
    for d in nextflow/pipelines/chipseq-ext nextflow/pipelines/rnaseq-ext nextflow/pipelines/atacseq-ext; do
        for f in "$d"/bin/*.py; do
            python3 -m py_compile "$f" || ok=0
        done
    done
    [ "$ok" = 1 ] || fail=1
fi

# 6. Helper script unit tests
if python3 -c "import pytest" 2>/dev/null; then
    note "pytest unit suites"
    ok=1
    python3 -m pytest nextflow/pipelines/chipseq-ext/tests/unit -q || ok=0
    python3 -m pytest nextflow/pipelines/rnaseq-ext/tests/unit -q || ok=0
    python3 -m pytest nextflow/pipelines/atacseq-ext/tests/unit -q || ok=0
    [ "$ok" = 1 ] || fail=1
else
    note "SKIP pytest unit suites (pytest not importable by this python)"
fi

if [ "$fail" = 0 ]; then
    note "OK (all runnable checks passed; skips listed above if any)"
else
    note "FAILED"
fi
exit "$fail"
