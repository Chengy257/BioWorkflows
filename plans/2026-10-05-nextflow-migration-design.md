# BioWorkflows: Snakemake to Nextflow Migration Design

- Date: 2026-10-05
- Status: design agreed with the user; implementation NOT started
- Scope: all five Snakemake workflows (rna-seq, chip_cuttag_atac_faire,
  seclip-seq, srna-seq, bs-seq), the shared/ layer, tests, CI, and docs
- Machine context: written on the HPC Linux server. Repository sync checked
  2026-10-05: local main == origin/main at 8a33551 (2026-09-12), working tree
  clean, no stash.
- Document class: development-process document. plans/ is local-by-default
  (git-ignored via plans/*); this design doc is explicitly whitelisted in
  .gitignore and tracked in the repository.

## 1. Current state (verified 2026-10-05)

The monorepo is 100% Snakemake 7.32.4 based (CI pins snakemake==7.32.4 plus
pulp==2.7.0). Five workflows share one architecture and a shared/ layer.

| Workflow   | Rules | Workflow logic LOC | Scripts (py/R) | Scripts LOC | run.sh LOC |
|------------|-------|--------------------|----------------|-------------|------------|
| rna-seq    | 25    | ~2,422             | 13             | 1,134       | 535        |
| chip       | 52    | ~4,342             | 12             | 1,228       | 541        |
| seclip-seq | 24    | ~1,844             | 4              | 516         | 459        |
| srna-seq   | 13    | ~1,479             | 6              | 533         | 459        |
| bs-seq     | 15    | ~1,408             | 3              | 310         | 460        |
| Total      | 129   | ~11,495            | 38             | ~3,721      | ~2,454     |

Plus: shared/ layer ~430 LOC (launcher.sh shell helpers, WorkflowSpec runtime
framework, software version collector); per-workflow lint/dry-run test suites;
five CI jobs. chip_cuttag_atac_faire is the largest single unit (52 rules, 17
feature flags, a 1,231-line common.smk).

Architecture characteristics that matter for migration:

- Zero per-rule conda directives anywhere: one pre-activated environment per
  run; executables resolved at runtime via software.yaml plus <PREFIX>_TOOL_*
  environment variables (WorkflowSpec framework). This model is
  framework-agnostic and migrates well.
- All 38 scripts are standalone argparse CLIs with zero snakemake-object
  coupling; rules pass params and shell invocations. Scripts port to Nextflow
  unchanged.
- Deeply Snakemake-specific layers that must be rebuilt, not translated:
  parse-time Python in common.smk (config validation aggregating all errors,
  species preset backfill, conditional rule inclusion driven by config flags),
  the 4-layer configfile merge chain, expand()-based I/O declarations,
  dry-run DAG job-count baselines as tests, snakemake --lint, and the
  460-540-line run.sh launchers with SGE/PBS/SLURM auto-detection.

## 2. Decisions (agreed with the user, 2026-10-05)

- D1 Motivation: reuse the nf-core ecosystem; become cloud/container ready.
- D2 Strategy: hybrid. Adopt nf-core pipelines where coverage exists; write
  custom DSL2 pipelines (reusing nf-core modules) where it does not.
- D3 Conventions: adopted pipelines are used as shipped (their samplesheet
  format, output layout, tool versions). No compatibility shims mimicking the
  legacy layout; downstream adaptation is our own work. Legacy Snakemake
  implementations remain available during the transition.
- D4 Software environment model: custom pipelines phase in (pre-activated
  environment, then per-process conda, then containers). Adopted nf-core
  pipelines start on containers directly: a single pre-activated environment
  cannot satisfy a 40-50 process nf-core pipeline, so "use the current env
  first" only applies to custom pipelines.
- D5 HPC reality: Apptainer/Singularity usable, network access available.
  The containerization target is directly feasible.
- D6 Legacy handling: Snakemake implementations are frozen (bugfix only) and
  retired per pipeline after the Nextflow replacement passes real-data
  acceptance.

## 3. nf-core landscape facts (checked 2026-10-05)

- No released nf-core/eclip exists. The CLIP pipeline is nf-core/clipseq
  (iCount / PureCLIP / Piranha / paraclu in parallel); it does not match
  seclip-seq's core design (UMI dedup plus input-control reproducible-peak
  logic; no CLIPper).
- nf-core modules exist for the key building blocks of the custom pipelines:
  seacr_callpeaks, tobias_* (footprint / bindetect), diffbind, homer_findpeaks,
  pureclip, among 2,100+ modules installable with nf-core tooling.
- The four chromatin assays are split across nf-core/chipseq, nf-core/atacseq,
  and nf-core/cutandrun. None offers this repo's unified four-assay design
  with per-assay flags, spike-in normalization, and the QC gate table.
- nf-core/rnaseq (mature, ~50 processes), nf-core/methylseq (mature, Bismark
  based), and nf-core/smrnaseq (released; animal miRNA / miRDeep2 focus)
  cover the remaining biology to varying degrees (see disposition matrix).

## 4. Target architecture

```
BioWorkflows/
|-- rna-seq/ ... bs-seq/        # five Snakemake workflows: frozen legacy,
|                               #   retired per pipeline after NF acceptance
|-- shared/                     # kept; WorkflowSpec and launcher.sh serve
|                               #   legacy and custom NF pipelines (phase 1)
|-- plans/                      # dev-process docs: local by default,
|                               #   whitelisted design docs are tracked
|-- nextflow/                   # NEW: everything of the Nextflow era
|   |-- conf/                   # institutional config: sge.config / pbs.config
|   |   |                       #   / slurm.config / apptainer.config / mirrors
|   |   `-- adopted/            # run configs for adopted pipelines
|   |                           #   (params.yml + genome config + launchers)
|   |-- pipelines/              # custom DSL2 pipelines, nf-core-style layout
|   |   |-- srna-nf/  seclip-nf/  chip-nf/  lncrna-nf/  dmr-nf/
|   |-- modules/                # local modules shared across custom pipelines
|   `-- docs/                   # migration matrix, acceptance records
`-- AGENTS.md / README.md       # governance updated; English / LF /
                                # conventional-commit rules all carry over
```

Adopted nf-core pipelines are run (not forked) with pinned versions, custom
params YAML, and the institutional config. Custom pipelines follow a
lightweight nf-core template: main.nf, workflows/, subworkflows/, modules/,
conf/, assets/ (samplesheet schema), tests/.

## 5. Per-workflow disposition

| Workflow | Route | Rationale |
|----------|-------|-----------|
| bs-seq   | Adopt nf-core/methylseq + custom dmr-nf add-on | Bismark engine matches exactly; the 13-round Bismark flag contract knowledge carries over. methylKit DMR has no nf-core home; existing run_dmr.R is a standalone CLI and is reused unchanged. |
| rna-seq  | Adopt nf-core/rnaseq (DE via differentialabundance) + custom lncrna-nf side pipeline | Upstream and DE are covered; optional StringTie assembly exists. CNCI / Pfam / NR lncRNA identification must be custom; scripts reused. |
| chip     | Custom chip-nf DSL2 + installed nf-core modules (seacr, tobias, diffbind, homer) | The unified four-assay design, spike-in normalization, and QC gate table do not exist in any nf-core pipeline; splitting into three pipelines would fragment them. Largest unit (52 rules), scheduled last. |
| seclip   | Custom seclip-nf (pureclip and related modules) | UMI + input-control reproducible-peak logic is in-house; clipseq does not match. Scripts reused. |
| srna     | Custom srna-nf | Plant small-RNA biology (rice osa presets; cascade filter with per-class counting) does not match smrnaseq's animal miRNA / miRDeep2 focus; only 13 rules, cheap to port. |

## 6. Software environment model (phased)

| Phase | Custom pipelines | Adopted pipelines |
|-------|------------------|-------------------|
| 1 (now) | Pre-activated environment + software.yaml; WorkflowSpec runtime resolution reused as-is | apptainer/singularity profile from day one |
| 2 | Per-process conda directives (auto-created envs) | already containerized |
| 3 | Apptainer containers | - |

Note on the HPC: conda lives at ~/soft/miniconda3 (base python 3.10); the WSL
project analysis envs are not provisioned on this server, and snakemake is not
in the base env. Phase 1 on this machine therefore requires either
provisioning envs here or running phase 1 validation on the WSL side; decide
at implementation start of each pipeline.

## 7. Implementation phases

### Phase 0 - Foundations (first PR-sized deliverable)

1. Create the nextflow/ directory skeleton and the custom-pipeline template
   (lightweight nf-core structure).
2. Institutional configs: SGE / PBS / SLURM executor configs aligned with the
   current run.sh policy (queue, h_vmem-style memory resource, runtime), plus
   apptainer.config and mirror settings; validate on this HPC.
3. Shared component: samplesheet converter (single-column sample_id CSV ->
   nf-core samplesheet format), preserving the aggregate-validation style of
   common.smk.
4. Testing basis: -stub-run DAG assertions plus nf-test for modules; CI gains
   a Nextflow job (Java + Nextflow install); the five legacy CI jobs stay as
   they are.
5. AGENTS.md governance additions for the Nextflow era; this design doc is
   the reference for the migration matrix.

### Phase 1 - Pilot A: the adoption path (bs-seq -> nf-core/methylseq)

- Samplesheet conversion, rice/osa and custom-reference genome config,
  -profile apptainer/singularity runs on WSL (or HPC) and on the HPC.
- dmr-nf small pipeline (phase 1 environment model): wires the existing
  run_dmr.R unchanged via WorkflowSpec.
- Acceptance: cross-compare against legacy bs-seq outputs on real data
  (Bismark coverage / methylKit DMR consistency).

### Phase 2 - Pilot B: the custom path (srna-nf)

- Port all 13 rules; establish the standard custom-pipeline idioms: channel
  model, includeConfig layered config chain, Groovy aggregate validation,
  species preset backfill, feature flags driving conditional processes.
- Reuse all 6 scripts unchanged (phase 1 environment model); pytest keeps
  covering the script layer.
- Deliverable: the finalized custom-pipeline template that the later
  pipelines copy.

### Phase 3 - Main custom builds (seclip-nf, then chip-nf)

- seclip-nf: UMI + STAR + pureclip module + in-house reproducible-peak logic
  (input control, filter by input).
- chip-nf: unified four-assay DSL2; all 17 feature flags parameterized;
  seacr_callpeaks / tobias_* / diffbind / homer modules installed and
  adapted; QC gate table via gates_summary.py reuse. Deliver in batches:
  ChIP main chain first, then spike-in, IDR/consensus, QC gates, footprint.
- Each pipeline: WSL real-run, then HPC run, then cross-acceptance against
  the legacy Snakemake outputs.

### Phase 4 - Adopted heavy pipeline (rna-seq)

- nf-core/rnaseq adoption (differentialabundance for DE/enrichment);
  lncrna-nf side pipeline custom (StringTie assembly -> CNCI -> Pfam/NR).
- Map batch_correction and other in-house capabilities to nf-core
  equivalents; record deliberate differences in nextflow/docs/.

### Phase 5 - Containerization sweep + legacy retirement

- Switch custom pipelines from the phase-1 environment model to apptainer
  profiles, one pipeline at a time.
- Per pipeline, after NF acceptance: freeze the Snakemake version, archive it
  outside the repo per AGENTS.md, update README / user-guide / CHANGELOG.
- Replace the dry-run job-count baselines in AGENTS.md with the Nextflow
  equivalents as each pipeline retires.

## 8. Testing and acceptance strategy

- Custom pipelines: -stub-run assertions that key processes appear / do not
  appear per feature-flag scenario (mirrors today's DAG assertions), nf-test
  for modules, pytest for scripts. One scenario per feature flag.
- Adopted pipelines: nf-core test profile on small data to prove plumbing,
  then one real-data full acceptance run.
- Cross-acceptance: every NF pipeline is compared against the legacy
  Snakemake run on the same samples (count tables / peak counts / DMR loci).
  Differences must be explainable before legacy retirement.
- CI evolves incrementally; legacy jobs stay green but frozen.

## 9. Risks and mitigations

- nf-core version drift: pin adopted pipeline versions; upgrades go through
  dedicated PRs with re-acceptance.
- Non-model reference genomes (rice/osa and custom): igenomes has no
  coverage; maintain genome configs centrally under nextflow/conf/adopted/.
- chip-nf size: 52 rules and 17 flags - deliver in batches (see Phase 3).
- R package environment differences (Bioconductor versions inside nf-core
  containers): check DESeq2 / DiffBind / methylKit numeric differences during
  cross-acceptance.

## 10. First milestone (upon implementation start)

Phase 0 in full, plus the start of Phase 1: directory skeleton, the three
institutional executor configs validated on this HPC, samplesheet converter,
AGENTS.md governance update, and the dmr-nf skeleton wiring run_dmr.R. The
five legacy workflows remain untouched in this phase.

## 11. References

- nf-core/clipseq: https://nf-co.re/clipseq/1.0.0
- nf-core modules browser: https://nf-co.re/modules
- homer_findpeaks module: https://nf-co.re/modules/homer_findpeaks
- pureclip module:
  https://github.com/nf-core/modules/blob/master/modules/nf-core/pureclip/meta.yml
