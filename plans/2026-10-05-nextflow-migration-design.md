# BioWorkflows: nf-core-first Nextflow Migration Mainline

- Date: 2026-10-05
- Status: mainline design revision; implementation NOT started
- Scope: migration of the current Snakemake workflows to a simpler
  nf-core-first Nextflow workflow suite
- Primary use case: personal/research-group use on local Linux/WSL and HPC
- Authority: this document defines the migration direction and architecture.
  Detailed implementation specifications are deferred until this mainline is
  reviewed and frozen.
- Legacy baseline: existing Snakemake workflows remain available and are
  bugfix-only during migration.

## 1. Project positioning

BioWorkflows is NOT intended to become a second nf-core framework or a generic
bioinformatics workflow platform.

The target is a **personal modular workflow suite built on the nf-core
ecosystem**:

- use mature nf-core pipelines directly when they already solve the main
  analysis problem;
- extend them with small BioWorkflows-owned downstream or side workflows only
  when necessary;
- build custom DSL2 pipelines only where no suitable nf-core pipeline exists;
- reuse nf-core modules and subworkflows instead of maintaining equivalent
  local wrappers;
- keep local infrastructure deliberately small.

The value of BioWorkflows is therefore not "another collection of standard
pipelines". Its value is the small layer of decisions and capabilities that are
specific to the user's real work:

- tested pipeline revisions and preferred parameters;
- rice and other non-default reference configurations;
- local/HPC execution profiles;
- assay-specific analysis choices not covered by an upstream pipeline;
- retained project-specific analysis scripts;
- reproducible handoffs between selected workflows;
- real-data acceptance against the known legacy workflows.

## 2. Why nf-core alone is not the complete project

nf-core should provide as much of the execution stack as possible, but it does
not eliminate the need for BioWorkflows.

Current BioWorkflows capabilities fall into three categories.

### 2.1 Directly adoptable

Examples include the main execution paths for:

- BS-seq via nf-core/methylseq;
- bulk RNA-seq via nf-core/rnaseq;
- differential abundance via nf-core/differentialabundance where suitable;
- ChIP-seq via nf-core/chipseq;
- CUT&Tag via nf-core/cutandrun;
- ATAC-seq via nf-core/atacseq.

BioWorkflows should not reimplement these pipelines merely to preserve old
Snakemake structure or output paths.

### 2.2 Adoptable core plus BioWorkflows extension

Some analyses have a strong upstream pipeline but also retained local
capabilities.

Examples:

- methylseq -> local methylKit DMR analysis where retained;
- rnaseq -> project-specific lncRNA discovery;
- adopted chromatin pipelines -> selected local downstream analysis when the
  upstream pipeline does not provide the required scientific behavior.

These extensions should consume explicit upstream outputs and remain bounded.
They are not justification for forking the whole upstream pipeline.

### 2.3 Custom workflow still justified

Some current workflows have important analysis logic that is not well
represented by an existing nf-core pipeline.

Current examples include:

- plant-oriented small-RNA cascade filtering/counting;
- the current seCLIP UMI/input-control/reproducible-peak design;
- FAIRE-seq;
- selected project-specific lncRNA and DMR analyses.

These become BioWorkflows-owned DSL2 pipelines built with nf-core conventions
and reusable nf-core components wherever practical.

## 3. Frozen architecture principles

### D1. nf-core first

For every workflow entry point, evaluate the routes in this order:

1. run an existing nf-core pipeline;
2. run an existing nf-core pipeline plus a bounded local extension;
3. build a custom DSL2 pipeline from nf-core modules/subworkflows.

Do not choose route 3 merely because the legacy Snakemake workflow was custom.

### D2. One repository, flat workflow entry points

BioWorkflows remains one repository.

Each sequencing/assay type is a first-class workflow entry point. There is no
formal RNA-family, chromatin-family, or other biological-family execution
layer.

Conceptually:

```
BioWorkflows
|-- rnaseq
|-- srnaseq
|-- bsseq
|-- chipseq
|-- cuttag
|-- atacseq
|-- faireseq
|-- seclip
|-- lncrna
|-- dmr
`-- future assay/workflow entries
```

Documentation may group related workflows for navigation, but family labels do
not determine code ownership or execution.

### D3. Workflow entries compose reusable capabilities

A workflow entry owns its scientific semantics, parameters, input requirements,
and acceptance criteria.

Shared modules/subworkflows are used only when the underlying scientific
operation is actually reusable.

The project must not recreate the legacy pattern of a single workflow with a
large `seqtype` switch and many assay-dependent feature flags.

### D4. Native nf-core interfaces are preferred

Adopted pipelines keep their native:

- samplesheet schema;
- parameter schema;
- output organization;
- container/process definitions;
- release/version semantics.

BioWorkflows does not define a universal samplesheet or universal run manifest
as a prerequisite for migration.

Small one-way input helpers may be added when they save repeated manual work,
but each helper targets a specific workflow and emits that workflow's native
input format.

### D5. Define cross-workflow contracts only where a handoff exists

Do not create a global artifact model.

Define a small explicit handoff only for workflows that actually consume the
outputs of another workflow.

Initial examples:

- rnaseq -> differentialabundance;
- rnaseq -> lncrna;
- methylseq -> dmr.

Additional handoffs are introduced only when a real workflow requires them.

### D6. Container-native execution

Adopted nf-core pipelines use their native Apptainer/Singularity support.

Custom BioWorkflows DSL2 pipelines use nf-core module software definitions or
explicit per-process containers from the start.

The legacy WorkflowSpec executable resolver is migration knowledge, not the new
Nextflow runtime architecture.

A developer Conda profile may exist where useful, but production/HPC
reproducibility is container-based.

### D7. Explicit version locking

Production and acceptance runs do not use an unpinned `latest`.

A small tracked version file records the tested baseline for:

- Nextflow;
- adopted nf-core pipeline revisions;
- nf-core tooling where relevant;
- nf-test;
- important local container/module revisions.

An upstream upgrade is an explicit maintenance change and requires relevant
re-testing.

### D8. Keep the shared BioWorkflows layer small

The shared layer should initially contain only proven cross-workflow needs:

- HPC/site execution configuration;
- container/cache/mirror configuration;
- version pins;
- reference presets or reference helper configuration where shared;
- truly reusable local modules/subworkflows;
- testing conventions;
- small input/output adapters that have demonstrated repeated use.

Do not build a generic orchestration framework, plugin system, global schema,
registry, workspace model, or unified CLI during this migration.

A higher-level convenience CLI can be considered later only if repeated usage
shows a real need.

### D9. Acceptance before retirement

Each workflow entry has its own migration acceptance criteria.

Acceptance covers:

1. successful execution;
2. expected output completeness;
3. scientific consistency with the legacy workflow or an explicitly approved
   changed method;
4. reproducibility on the target execution environment.

Exact metrics and thresholds belong in the implementation specification for
that workflow.

### D10. Git-native legacy preservation

Before a legacy Snakemake implementation is removed from the active tree, its
final accepted state is preserved with an annotated Git tag/release.

Machine-local archives are optional extra backups, not the formal retirement
mechanism.

## 4. Simplified target repository organization

The exact directory names are not frozen, but the intended ownership model is:

```
BioWorkflows/
|-- <legacy Snakemake workflows>/       # retained during migration
|-- nextflow/
|   |-- conf/                           # HPC/site/container profiles
|   |-- versions.yml                    # tested tool/pipeline revisions
|   |-- adopted/                        # configs/docs/helpers only
|   |   |-- methylseq/
|   |   |-- rnaseq/
|   |   |-- differentialabundance/
|   |   |-- chipseq/
|   |   |-- cutandrun/
|   |   `-- atacseq/
|   |-- pipelines/                      # BioWorkflows-owned DSL2 pipelines
|   |   |-- srnaseq/
|   |   |-- seclip/
|   |   |-- faireseq/
|   |   |-- lncrna/
|   |   `-- dmr/
|   |-- modules/                        # local modules only when needed
|   |-- subworkflows/                   # local reusable compositions when needed
|   `-- docs/                           # migration and acceptance records
|-- plans/
`-- README.md
```

The `adopted/` directories are not forks of nf-core pipelines. They contain
only the BioWorkflows-owned material required to run a pinned upstream pipeline,
such as params, reference settings, usage notes, acceptance records, or a small
input helper.

## 5. Workflow disposition

| Workflow entry | Mainline route | BioWorkflows ownership |
|---|---|---|
| bsseq | Adopt nf-core/methylseq | params/reference/HPC configuration and acceptance |
| dmr | Custom downstream DSL2 | retained DMR logic and methylseq handoff |
| rnaseq | Adopt nf-core/rnaseq | params/reference/HPC configuration and acceptance |
| differentialabundance | Adopt nf-core/differentialabundance where suitable | explicit rnaseq handoff and preferred analysis configuration |
| lncrna | Custom side workflow | retained project-specific lncRNA discovery |
| srnaseq | Custom DSL2 | plant-oriented cascade filtering/counting and retained analysis |
| seclip | Custom DSL2 | UMI/input-control/reproducible-peak design |
| chipseq | Prefer direct nf-core/chipseq adoption | only bounded missing downstream capabilities remain local |
| cuttag | Prefer direct nf-core/cutandrun adoption | only bounded missing CUT&Tag-specific capabilities remain local |
| atacseq | Prefer direct nf-core/atacseq adoption | only bounded missing ATAC-specific capabilities remain local |
| faireseq | Custom DSL2 | FAIRE-specific workflow using reusable components where valid |

"Prefer direct adoption" is a mainline default, not a claim of exact
feature-equivalence. The implementation review must map retained legacy
capabilities to the selected upstream revision before retirement.

## 6. Decomposition of the legacy chip_cuttag_atac_faire workflow

The current combined workflow is a migration source, not a target architecture.

It is decomposed into four independent workflow entries:

```
legacy chip_cuttag_atac_faire
        |
        +--> chipseq  --> nf-core/chipseq
        +--> cuttag   --> nf-core/cutandrun
        +--> atacseq  --> nf-core/atacseq
        `--> faireseq --> BioWorkflows custom DSL2
```

The new design intentionally removes:

- one samplesheet mixing four assay types;
- one global `seqtype` router;
- global flags whose meaning changes by assay;
- assay-inapplicable stages that silently schedule nothing;
- one large mixed-assay DAG.

Potentially reusable capabilities such as annotation, motif analysis,
provenance, selected QC calculations, or output conversion may be shared only
after confirming identical semantics.

Examples of logic that remains owned by the relevant workflow include:

- ChIP control/replicate/IDR behavior;
- CUT&Tag spike-in and SEACR/MACS2 behavior;
- ATAC TSS/accessibility QC and optional footprinting;
- FAIRE-specific analysis choices.

A project may of course contain several assay types. They are run through their
respective workflow entries and compared/integrated downstream when needed.
They are not forced into one execution DAG.

The legacy combined workflow is retired only when all replacement workflow
entries required for its supported use cases have passed acceptance.

## 7. What BioWorkflows deliberately does NOT build

To keep the project appropriate for personal/research use, the migration does
not initially build:

- a universal metadata database;
- a universal run/workspace object;
- a generic artifact registry;
- a plugin/capability registry;
- a workflow-family hierarchy;
- a central orchestration service;
- a unified GUI;
- a custom replacement for nf-core schema/module tooling;
- a mandatory `bioworkflows` CLI;
- compatibility copies of all legacy output layouts.

Any of these may be reconsidered later only if a repeated concrete need emerges.

## 8. Testing and acceptance

### 8.1 Adopted pipelines

For each adopted nf-core pipeline:

- pin a tested release/revision;
- validate the local params/reference/HPC configuration;
- use upstream test profiles for basic execution;
- run BioWorkflows-owned real-data acceptance where migration equivalence
  matters;
- document deliberate differences from the legacy method.

Do not duplicate upstream module-level tests.

### 8.2 Custom pipelines

For BioWorkflows-owned DSL2 pipelines:

- follow the nf-core custom-pipeline structure where useful;
- install/reuse nf-core modules and subworkflows where possible;
- use nf-test for local modules/workflow composition;
- retain focused tests for reused Python/R scripts;
- use real-data acceptance before legacy retirement.

### 8.3 Acceptance philosophy

The purpose of migration testing is not byte-for-byte reproduction when the
upstream method has deliberately changed.

The acceptance record must distinguish:

- accidental regression;
- parameter/default differences;
- tool-version differences;
- intentional methodological changes.

Only unexplained or unacceptable differences block retirement.

## 9. Migration sequence

### Stage 0 - Minimal common foundation

Before migrating biology:

- create the minimal `nextflow/` layout;
- pin the tested Nextflow/nf-core/nf-test baseline;
- establish Apptainer and HPC profiles;
- define the lightweight adopted-pipeline directory convention;
- define the custom-pipeline convention based on nf-core tooling;
- add minimal Nextflow CI;
- define the acceptance-record template.

Do not build a universal manifest, CLI, registry, or orchestration layer.

### Stage 1 - Adoption pilot: BS-seq

Use nf-core/methylseq as the first direct-adoption case and connect the retained
DMR analysis only through the outputs it actually needs.

This validates the adopted-pipeline pattern.

### Stage 2 - Custom pilot: small-RNA

Build srnaseq as the first BioWorkflows-owned DSL2 workflow using nf-core
components where appropriate.

This validates the custom-pipeline pattern.

### Stage 3 - Remaining workflow migrations

After the two patterns above are stable, migrate the remaining workflow entries
independently.

Planned routes:

- seclip -> custom;
- chipseq -> adopted;
- cuttag -> adopted;
- atacseq -> adopted;
- faireseq -> custom;
- rnaseq -> adopted;
- differentialabundance -> adopted handoff;
- lncrna -> custom side workflow.

Their exact implementation order is operational, not architectural, and can be
chosen according to current research need and migration risk.

### Stage 4 - Legacy retirement and cleanup

Retire each legacy implementation only after its replacement has passed
acceptance.

For the combined chip_cuttag_atac_faire legacy workflow, retirement waits until
all required replacement assay entries are accepted.

Create the final legacy Git tag/release, update user documentation, and remove
obsolete Snakemake-era shared machinery only when it no longer has active
consumers.

## 10. Mainline decision rule for future workflows

When a new sequencing or assay type is needed:

```
Is there a suitable nf-core pipeline?
        |
       yes
        |
Can it be used directly with params/config?
        | yes --> adopt it
        |
        no
        |
Can the missing capability be a bounded side/downstream workflow?
        | yes --> adopt + extend
        |
        no
        v
Build a custom DSL2 workflow using nf-core components
```

This decision rule is the main architectural safeguard against BioWorkflows
growing into another monolithic or duplicative framework.

## 11. Explicitly deferred

The mainline intentionally does not freeze:

- exact Nextflow/nf-core/nf-test version numbers;
- exact per-pipeline params;
- exact reference paths;
- exact local module boundaries;
- exact acceptance thresholds;
- exact execution order after the two pilot migrations;
- whether a convenience launcher/CLI is ever needed.

These belong to later implementation specifications or to evidence from real
use.

## 12. Freeze criterion

This mainline can be frozen when an independent consistency review confirms:

1. every legacy capability has an intended migration route or an explicit
   decision to drop it;
2. no new shared abstraction duplicates an existing nf-core capability without
   a demonstrated need;
3. adopted and custom workflows remain independently runnable;
4. the legacy chip_cuttag_atac_faire workflow has been fully decomposed at the
   architecture level;
5. Stage 0 remains a minimal foundation rather than a new platform-building
   project.

Only after mainline freeze should Stage 0 implementation specifications be
written.
