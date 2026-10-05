# BioWorkflows: Snakemake to Nextflow Migration Mainline Design

- Date: 2026-10-05
- Status: mainline design revision; implementation NOT started
- Scope: migration of the five legacy Snakemake workflow families, the shared
  execution layer, testing, CI, documentation, and legacy retirement
- Authority: this document defines migration direction and architectural
  boundaries. Detailed implementation specifications are intentionally deferred
  until the mainline is reviewed and frozen.
- Legacy baseline: current Snakemake implementations remain available and are
  bugfix-only during migration.

## 1. Purpose

BioWorkflows is moving from a Snakemake-centered monorepo to a Nextflow DSL2
workflow collection with three goals:

1. reuse mature nf-core pipelines and modules instead of rebuilding solved
   infrastructure;
2. make HPC, Apptainer/container, and future cloud execution first-class; and
3. use migration to correct legacy architectural coupling rather than translate
   the old DAGs mechanically.

The migration is therefore a scientific and architectural rebaseline, not a
line-by-line Snakemake-to-Nextflow translation.

## 2. Current-state observations that constrain the migration

The repository currently contains five first-class Snakemake workflows:
rna-seq, chip_cuttag_atac_faire, seclip-seq, srna-seq, and bs-seq. They share
runtime resolution, launcher conventions, scheduler handling, configuration
patterns, and standalone Python/R scripts.

Important migration properties:

- The existing Python/R analysis scripts are mostly standalone CLI programs and
  can usually be reused without Snakemake coupling.
- Parse-time Python, config merging, conditional rule inclusion, target
  expansion, scheduler launchers, and dry-run DAG baselines are Snakemake
  architecture and must be redesigned rather than translated.
- The current chip_cuttag_atac_faire workflow is a legacy aggregation of four
  biologically distinct assay types. Its common sample table routes
  ChIP-seq, CUT&Tag, ATAC-seq, and FAIRE-seq through one DAG, while many
  downstream stages are already assay-specific. This coupling must not be
  reproduced in the Nextflow architecture.
- AGENTS.md is machine-local and git-ignored. Repository-level migration
  governance must therefore live in tracked project documentation, not in
  AGENTS.md.

## 3. Frozen migration principles

### D1. Hybrid adoption strategy

Prefer an existing nf-core pipeline when it covers the biological execution
path sufficiently well. Use custom DSL2 only for capabilities that are genuinely
project-specific or not covered by a suitable upstream pipeline.

### D2. Adopt upstream pipelines; do not fork them by default

Adopted nf-core pipelines are run at pinned revisions with BioWorkflows-owned
input adapters, parameter sets, reference configuration, executor profiles, and
acceptance records.

Do not fork an nf-core pipeline merely to preserve the legacy BioWorkflows
output layout or CLI. Deliberate downstream adaptation is preferred over
compatibility shims.

### D3. Assay-first decomposition

The legacy chip_cuttag_atac_faire workflow will NOT migrate to one monolithic
`chip-nf` pipeline.

ChIP-seq, CUT&Tag, ATAC-seq, and FAIRE-seq become independent execution units
with explicit assay contracts. They may reuse shared modules and subworkflows,
but assay-specific decisions remain in the assay that owns them.

A new Nextflow run has one primary assay contract. Mixed-assay projects may
contain several assay runs, but cross-assay coordination occurs above the
individual execution DAGs rather than by mixing all assays in one sample table.

### D4. Share capabilities, not giant workflows

Reusable capabilities should be factored into shared DSL2 modules/subworkflows
or reusable post-processing components where scientifically appropriate.

Examples include:

- reference and input validation;
- alignment/QC primitives where parameters are genuinely compatible;
- peak-format normalization and annotation;
- MultiQC integration;
- motif analysis;
- differential-binding primitives;
- provenance/version capture;
- selected QC summaries and reporting.

Assay-specific logic must not be generalized merely to remove duplication.

Examples that stay assay-specific include:

- ATAC-specific TSS and accessibility QC;
- ATAC/FAIRE footprinting;
- CUT&Tag/CUT&RUN-style SEACR and spike-in semantics;
- ChIP control and replicate/IDR policy;
- assay-specific deduplication and peak-calling defaults.

### D5. Container-native Nextflow target

New Nextflow execution should be container-native from the start wherever
practical.

- Adopted nf-core pipelines use their native Apptainer/Singularity support.
- Custom DSL2 processes use nf-core module software definitions or explicit
  per-process containers.
- A local Conda/developer profile may exist as a secondary execution option.
- The legacy WorkflowSpec executable resolver is not a long-term Nextflow
  runtime contract. It may be consulted for migration knowledge and temporarily
  used for isolated in-house scripts only when a container is not yet available.

This removes the previous plan to rebuild custom pipelines first around one
pre-activated environment and containerize them later.

### D6. Explicit version locking

No production or acceptance run uses an unpinned `latest`.

Phase 0 establishes a tracked version manifest covering at minimum:

- tested Nextflow version;
- nf-core pipeline revisions;
- nf-core tooling revision where relevant;
- nf-test revision;
- locally maintained module/container revisions when applicable.

Upgrades are explicit changes followed by re-validation.

### D7. Normalized input contract plus pipeline-specific adapters

BioWorkflows will not rely on a universal
`sample_id -> nf-core samplesheet` converter.

Instead, a normalized BioWorkflows run manifest captures the information
actually needed to describe a run, including as applicable:

- sample identity;
- assay;
- FASTQ paths and layout;
- biological group/condition/batch;
- treatment/control or control linkage;
- strandedness or assay-specific metadata;
- reference selection and optional run metadata.

Pipeline-specific adapters validate this manifest and emit the exact
samplesheet/parameter contract required by the selected adopted or custom
pipeline.

The normalized manifest is a BioWorkflows boundary. The upstream nf-core input
formats remain unchanged.

### D8. Acceptance precedes retirement

Every migrated execution unit must have a pipeline-specific acceptance contract
before real-data cross-comparison begins.

"Results are explainable" is not sufficient as a retirement gate.

Acceptance must cover four dimensions:

1. execution completeness;
2. expected output/artifact completeness;
3. scientific consistency with the legacy baseline or an explicitly approved
   changed method;
4. reproducibility and target-HPC execution.

Exact metrics and tolerances are defined in the later implementation
specification for each execution unit.

### D9. Legacy preservation is Git-native

The final accepted Snakemake baseline is preserved by Git tag/release before its
removal from the active tree.

Machine-local archives may exist as extra backups, but they are not the formal
project retirement mechanism.

A combined legacy workflow such as chip_cuttag_atac_faire is not retired until
all execution units required to replace its supported assays have independently
passed acceptance.

## 4. Target repository architecture

The exact file layout may be refined during Phase 0, but the ownership model is
frozen as follows:

```
BioWorkflows/
|-- rna-seq/ ... bs-seq/              # frozen legacy Snakemake during transition
|-- shared/                            # legacy shared code; reusable generic code
|-- plans/
|   `-- 2026-10-05-nextflow-migration-design.md
|-- nextflow/
|   |-- conf/                          # executor/container/site configuration
|   |-- contracts/                     # normalized run manifest + artifact contracts
|   |-- versions/                      # tested/pinned tool and pipeline revisions
|   |-- adopted/                       # BioWorkflows configs/adapters for nf-core runs
|   |   |-- methylseq/
|   |   |-- rnaseq/
|   |   |-- chipseq/
|   |   |-- atacseq/
|   |   `-- cutandrun/
|   |-- pipelines/                     # BioWorkflows-owned DSL2 pipelines
|   |   |-- srna-nf/
|   |   |-- seclip-nf/
|   |   |-- faire-nf/
|   |   |-- lncrna-nf/
|   |   `-- dmr-nf/
|   |-- modules/                       # local reusable processes only when needed
|   |-- subworkflows/                  # shared compositions / post-processing blocks
|   `-- docs/                          # migration matrix + acceptance records
`-- README.md
```

The directory names above express component ownership, not a requirement that
every listed directory be created in Phase 0.

## 5. Migration disposition by biological execution unit

| Legacy capability | Mainline Nextflow route | Notes |
|---|---|---|
| BS-seq core | Adopt nf-core/methylseq | High functional overlap, but legacy Bismark parameters and output semantics require explicit mapping and real-data acceptance. |
| BS-seq DMR | Custom dmr-nf | Reuse the existing methylKit analysis logic where scientifically retained; consume a defined methylation artifact contract rather than legacy paths. |
| Bulk RNA-seq core | Adopt nf-core/rnaseq | Use pinned upstream behavior and BioWorkflows input/reference adapters. |
| RNA differential analysis | Adopt nf-core/differentialabundance where suitable | Treat rnaseq -> differentialabundance as an explicit artifact handoff, not an internal rnaseq stage. Preserve custom analysis only where upstream capability is insufficient. |
| lncRNA discovery | Custom lncrna-nf | CNCI/Pfam/NR and other retained project-specific logic remain a side pipeline consuming explicit RNA-seq artifacts. |
| small-RNA | Custom srna-nf | Plant-oriented cascade filtering and counting remain project-specific; this is the first custom DSL2 pilot. |
| seCLIP | Custom seclip-nf | Preserve UMI/input-control/reproducible-peak logic while reusing upstream modules where appropriate. |
| ChIP-seq | Prefer adopted nf-core/chipseq + BioWorkflows add-ons | Do not inherit unrelated ATAC/FAIRE/CUT&Tag branches. Missing retained capabilities are added as bounded downstream components rather than by forking upstream by default. |
| CUT&Tag | Prefer adopted nf-core/cutandrun + BioWorkflows add-ons | CUT&Tag support, controls, spike-in and peak calling make this a strong adoption route; retained BioWorkflows behavior still requires explicit mapping/acceptance. |
| ATAC-seq | Prefer adopted nf-core/atacseq + BioWorkflows add-ons | ATAC-specific QC and optional footprinting remain explicitly owned by the ATAC route. |
| FAIRE-seq | Custom faire-nf | No mature upstream replacement is assumed. Reuse appropriate accessibility/chromatin modules without pretending FAIRE is ATAC. |

"Prefer adopted" means adoption is the mainline design. A later implementation
review may fall back to a thin custom DSL2 pipeline only if real contract gaps
make upstream adoption scientifically or operationally unsuitable.

## 6. Chromatin-family decomposition

### 6.1 What is being removed from the legacy design

The following legacy behavior is intentionally not preserved as a core
execution contract:

- one sample sheet containing ChIP, CUT&Tag, ATAC, and FAIRE rows;
- one `ASSAYS=(chip, cuttag, atac, faire)` branch table driving a single DAG;
- global feature flags whose meaning changes by assay;
- optional stages silently producing jobs for only a subset of assay groups;
- a single 50+ rule workflow whose complexity is dominated by combinations of
  assay and feature flags.

This is an architectural breaking change and is intentional.

### 6.2 New ownership model

Each assay route owns:

- its required input metadata;
- assay-specific validation;
- upstream execution;
- default QC;
- peak/signal semantics;
- its own acceptance baseline.

Shared components may then be composed where appropriate.

A conceptual model is:

```
normalized project/run manifest
        |
        +-- ChIP run ------> chipseq route -----+
        +-- CUT&Tag run ---> cutandrun route ---+--> shared compatible post-processing
        +-- ATAC run ------> atacseq route -----+    and normalized artifacts
        +-- FAIRE run -----> faire-nf route ----+
                                                   |
                                                   v
                                      optional cross-assay integration
```

Cross-assay integration is downstream of assay execution. It is not a reason to
recombine assay-specific processing into one pipeline.

### 6.3 Shared chromatin/accessibility capabilities

Candidates for reusable local modules/subworkflows include:

- common artifact metadata and provenance;
- peak BED/narrowPeak/broadPeak normalization;
- peak annotation;
- generic FRiP-like calculations where definitions are harmonized;
- motif analysis;
- selected differential-binding preparation;
- MultiQC/custom summary integration;
- common reference utilities.

Reuse is conditional on identical scientific semantics. Similar-looking steps
with different assay assumptions remain separate.

### 6.4 Assay-specific retained capabilities

ChIP-seq:
- treatment/control semantics;
- narrow/broad peak behavior;
- replicate-aware/IDR behavior where retained;
- differential binding and motif routes as applicable.

CUT&Tag:
- CUT&Tag-appropriate peak calling;
- SEACR/MACS2 mapping where retained;
- spike-in/control behavior;
- assay-appropriate duplicate handling.

ATAC-seq:
- accessibility-specific peak behavior;
- TSS enrichment;
- organelle fraction reporting where relevant;
- optional TOBIAS footprinting;
- accessibility-specific QC interpretation.

FAIRE-seq:
- independent assay contract;
- only reuse ATAC components when the scientific assumptions are genuinely
  shared;
- no automatic inheritance of ATAC-only TSS/footprinting semantics.

## 7. Configuration and site execution boundary

Nextflow configuration is split conceptually into:

1. portable workflow/pipeline parameters;
2. portable executor resource labels;
3. institution/site overlays for SGE/PBS/SLURM details, queue names, memory
   syntax, container cache/mirror paths, and filesystem policy.

The old run.sh scheduler auto-detection behavior is migration evidence, not the
new configuration architecture.

A pipeline should not embed one HPC site's queue or memory syntax in its
scientific configuration.

## 8. Testing and acceptance model

### 8.1 Custom DSL2 pipelines

Use:

- Nextflow stub runs for workflow composition;
- nf-test for processes/modules/subworkflows where appropriate;
- existing pytest/script tests for reused CLI scripts;
- scenario tests for feature behavior;
- selected interaction scenarios for coupled options.

For complex pipelines, especially chromatin routes, testing every flag in
isolation is not sufficient. Known interacting options receive explicit
combination scenarios.

### 8.2 Adopted nf-core pipelines

Use:

- upstream test profile to validate execution plumbing;
- BioWorkflows adapter/config validation;
- one or more project-owned real-data acceptance runs;
- pinned revision records.

### 8.3 Real-data cross-acceptance

The later per-pipeline specification defines objective metrics.

Examples of metric classes:

- RNA/small-RNA: mapping, quantification/count concordance, retained feature
  completeness and differential-analysis consistency;
- BS-seq: mapping, CpG coverage/methylation agreement and DMR consistency;
- ChIP/CUT&Tag: mapping/QC, peak-set and signal agreement, replicate/control
  behavior and optional spike-in behavior;
- ATAC/FAIRE: mapping/QC, accessibility peak/signal agreement, assay-specific
  QC, and optional footprinting outputs;
- seCLIP: deduplication, peak/control behavior and reproducible-peak outputs.

A changed upstream method may intentionally produce non-identical output. Such
changes require an explicit migration decision and acceptance rationale rather
than being hidden behind a compatibility threshold.

## 9. Migration phases

### Phase 0 - Migration foundation and contracts

Freeze and implement the common migration substrate:

- Nextflow repository skeleton;
- tested Nextflow/nf-core/nf-test version manifest;
- normalized run-manifest contract;
- pipeline-specific adapter pattern;
- executor/site-overlay configuration pattern;
- Apptainer/container baseline;
- custom DSL2 template;
- initial Nextflow CI;
- migration matrix and acceptance-record format.

The five legacy workflows stay unchanged except for necessary bug fixes.

### Phase 1 - Adoption pilot: BS-seq

Migrate the BS-seq core to pinned nf-core/methylseq and implement the small
dmr-nf downstream component.

This phase proves:

- adopted-pipeline execution;
- input/reference adapters;
- container/HPC profiles;
- artifact handoff;
- real-data acceptance.

### Phase 2 - Custom DSL2 pilot: small-RNA

Build srna-nf.

This phase proves the BioWorkflows-owned DSL2 template, local module policy,
validation pattern, custom container execution, and nf-test strategy.

The finalized custom pattern becomes the reference for later custom pipelines.

### Phase 3 - seCLIP custom migration

Build seclip-nf using the Phase 2 template and reusable nf-core modules where
appropriate. Preserve the in-house UMI/input-control/reproducible-peak
scientific contract.

### Phase 4 - Chromatin-family decomposition and migration

Migrate the legacy chip_cuttag_atac_faire capabilities as independent assay
units.

The phase begins by freezing the shared chromatin artifact boundaries, then
migrates and accepts the assay routes independently:

- CUT&Tag -> nf-core/cutandrun adoption route;
- ATAC-seq -> nf-core/atacseq adoption route;
- ChIP-seq -> nf-core/chipseq adoption route;
- FAIRE-seq -> custom faire-nf.

BioWorkflows-specific downstream capabilities are attached as bounded reusable
components where needed.

The legacy combined Snakemake workflow remains available until all four
replacement assay routes required for feature parity have passed their own
acceptance gates.

### Phase 5 - RNA-seq adoption and side pipelines

Adopt nf-core/rnaseq, establish the explicit handoff to
nf-core/differentialabundance where retained, and build lncrna-nf for the
project-specific lncRNA discovery path.

RNA-seq is scheduled after the two migration patterns and chromatin
decomposition have stabilized because it has a broad downstream surface.

### Phase 6 - Consolidation and legacy retirement

Per accepted execution unit:

- finalize user documentation and examples;
- record the accepted version/configuration baseline;
- mark the legacy implementation as superseded.

Before deleting a legacy implementation from the active tree:

- create an annotated Git tag/release for the final Snakemake state;
- confirm all replacement units required for that legacy workflow have passed
  acceptance;
- keep migration/acceptance records tracked.

There is no separate late "containerization sweep": container-native execution
is part of the new implementation from the beginning.

## 10. Mainline risks and controls

### Upstream nf-core drift

Control: pin revisions; upgrade through dedicated changes and re-acceptance.

### Non-model references

Control: BioWorkflows-owned reference adapters/configuration for rice and custom
genomes; do not depend on iGenomes availability.

### Over-generalization of chromatin assays

Control: assay-first top-level contracts and explicit ownership of
assay-specific QC/peak semantics.

### Duplication after assay split

Control: share only scientifically identical modules/subworkflows; accept some
assay-local duplication when semantics differ.

### Local custom code becoming a second nf-core fork

Control: keep adopted pipelines upstream-owned; put BioWorkflows-specific
capabilities behind explicit artifact boundaries instead of patching upstream
internals by default.

### Runtime inconsistency

Control: container-native target plus version manifest and tested HPC site
profiles.

### False parity

Control: objective per-pipeline acceptance specifications before retirement.

## 11. Explicitly deferred from this mainline document

The following belong to later implementation specifications and are not frozen
here:

- exact Nextflow/nf-core/nf-test version numbers;
- exact directory/file names beyond the ownership model;
- exact samplesheet schemas for each adapter;
- exact container images for in-house scripts;
- detailed per-process resource labels;
- exact acceptance thresholds;
- exact chromatin post-processing component boundaries;
- exact implementation order within the four chromatin assay routes when work
  can proceed independently.

These are implementation decisions constrained by this mainline, not missing
mainline decisions.

## 12. First implementation milestone after mainline freeze

Only after this document passes a final consistency review:

1. write the Phase 0 implementation specification;
2. establish the version and execution baselines;
3. create the normalized run-manifest and adapter contracts;
4. create the minimal custom DSL2 template and CI;
5. start the BS-seq adoption pilot.

No assay migration implementation should begin before Phase 0 contracts are
frozen.
