# BioWorkflows: Scientific-First nf-core-Based Nextflow Rebaseline

- Date: 2026-10-05
- Status: FROZEN after independent mainline review
- Scope: scientific and architectural rebaseline of the current Snakemake
  workflows into a simpler nf-core-based Nextflow workflow suite
- Primary use case: personal/research-group use on local Linux/WSL and HPC
- Legacy policy: existing Snakemake workflows remain available during the
  transition, but they are references and historical baselines, NOT migration
  compatibility targets
- Authority: this document freezes the project-level direction. Detailed
  scientific choices, exact pipeline revisions, parameters, and implementation
  specifications are decided per workflow after a current-state review.

## 1. Project positioning

BioWorkflows is a **personal modular workflow suite built on the nf-core
ecosystem**.

It is not intended to become:

- a second nf-core framework;
- a generic workflow platform;
- a compatibility layer around the old Snakemake repository;
- a framework that forces all assays through one local abstraction.

The project should use as much maintained community infrastructure as possible
and keep BioWorkflows-owned code focused on real scientific or operational gaps.

The expected value of BioWorkflows is:

- selecting and pinning scientifically appropriate upstream workflows;
- keeping useful personal/reference/HPC configuration;
- adding analysis that is genuinely missing upstream;
- composing a small number of real downstream handoffs;
- validating workflows on the user's actual biological use cases;
- retaining custom analysis only when it remains scientifically justified.

## 2. The migration is an upgrade, not a reproduction exercise

The central rule of this rebaseline is:

> **Scientific quality and useful functionality take precedence over legacy
> compatibility, implementation preservation, output-layout preservation, and
> numerical reproduction of the old workflow.**

The old Snakemake implementation provides evidence about:

- the data types the user analyzes;
- previously useful features;
- known edge cases;
- historical outputs that can help detect regressions;
- locally validated operating knowledge.

It does NOT define what the new workflow must contain.

During rebaseline, an old feature may be:

1. **adopted** from a maintained upstream workflow;
2. **replaced** by a scientifically better current method;
3. **expanded** using richer upstream/community functionality;
4. **retained locally** because it still fills a real gap;
5. **dropped** because it is obsolete, redundant, weakly justified, or no
   longer useful.

Feature parity with the legacy implementation is therefore NOT a freeze or
retirement requirement.

## 3. Scientific-first decision order

For every sequencing/assay workflow, decisions are made in this order.

### Step 1 - Define the biological objective

State what biological question and output the workflow should support today.

Do not start from the old rule graph.

### Step 2 - Review the current ecosystem

Review, at implementation time:

- current stable nf-core pipelines;
- relevant nf-core modules/subworkflows;
- maintained community workflows;
- current commonly accepted analysis methods;
- known assay-specific guidance and limitations.

The review must be refreshed before implementation because upstream pipelines
and best practices change.

### Step 3 - Choose the best route

Use the simplest scientifically adequate route:

1. direct upstream adoption;
2. upstream adoption plus a bounded local extension;
3. composition of community/nf-core components;
4. custom BioWorkflows DSL2 only when a genuine gap remains.

### Step 4 - Choose methods independently of the legacy implementation

Do not retain:

- a tool merely because the Snakemake workflow used it;
- an old parameter merely because it was previously validated;
- a downstream step merely to preserve historical output;
- a local script when a maintained community implementation is better.

### Step 5 - Validate the new scientific workflow

Validation asks whether the new workflow is scientifically credible, complete,
reproducible, and useful.

Comparison with the legacy workflow is supporting evidence, not the acceptance
target.

## 4. Frozen architecture principles

### D1. nf-core first, community-aware

Use maintained nf-core pipelines and components as the default ecosystem.

When nf-core does not cover a requirement well, also consider maintained
community workflows before deciding to implement a custom pipeline.

Custom code is the last option, not the default destination of old custom code.

### D2. One repository, flat workflow entries

BioWorkflows remains one repository.

Each sequencing/assay type is an independent first-class workflow entry. There
is no formal RNA-family, chromatin-family, or other family execution layer.

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
|-- optional downstream analyses
`-- future workflow entries
```

Documentation may group related assays for navigation only.

### D3. No giant multi-assay workflow

The legacy `chip_cuttag_atac_faire` architecture is explicitly retired as a
design pattern.

ChIP-seq, CUT&Tag, ATAC-seq, and FAIRE-seq receive separate workflow entry
points and separate scientific reviews.

Shared operations may use common modules/subworkflows only when their scientific
semantics are genuinely the same.

### D4. Native upstream interfaces are preferred

For an adopted upstream pipeline, prefer its native:

- samplesheet/input schema;
- parameters;
- outputs;
- containers;
- release/version model;
- testing behavior.

BioWorkflows does not require a universal samplesheet, universal run manifest,
universal artifact model, or legacy output compatibility layer.

Small adapters may be written only when they solve repeated real work.

### D5. Cross-workflow contracts are local and evidence-driven

Do not define a global artifact framework.

Define a handoff only when one workflow actually consumes another workflow's
outputs.

Examples that may justify explicit handoffs include:

- rnaseq -> differential analysis;
- rnaseq -> lncRNA analysis;
- methylation calling -> differential methylation.

The exact downstream method is not frozen by this mainline.

### D6. Container-native execution

Adopted nf-core workflows use their native container support.

BioWorkflows-owned Nextflow processes should use nf-core module software
definitions or explicit per-process containers.

The legacy WorkflowSpec executable resolver is not part of the future runtime
architecture.

A developer Conda profile may exist, but HPC/reproducible use is
container-oriented.

### D7. Explicit version locking, without freezing versions in the mainline

Production and scientific validation runs must use explicit tested revisions.

A small tracked baseline records at least:

- Nextflow version;
- adopted workflow revisions;
- nf-core tooling where relevant;
- nf-test where used;
- local module/container revisions where needed.

This mainline does not freeze exact versions because those must be selected from
the current ecosystem at implementation time.

### D8. Minimal local shared layer

Only share things with demonstrated cross-workflow value, for example:

- HPC/site profiles;
- Apptainer/cache/mirror configuration;
- version pins;
- reusable reference configuration where appropriate;
- genuinely reusable modules/subworkflows;
- focused testing conventions.

Do not build during this migration:

- a central orchestration service;
- a workspace/project object model;
- a plugin registry;
- a global capability registry;
- a global metadata database;
- a universal workflow schema;
- a mandatory BioWorkflows CLI;
- a custom replacement for nf-core tooling.

### D9. Scientific acceptance, not parity acceptance

A workflow is accepted when it is:

1. scientifically appropriate for its stated biological objective;
2. functionally complete for the intended use;
3. reproducible;
4. successfully validated on relevant real data;
5. operationally usable on the target environment.

Legacy numerical/output similarity may be checked when informative, but:

- exact parity is not required;
- changed tools/defaults are allowed;
- richer new functionality is allowed;
- obsolete legacy functionality may be removed.

Any important scientific difference should be understood and documented, not
necessarily eliminated.

### D10. Git-native legacy preservation

Before removing a legacy implementation from the active tree, preserve its
final state with an annotated Git tag/release.

This is historical reproducibility only. It does not impose compatibility on
the replacement.

## 5. Independent ecosystem review findings (2026-10-05)

These observations informed the freeze decision. They are evidence, not pinned
future versions.

### 5.1 RNA-seq

Current nf-core/rnaseq is substantially richer and more actively maintained than
the legacy local RNA-seq implementation and supports multiple alignment and
quantification strategies plus extensive QC.

Mainline conclusion:

- nf-core/rnaseq is the primary candidate;
- do not preserve old STAR/featureCounts/StringTie/parameter choices merely for
  parity;
- downstream differential and lncRNA analysis must be separately re-evaluated.

### 5.2 Differential analysis

nf-core/differentialabundance provides a maintained matrix-based framework for
differential statistics, plots, gene-set analysis, reports, and multiple
analysis profiles.

Mainline conclusion:

- evaluate it as the default downstream analysis route;
- retain local DE/enrichment scripts only where they offer a scientifically
  necessary capability not covered adequately upstream.

### 5.3 DNA methylation

Current nf-core/methylseq supports several methylation analysis paths rather than
only the legacy Bismark path, with modern QC and additional capabilities.

Mainline conclusion:

- nf-core/methylseq is the primary upstream candidate;
- Bismark is not preserved as a requirement;
- the appropriate aligner/caller route is selected during scientific rebaseline;
- differential methylation is reviewed separately;
- methylKit is not frozen as the required DMR method.

### 5.4 Small RNA

Current nf-core/smrnaseq includes considerably more functionality than assumed
in the earlier migration draft, including UMI handling, multiple contamination
classes, miRNA/isomiR analysis, genome quantification, and novel miRNA analysis.

However, some upstream documentation explicitly notes limited validation of
contamination filtering outside human data.

Mainline conclusion:

- do NOT pre-classify small-RNA as a custom workflow;
- first evaluate nf-core/smrnaseq on the actual plant/rice requirements;
- use direct adoption or a bounded extension if scientifically adequate;
- build a custom plant-oriented workflow only for proven gaps.

### 5.5 ChIP-seq

nf-core/chipseq provides a maintained ChIP-seq-specific workflow with peak
calling, QC, and differential-analysis capabilities.

Mainline conclusion:

- direct adoption is the default starting point;
- legacy ChIP-specific code is retained only if a current capability gap is
  demonstrated.

### 5.6 CUT&Tag

nf-core/cutandrun explicitly supports CUT&Tag and includes spike-in support,
controls, SEACR/MACS2 peak calling, consensus processing, and downstream QC.

Mainline conclusion:

- direct adoption is the default starting point;
- legacy CUT&Tag branches from the mixed workflow are not migration targets.

### 5.7 ATAC-seq

nf-core/atacseq is a dedicated maintained ATAC-seq workflow.

Mainline conclusion:

- direct adoption is the default starting point;
- compare its current QC/peak/accessibility functionality with actual research
  needs before creating any local extension;
- old TOBIAS or other optional steps are retained only if still scientifically
  useful and not already better covered by current community tools.

### 5.8 CLIP/eCLIP

The released nf-core/clipseq pipeline is old DSL1 and is incompatible with
modern Nextflow releases. It therefore cannot simply be adopted as the modern
replacement.

At the same time, maintained community DSL2 CLIP workflows exist and should be
reviewed before rebuilding the legacy local implementation.

Mainline conclusion:

- perform a fresh CLIP/eCLIP ecosystem review;
- do not port the old seCLIP pipeline by default;
- select a maintained community route or compose modern modules where possible;
- write a new local DSL2 workflow only for remaining scientifically necessary
  gaps.

### 5.9 lncRNA

nf-core/lncpipe is under active modernization but remains under development.

Mainline conclusion:

- do not freeze the legacy CNCI/Pfam/NR chain as the target design;
- review current lncRNA identification/coding-potential/community methods when
  this workflow is reached;
- adopt, compose, or implement only the scientifically justified final route.

### 5.10 FAIRE-seq

No mature nf-core FAIRE-seq replacement is assumed by this mainline.

Mainline conclusion:

- perform a current community review before implementation;
- reuse compatible chromatin/accessibility components where scientifically
  valid;
- create a small custom workflow only if no maintained alternative fits.

## 6. Reframing the legacy workflow inventory

The old workflow names no longer determine the new implementation route.

The current working disposition is:

| Biological workflow | Primary candidate | Mainline status |
|---|---|---|
| RNA-seq | nf-core/rnaseq | adopt-first; scientific rebaseline required |
| Differential abundance | nf-core/differentialabundance | adopt-first; only extend for proven gaps |
| BS-seq / methylation calling | nf-core/methylseq | adopt-first; method choice re-evaluated |
| Differential methylation | current community methods | open; do not freeze methylKit |
| small-RNA | nf-core/smrnaseq | re-evaluate first; custom only if plant gaps remain |
| ChIP-seq | nf-core/chipseq | adopt-first |
| CUT&Tag | nf-core/cutandrun | adopt-first |
| ATAC-seq | nf-core/atacseq | adopt-first |
| FAIRE-seq | community review / custom if needed | open |
| seCLIP/eCLIP | community review / composition / custom if needed | open |
| lncRNA analysis | current nf-core/community methods | open; do not preserve legacy chain automatically |

This table is a routing baseline, not a promise of feature parity with the old
repository.

## 7. Decomposition of chip_cuttag_atac_faire

The combined legacy workflow is split because the assays have different
scientific semantics and now have stronger dedicated ecosystem support.

```
legacy chip_cuttag_atac_faire
        |
        +--> ChIP-seq  --> review/adopt nf-core/chipseq
        +--> CUT&Tag   --> review/adopt nf-core/cutandrun
        +--> ATAC-seq  --> review/adopt nf-core/atacseq
        `--> FAIRE-seq --> current ecosystem review; custom only if needed
```

The new design does not preserve:

- the shared `seqtype` switch;
- a single mixed-assay samplesheet;
- global flags that change meaning by assay;
- the old feature matrix merely for compatibility.

If an upstream workflow has a more scientifically appropriate QC, peak-calling,
normalization, replicate, or reporting design, the upstream/current design may
replace the old local implementation.

A project can still contain several assay types; they simply run as independent
assay workflows.

## 8. Simplified repository direction

Exact paths are implementation details, but the repository should remain
conceptually simple:

```
BioWorkflows/
|-- legacy Snakemake workflows/        # temporary during transition
|-- nextflow/
|   |-- conf/                          # local/HPC/container configuration
|   |-- versions.yml                   # tested revisions
|   |-- adopted/                       # thin configs/docs/recipes for upstream workflows
|   |-- pipelines/                     # only genuinely necessary local DSL2 workflows
|   |-- modules/                       # local modules only when upstream lacks one
|   |-- subworkflows/                  # only demonstrated reusable compositions
|   `-- docs/                          # scientific decisions + validation records
|-- plans/
`-- README.md
```

An `adopted/` entry is not a fork. It may contain only:

- a params file;
- rice/reference configuration;
- an HPC profile/example command;
- a short scientific decision record;
- validation notes.

If upstream can be run cleanly without even that directory, do not create one.

## 9. Per-workflow scientific rebaseline record

Before implementing or adopting each workflow, create one concise scientific
decision record answering:

1. What biological use case is required?
2. What inputs and scientifically meaningful outputs are needed?
3. What does the current nf-core/community ecosystem provide?
4. Which current methods are preferable and why?
5. Which legacy functions are:
   - replaced,
   - retained,
   - expanded,
   - dropped?
6. Is any local code still necessary?
7. What real-data validation demonstrates scientific usability?

This record comes before implementation specifications.

It prevents both blind migration and uncontrolled feature accumulation.

## 10. Testing and validation philosophy

### 10.1 Adopted workflows

Do not duplicate upstream module-level tests.

BioWorkflows validation focuses on:

- local/HPC execution;
- reference configuration;
- real biological datasets;
- expected scientific outputs;
- chosen optional functionality;
- important upstream version changes.

### 10.2 Local workflows/extensions

Use nf-test or appropriate focused tests for locally owned logic.

Reuse upstream-tested modules rather than reproducing their tests locally.

### 10.3 Legacy comparison

Legacy comparison is optional and targeted.

Use it when it helps answer questions such as:

- did a major biological signal disappear unexpectedly?
- is a new default producing a surprising systematic difference?
- did a previously important capability disappear?

Do not fail a new workflow because:

- file names changed;
- output directories changed;
- tools changed;
- defaults changed;
- numerical values differ for understood scientific reasons;
- a weak/obsolete legacy feature was intentionally removed.

## 11. Migration sequence

### Stage 0 - Minimal shared foundation

Implement only the common infrastructure that is already clearly necessary:

- Nextflow baseline;
- Apptainer/HPC execution;
- tested version recording;
- thin adopted-workflow convention;
- local custom-pipeline convention;
- minimal CI;
- concise scientific decision/validation record templates.

Do NOT build a universal manifest, registry, CLI, artifact platform, or
orchestration layer.

### Stage 1 - First adoption case

Use one high-confidence upstream workflow as the first real migration.

BS-seq/nf-core/methylseq remains a good candidate, but the selected internal
methylation route is chosen scientifically rather than inherited from legacy.

The purpose is to validate the thin adoption model, not to reproduce bs-seq
v0.2.

### Stage 2 - Continue by research value and scientific confidence

After the first adoption case, migrate workflows based on:

- current research need;
- confidence in the upstream/community route;
- expected scientific benefit;
- maintenance reduction.

There is no mandatory "custom pipeline pilot".

If small-RNA can be solved well with nf-core/smrnaseq plus a small plant
extension, do that instead of creating srnaseq solely to demonstrate a custom
pipeline.

### Stage 3 - Build local DSL2 only where the review proves it is necessary

Likely candidates may include parts of CLIP/eCLIP, FAIRE-seq, lncRNA, or
differential methylation, but none is pre-committed to custom implementation.

### Stage 4 - Legacy retirement

Retire old code when the user's real analysis needs are covered by accepted new
routes.

Retirement does NOT require every historical option to have a replacement.

Before deletion:

- document deliberately dropped capabilities;
- preserve the final legacy state by Git tag/release;
- confirm that currently required biological use cases are covered.

## 12. Future-workflow rule

For any new assay or analysis:

```
Define biological need
        |
Review current nf-core + maintained community workflows
        |
        +--> suitable upstream exists
        |       |
        |       +--> use directly
        |       `--> add a small extension only if needed
        |
        `--> no suitable upstream
                |
                +--> compose existing modules/subworkflows
                |
                `--> custom DSL2 only for the remaining gap
```

The project should prefer deleting local code over owning equivalent code when a
maintained upstream solution becomes better.

## 13. Freeze-review conclusions

Independent review identified and corrected five remaining migration biases in
the previous draft:

1. **legacy parity bias** - removed;
2. **premature custom small-RNA decision** - removed;
3. **premature methylKit / legacy lncRNA method retention** - removed;
4. **mandatory custom-pipeline pilot** - removed;
5. **legacy output/scientific consistency as acceptance target** - replaced by
   scientific usability and current best-practice validation.

The mainline is therefore frozen with the following interpretation:

> BioWorkflows will evolve by adopting and composing the best maintained
> nf-core/community workflows for the user's real scientific needs, while
> keeping only the minimum amount of local workflow code required for genuine
> gaps.

## 14. What is explicitly deferred

The frozen mainline does NOT decide:

- exact nf-core pipeline versions;
- exact tool/aligner/caller choices inside adopted pipelines;
- exact parameters;
- exact DMR method;
- exact lncRNA identification strategy;
- exact small-RNA route;
- exact CLIP/eCLIP route;
- exact FAIRE implementation;
- exact local module boundaries;
- exact validation thresholds;
- implementation order after the first adoption case.

Those decisions require a fresh scientific review at the time each workflow is
implemented.
