# BioWorkflows: Nextflow Migration Implementation Guidance

- Date: 2026-10-05
- Status: GUIDANCE
- Depends on: `plans/2026-10-05-nextflow-migration-design.md`
- Purpose: provide a practical implementation arrangement for the frozen
  scientific-first nf-core-based mainline without creating rigid per-stage
  implementation specifications.
- Scope: migration and optimization of the current BioWorkflows analyses.
- Rule: implementation remains adaptable to the current nf-core/community
  ecosystem and to scientific evidence found during each workflow review.

## 1. How to use this document

This document is intentionally lighter than a traditional implementation
specification.

It defines:

- the preferred implementation direction for each current analysis workflow;
- the minimal shared work that should be done once;
- the order in which implementation decisions should be made;
- what should remain flexible until real implementation begins;
- what should NOT be preserved merely for compatibility.

It does NOT freeze:

- exact Nextflow process graphs;
- exact tool versions or parameters;
- exact nf-core revisions;
- exact file layouts;
- exact acceptance thresholds;
- detailed task-by-task schedules.

The frozen mainline remains authoritative for architecture. This document is a
working guide for implementation.

## 2. Global implementation rule

The migration is an optimization and upgrade, not a translation exercise.

For every workflow:

1. start from the current biological objective;
2. review the best maintained nf-core/community options;
3. adopt scientifically better functionality even when it differs from the
   legacy Snakemake implementation;
4. retain local code only when it still provides useful capability;
5. add new upstream/community functionality when it improves analysis quality
   or practical usability;
6. do not preserve old tools, parameters, output paths, or optional stages only
   for compatibility.

Legacy results are useful for troubleshooting and comparison, but are not the
design target.

## 3. Minimal common foundation

Only a small shared Nextflow foundation is required.

Implement once:

- a tested Nextflow baseline;
- Apptainer/Singularity execution for local/HPC use: the apptainer profile is
  the only supported execution profile on the HPC; image seeding follows the
  validated per-source routing recipe in section 3.1; machine-specific paths
  stay in machine-local git-ignored configuration;
- SGE/PBS/SLURM site configuration where actually needed;
- a tracked record of tested workflow/tool revisions;
- a simple convention for BioWorkflows-owned custom DSL2 pipelines;
- a simple convention for thin adopted-workflow configuration;
- minimal CI for locally owned code/configuration, scoped to static checks
  only (YAML/shell lint, Nextflow config syntax validation, documentation
  link checks); upstream test profiles are not re-run locally, and nf-test is
  added when the first locally owned module or pipeline appears;
- concise scientific decision and real-data validation records.

Avoid introducing:

- a universal samplesheet;
- a global run manifest;
- a global artifact model;
- a workflow-family hierarchy;
- a mandatory wrapper CLI;
- a plugin/registry/workspace framework;
- compatibility output layers.

Add shared utilities only after repeated real use demonstrates that they are
worth sharing.

### 3.1 Image acquisition recipe (validated 2026-10-05 on the HPC)

For each adopted pipeline, the seeding procedure is:

1. Pin the pipeline revision (branch or commit SHA) compatible with the
   pinned Nextflow driver. Tagged releases may lag the active branch by a
   year or more and may use config syntax rejected by modern Nextflow
   (chipseq 2.1.0 is the recorded example).
2. Resolve the authoritative image list with `nextflow inspect -format json
   -profile <execution profiles> main.nf`. Never grep module files: unused
   module variants and blob-path style references pollute the result.
3. Pull images sequentially, one at a time, routed per source: direct
   connections where reachable and stable (wave, depot, quay), the layered
   proxy environment where required (shell/Go read https_proxy; the JVM only
   reads NXF_JVM_ARGS system properties), and - for images unreachable by
   either - local builds from the USTC bioconda mirror using the
   version/build string encoded in the image tag (a fresh conda env, not the
   base env, to avoid python freeze conflicts).
4. Where a source is rewritten (docker.io-hosted biocontainers images pulled
   from quay.io/biocontainers - DaoCloud rejects biocontainers by allowlist)
   or built locally, record a machine-local container-override configuration
   mapping process names to the cached image URI or local SIF path, and pass
   it with `nextflow -c`.
5. Validate with the pipeline's test profile plus the apptainer profile on
   the HPC, throttling the local executor on the shared head node (full
   parallelism hits fork EAGAIN alongside lab jobs).
6. Record the tested driver, pipeline revision, and image set in the version
   baseline (BASELINE.yml machine-local today; the tracked versions.yml when
   the nextflow/ tree exists).

The machine-local tooling on ginpie (`nf-pull-seq.sh`, `nf-pull-env.sh`,
`nf-build-depot.sh`) implements this recipe and is the reference
implementation for seeding future pipelines.

## 4. Workflow implementation directions

### 4.1 RNA-seq — community-first adoption

Primary direction:

- use mature community/nf-core RNA-seq workflows as the main execution path;
- nf-core/rnaseq is the default candidate;
- use upstream-supported aligners, quantifiers, QC, reporting, and optional
  functionality based on scientific need rather than preserving the legacy
  STAR/featureCounts/StringTie configuration.

BioWorkflows should own only:

- preferred parameter/reference configuration;
- rice/non-default reference handling where needed;
- HPC execution configuration;
- missing downstream functionality that is genuinely useful;
- validation on representative real datasets.

Do not rebuild standard RNA-seq execution locally.

### 4.2 lncRNA — retain local foundation, then scientifically upgrade

Primary direction:

- use the existing BioWorkflows lncRNA analysis as a functional starting point;
- do not freeze the old CNCI/Pfam/NR chain as the final method;
- review current coding-potential, transcript-quality, annotation, homology,
  domain, and other relevant community approaches;
- replace or expand legacy components where better maintained/scientifically
  stronger methods exist;
- reuse outputs from the adopted RNA-seq path where appropriate.

The target is therefore a BioWorkflows-owned lncRNA workflow optimized from the
existing implementation, not a compatibility port.

Local code should remain only where it adds value after current-method review.

### 4.3 small-RNA — retain local plant-oriented foundation, then optimize

Primary direction:

- use the current BioWorkflows small-RNA workflow as the scientific starting
  point because it contains plant/rice-oriented cascade filtering and counting
  logic that reflects actual use;
- review nf-core/smrnaseq and other maintained community workflows for useful
  capabilities;
- absorb better modules, QC, UMI handling, contamination handling, miRNA/isomiR
  analysis, novel miRNA analysis, reporting, and other useful features where
  scientifically appropriate;
- remove local implementations that become redundant;
- preserve plant-specific behavior only where it remains useful and supported by
  evidence.

The target is a scientifically upgraded BioWorkflows small-RNA workflow, not a
strict port of the existing Snakemake DAG and not an automatic full replacement
by nf-core/smrnaseq.

### 4.4 BS-seq / DNA methylation — community replacement route

Primary direction:

- replace the legacy local core workflow with a stronger maintained
  community/nf-core methylation workflow;
- nf-core/methylseq is the default candidate, but the exact aligner/caller path
  should be selected based on current scientific suitability;
- do not preserve Bismark solely because the old workflow used it;
- evaluate modern QC, extraction/calling, reporting, and reference handling from
  the selected upstream workflow.

Differential methylation:

- review current DMR/DMC methods independently;
- methylKit may be retained if still appropriate, but it is not the default by
  legacy inheritance;
- use a local downstream component only when the selected community workflow
  does not adequately cover the required analysis.

This is a replacement-and-upgrade route, not a port.

### 4.5 ChIP-seq — mature community workflow first

Primary direction:

- adopt a mature community/nf-core ChIP-seq workflow as the main execution path;
- nf-core/chipseq is the default candidate;
- use its maintained alignment, QC, peak calling, replicate/control handling,
  reporting, and other current functionality where scientifically appropriate.

BioWorkflows only supplements missing capabilities that are actually needed.

Examples of retained/local additions are allowed only after checking whether
the selected upstream workflow already provides an equal or better route.

Do not recreate the old mixed `chip_cuttag_atac_faire` architecture.

### 4.6 CUT&Tag — mature community workflow first

Primary direction:

- use a mature community workflow, with nf-core/cutandrun as the default
  candidate for CUT&Tag;
- prefer its maintained treatment/control, spike-in, peak calling, consensus,
  QC, and reporting functionality where appropriate;
- add local functionality only for demonstrated gaps.

Do not preserve the old CUT&Tag branch merely because it was previously
validated locally.

### 4.7 ATAC-seq — mature community workflow first

Primary direction:

- use a mature dedicated ATAC-seq workflow, with nf-core/atacseq as the default
  candidate;
- adopt stronger community QC/accessibility/peak/reporting functionality when
  available;
- review whether optional legacy analyses such as TOBIAS footprinting remain
  useful and whether the community ecosystem provides a better implementation.

Only missing, scientifically useful analysis should remain BioWorkflows-owned.

### 4.8 FAIRE-seq — current ecosystem review, then minimal custom work if needed

Primary direction:

- do not inherit ATAC behavior automatically;
- review maintained community workflows and compatible modules at implementation
  time;
- if no suitable maintained workflow exists, build a small FAIRE-specific DSL2
  workflow using reusable community components where scientifically valid.

Keep this workflow independent from ATAC-seq even when selected modules overlap.

### 4.9 seCLIP/eCLIP — current ecosystem review before deciding ownership

Primary direction:

- do not automatically port the current local seCLIP pipeline;
- review maintained modern CLIP/eCLIP community workflows and available modules;
- compare them against the actual requirements: UMI handling, controls,
  reproducible peaks, annotation, QC, and downstream use;
- adopt or compose community solutions if they are stronger;
- retain/build local DSL2 only for genuine missing capability.

The existing pipeline remains a source of tested requirements and edge cases,
not the required future implementation.

## 5. Decomposition of the legacy mixed chromatin workflow

`chip_cuttag_atac_faire` is not migrated as one workflow.

Its roles become separate entry points:

```
ChIP-seq  -> community-first ChIP workflow + optional local extensions
CUT&Tag   -> community-first CUT&Tag workflow + optional local extensions
ATAC-seq  -> community-first ATAC workflow + optional local extensions
FAIRE-seq -> community review, then minimal custom workflow if required
```

Shared code is extracted only when the same implementation and scientific
semantics are genuinely reusable.

A shared peak utility or reporting component may be reasonable.

A shared "chromatin workflow" layer is not.

## 6. Suggested implementation arrangement

The work does not need rigid numbered implementation specifications. A practical
arrangement is:

### Foundation work

Establish the small common Nextflow/HPC/container/version baseline first.

### Community-adoption work

Prefer early migration of workflows with strong maintained upstream coverage:

- BS-seq / methylation;
- RNA-seq;
- ChIP-seq;
- CUT&Tag;
- ATAC-seq.

These should usually require less locally owned workflow code and therefore
provide the fastest reduction in maintenance burden.

### Local optimization work

Develop the workflows where BioWorkflows-specific scientific logic remains
valuable:

- small-RNA;
- lncRNA.

These should be treated as upgrades of the current local scientific workflows,
while actively replacing local pieces with stronger community components.

### Open-review work

Handle workflows whose best route is not yet clear:

- seCLIP/eCLIP;
- FAIRE-seq;
- differential methylation.

For these, perform the ecosystem/scientific review immediately before
implementation rather than freezing a route now.

The exact order inside and between these groups may follow current research
priority.

## 7. Workflow-level implementation record

Before changing one workflow, create only a concise decision note, not a full
implementation specification.

It should answer:

1. What biological use case are we supporting now?
2. What current upstream/community workflow(s) were reviewed?
3. Which route was selected and why?
4. Which legacy functions will be retained, replaced, expanded, or dropped?
5. What local code is still necessary?
6. What representative real data will be used for validation?
7. What would make the new workflow scientifically acceptable?

This note can evolve during implementation as evidence is collected.

## 8. Validation guidance

Validation should focus on whether the new workflow is scientifically and
practically better for current use.

Evaluate, as relevant:

- data/QC completeness;
- expected biological signal;
- robustness of alignment/quantification/peak/methylation outputs;
- replicate/control behavior;
- reference/species compatibility;
- useful reporting;
- reproducibility on the target HPC;
- performance and maintainability where meaningful.

Legacy comparison may be used diagnostically, but no workflow should be forced
back to an inferior method to improve numerical parity.

The first concrete validation record followed exactly this shape on
2026-10-05: an nf-core/chipseq dev-HEAD smoke run on the HPC with the test
and apptainer profiles, a throttled local executor, verified outputs
(MultiQC report, peak calls, bigWigs), and a machine-local baseline record
(BASELINE.yml). It is the template for per-workflow validation records.

## 9. Retirement guidance

A legacy workflow or feature can be retired when:

- current real research use cases are covered;
- the replacement route has passed scientific/operational validation;
- deliberately dropped functionality has been documented;
- the final legacy state has been preserved with a Git tag/release.

It is NOT necessary to implement every historical option before retirement.

For `chip_cuttag_atac_faire`, the combined workflow can be retired after the
currently required ChIP/CUT&Tag/ATAC/FAIRE use cases have adequate independent
replacement routes; historical unused options do not require one-to-one
replacement.

## 10. Maintenance rule after migration

Prefer upstream maintenance over local ownership.

When a community/nf-core workflow or module later becomes clearly better than a
BioWorkflows-local implementation:

- evaluate it;
- migrate to it when practical;
- delete the redundant local code.

BioWorkflows should remain thin enough that its unique code corresponds mainly
to unique scientific needs.
