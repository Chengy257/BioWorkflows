# BioWorkflows — Agent-Ready Bioinformatics Infrastructure
## Architecture Exploration & Unified Discussion Baseline

- **Date:** 2026-10-31
- **Status:** **WORKING / DISCUSSION / NOT FROZEN**
- **Repository:** Chengy257/BioWorkflows
- **Discussion branch:** docs/agent-ready-infrastructure-discussion-20261031
- **Repository baseline at discussion entry:** main @ 8839e877a451b51eca460f937a07b34dbe65b107
- **Document role:** Capture the agreed direction, current repository reality, candidate architectural models, non-goals, risks, and unresolved decisions. This is **not** an implementation specification or authorization to refactor code.
- **Change policy:** Amend this document as discussion advances; keep competing proposals and rationale visible. No architecture is frozen solely by inclusion here.

---

## 0. Executive position

BioWorkflows is proposed to evolve from a personal Nextflow/nf-core-based modular workflow suite into an **agent-ready infrastructure for general bioinformatics analysis**.

The intended distinguishing capability is **not** unrestricted LLM generation of Nextflow code. It is the ability for human researchers and general-purpose AI agents to **discover, understand, select, compose, execute, inspect, reproduce, and scientifically evaluate** bioinformatics capabilities through validated, auditable interfaces.

**Working positioning statement (candidate, not frozen):**

> BioWorkflows provides a reusable, scientifically grounded and reproducible bioinformatics capability infrastructure through which humans and general-purpose AI agents can discover, compose, execute, and inspect validated analyses, using Nextflow and community workflows as the computational foundation.

Guiding tension: **analysis composition should remain extensible and flexible; accepted execution plans and their evidence must remain explicit, controlled, and traceable**.

BioWorkflows is a scientific execution/capability infrastructure, **not itself a universal agent framework**. It is expected to work with multiple agents, models, operating environments, and scientific domains without claiming that any model can independently guarantee biological correctness.

## 1. Current repository reality and retained assets

The following is a **repository/docs-based inventory**, not a new end-to-end runtime acceptance:

| Existing surface | Current route | Evidence and status boundary | Reuse direction |
|---|---|---|---|
| RNA-seq | adopted nf-core/rnaseq + rnaseq-ext | Upstream test-profile HPC smoke; local extension wiring/unit/stub evidence, but some real-container/real-data follow-up deferred | Retain as initial complete-pipeline and downstream handoff candidate |
| Differential analysis | adopted nf-core/differentialabundance | Existing adoption/configuration; validation notes explicitly identify missing validated end-to-end DA run at that point | Evaluate before treating as production-capable |
| ChIP-seq / CUT&Tag | adopted nf-core/chipseq + chipseq-ext | Real-container fixture and PBS execution recorded for extension stages; real biological case remains distinct | Reuse as multi-stage scientific case |
| ATAC-seq / FAIRE-seq | adopted nf-core/atacseq + atacseq-ext | Configuration/unit/stub evidence; full real-container extension validation deferred | Reuse after qualified execution tests |
| lncRNA identification | local DSL2 lncrna-nf | Real-container PBS fixture chain and planted positive/negative scientific gates recorded 2026-10-26; real biological data validation is separate | Reuse as custom-pipeline composition case |
| small RNA, methylation, CLIP | historical Snakemake routes retired; no active Nextflow replacement | Not yet migrated; historical archive/tag retained | Future candidates, not current platform claims |

Other assets:
- nextflow/adopted/: thin wrappers/launch examples around community pipelines.
- nextflow/pipelines/: genuine local extension pipelines, process modules and subworkflows.
- nextflow/docs/: assay-specific scientific decision records and validation histories.
- nextflow/versions.yml: tested-revision ledger, not yet a run-specific immutable lock record.
- nextflow/conf/: common execution configuration, with site/machine paths outside Git.
- nextflow/tests/: local static/unit/stub verification, with hosted GitHub Actions intentionally retired.

This project should **preserve functioning scientific code and verified evidence**, not rebuild it merely to fit an abstract platform design.

**Audit caveat:** published documentation can lag source behavior; before implementation, verify every production claim against current code, executable tests, and real data.

## 2. Relationship to the 2026-10-05 migration baseline

The document plans/2026-10-05-nextflow-migration-design.md was frozen for the **earlier scientific-first migration**, not for the newly proposed infrastructure scope. It explicitly rejected a generic workflow platform, global capability registry, mandatory CLI, global artifact model, and universal run manifest (particularly D4/D5/D8).

That earlier policy conflicts with the current proposed objective. The contradiction is **acknowledged**, not silently ignored.

**Proposed authority handling (discussion only):**
1. Preserve the October migration document as an authoritative historical record for decisions made during migration and for scientific-first reuse policy.
2. Establish a **new, separately reviewed architecture authority** before implementing platform-wide functionality.
3. Carry forward compatible principles: nf-core/community-first, local code only for genuine gaps, native upstream interfaces when sensible, container-native HPC, explicit version pinning, real scientific validation, and Git-native history.
4. Reopen only necessary former prohibitions: common capability metadata, minimal cross-workflow contracts, plan/run identities, a unified human/agent interface, and provenance.
5. Do **not** interpret this working baseline as automatically superseding the earlier frozen document. Formal authority and migration policy remain an explicit decision for a subsequent architectural closeout.

A new common *envelope* should not force every upstream nf-core pipeline into a universal samplesheet or erase its native parameters. The platform must avoid replacing nf-core/Nextflow with an in-house workflow engine.

## 3. Candidate system boundary

### 3.1 In scope

- Discovery and documentation of reusable scientific capabilities.
- Scientific capability declarations, requirements, versioning, validation status, and compatibility checks.
- Composition at the level of complete pipelines, subworkflows, and individual modules.
- Plans built interactively from scientific intent, including method choices and constraints.
- Deterministic validation/compilation of accepted plans to executable Nextflow assets.
- Launch/monitor/resume/cancel/inspect where the chosen environment supports these operations.
- Immutable or append-only run/attempt, parameter, reference, artifact, and provenance evidence.
- Scientific and operational QC gates, explicit limitations, and reportable failed/abstained states.
- Human-first CLI/API and adapter interfaces for general agents (skills/MCP considered), independent of any one LLM.
- Local Linux/HPC first; portability to additional execution environments is a design consideration, not an initial implementation promise.

### 3.2 Out of scope / anti-goals (working)

- A new scheduler/workflow execution engine to replace Nextflow.
- Rewriting mature nf-core pipelines as bespoke BioWorkflows implementations.
- A required proprietary AI model, agent harness, hosted SaaS, or autonomous multi-agent runtime.
- Unrestricted shell execution or unrestricted code generation as the sole agent interface.
- Automatic acceptance of arbitrary scientific conclusions simply because a computational workflow succeeded.
- Global ontology/database/universal schema encompassing every biological concept in v1.
- Mandatory integration of other Chengy257 repositories or absorbing their domain-specific logic.
- Automatic GitHub Actions or heavyweight hosted CI; favor reproducible local/manual gates unless policy changes explicitly.

## 4. Candidate architecture: five layers, with responsibilities

~~~text
Human / External General Agent (ChatGPT, Codex, other harnesses)
                         |
         Human/Agent Capability Interface
       discover / inspect / plan / validate /
        prepare / launch / monitor / inspect
                         |
      Scientific Capability & Validation Layer
      scientific semantics | I/O contracts | QC
      alternatives | eligibility | compatibility
                         |
          Workflow Composition Layer
  complete pipeline | subworkflow | process module
    native adoption | bounded adapter | custom DSL2
                         |
            Execution/Operations Layer
   Nextflow | Apptainer | profiles | HPC executor
        retries/resume | resource policies
                         |
         Evidence / Provenance / Artifacts
    plan revisions | attempt states | versions |
    results | logs | QC | scientific decisions
~~~

The five layers are **responsibility views**, not necessarily five services, deployables, databases, or packages. The initial implementation should favor a lightweight Python-facing interface and Nextflow DSL2, with process isolation determined by real requirements.

The evidence system is cross-cutting: every run transition must be explainable from durable records, without relying on a still-running LLM conversation.

## 5. Scientific Capability Model — working direction

### 5.1 Capability is not identical to a Nextflow module

A **Scientific Capability** describes a scientifically meaningful task and its admissible inputs, outputs, assumptions, limitations, and verification status. A specific Nextflow implementation is one possible **realization** of that capability.

Examples:
- expression quantification (realized by an adopted upstream pipeline or specific subworkflow);
- differential expression analysis (requires valid study design and contrasts);
- peak calling and chromatin accessibility analysis (assay-specific assumptions);
- lncRNA identification (species/annotation-dependent filters and models).

An implementation may be a complete pipeline, subworkflow, process module, or a documented bundle of components. Multiple alternative implementations can serve one capability **without assuming scientific equivalence**.

### 5.2 Minimum proposed descriptor questions (not a schema freeze)

- What scientific question, assay, and result does this capability support?
- Which upstream implementation and exact revision realizes it?
- What are its native input/output contracts and required metadata?
- Which organism, reference assembly/annotation, library and sample-design constraints apply?
- Which choices are scientifically meaningful parameters versus operational settings?
- What artifacts are created and what constitutes success, warning, failure, or insufficient evidence?
- Which unit, fixture, real-container, and real biological-data validations exist?
- Which known assumptions, exclusions, compatibility constraints, and alternatives exist?
- What resource/security requirements and provenance are needed?

Prefer existing nf-core component metadata, pipeline schemas, and native Nextflow semantics as primary technical facts, with a **minimal BioWorkflows scientific annotation** where indispensable.

### 5.3 Compatibility is broader than path/format

Proposed compatibility dimensions:
- Technical: file format, index, record structure, compression, channel/cardinality/collection semantics.
- Biological: species, genome assembly, annotation version, coordinate convention, assay and sample context.
- Experimental: grouping, replicates, pairing, treatments, controls, batches, strandedness and design metadata.
- Statistical: contrast, normalization, assumptions, uncertainty, multiple-testing/independence concerns.
- Operational: container, resources, executor, storage/access constraints.
- Evidence: validation grade, known pitfalls, and whether a result is appropriate for a scientific claim.

Proposed outcomes are **compatible / needs-explicit-adapter / conditional-review / incompatible**; do not quietly convert, coerce, or ignore mismatches. A syntactically valid DAG is not proof of scientific validity.

Unresolved: how much semantic typing should live in a core generic contract versus domain-specific capability packs.

## 6. Workflow Composition Model — competing modes retained

We favor **hybrid, graded composition**, not one universal pattern:

1. **Native pipeline adoption**: preferred when a mature full pipeline directly answers the scientific question.
2. **Capability/subworkflow composition**: preferred for reusable cross-analysis chains; reuse intermediate results and explicit handoffs.
3. **Module-level composition**: available for nonstandard analysis, but with stronger type/semantic validation and resource/security constraints.
4. **New module/pipeline authoring**: exceptional escape hatch when existing components do not satisfy a justified scientific need; requires review and validation before trusted publication.

A *framework* should be a reusable analysis pattern or suggested template, **not a required immutable workflow topology**. Agents may propose new patterns, variants, and method branches, but their executable realization must pass the same controls as human-designed plans.

Composition mechanisms to investigate:
- External staged orchestration of independent upstream pipeline runs with validated published-output handoffs.
- Native Nextflow DSL2 subworkflow composition.
- New Nextflow 26.10 typed workflow modules and pipeline composition.

**Important compatibility caution:** Nextflow 26.10 stable typing and workflow-module support do not mean that existing nf-core pipelines are already type-enabled or directly composable. Official documentation currently marks pipeline composition experimental. Retain Nextflow 26.04.6 as the known repository baseline until isolated compatibility experiments justify an upgrade. The Nextflow agent primitive is also preview and is **not** proposed as a required external agent-control mechanism.

Decision needed: preferred default composition unit, when to cross whole-pipeline execution boundaries, and how to turn generated DSL2 into a reviewed, stable bundle rather than a transient code blob.

## 7. Human/Agent control model

The same public scientific operations should serve both CLI users and external agents. Illustrative operation vocabulary (not promised API):
- discover / inspect capability
- describe inputs, reference and study design
- propose analysis plan
- validate plan, alternatives and compatibility
- prepare a locked execution bundle
- authorize and launch a run
- monitor/inspect run attempts, errors and QC
- retrieve artifacts, provenance and explanations
- propose a plan revision rather than silently mutate an accepted plan

Architectural principle: **LLM for semantic interpretation and proposing choices; deterministic components for binding, validation, execution authority, provenance and acceptance**. No unique dependence on a model vendor.

Do not give an agent blanket shell, filesystem, or network authority under the guise of a capability interface. Design scoped privileges, allowlisted execution targets, preview/approval for expensive or potentially destructive actions, and policy-aware access to sensitive datasets.

Human oversight is not necessarily a manual step for every small run: approval policy should be consequence- and environment-sensitive. For example, in-memory planning and dry validations may be automatic; resource-intensive or novel workflows may require explicit authorization.

Agent-produced code, prompts, method rationale, model/provider metadata (where available), and any agent-mediated changes should be recorded when they materially affect accepted results. Prefer to keep LLMs **off the numerical scientific data path** when feasible; if an LLM is on-path, record that dependency and accept the associated reproducibility limits.

## 8. Planning, immutable execution, and evidence

Proposed lifecycle:

~~~text
Biological Question / User Intent
          |
  Candidate Plan & Alternatives
          |
 Scientific/Technical Validation
          |
  Accepted Plan Revision + Lock
          |
 Prepared Executable Nextflow Bundle
          |
     Run -> Attempt(s)
          |
  Outputs / QC / Run Provenance
          |
 Inspect / Compare / Report / Revise
~~~

Separate:
- **Analysis intent:** what the researcher asks, including ambiguities.
- **Plan revision:** chosen capabilities, method/parameter rationale, design and dependency graph.
- **Resolved implementation/lock:** exact workflow/module sources, revisions, effective params, containers, references and execution profile.
- **Run / attempt:** a concrete execution identity, attempt state, logs, hardware/environment facts, retries and failures.
- **Artifact:** materialized output plus integrity, producer, context and retention location.
- **Validation/evidence:** operational success, QC outcomes, scientific sufficiency, method limitations and explicit rejection/abstention.

Suggested record contents: original question, sample and project metadata, input/reference identities and digests, immutable parameters, Nextflow revision, Git commit(s), container digests, resource/executor settings, command, task/trace/timeline/report, outputs and their integrity, QC and human/agent decisions. Capture access and retention limitations separately.

**Do not equate reproducible execution with byte-for-byte numerical identity.** Define reproducibility tolerance appropriate to nondeterministic numerical methods and external data resources. Distinguish *run succeeded* from *scientific result qualified*.

Reuse Nextflow's built-in trace/report facilities and evaluate nf-prov exports (Workflow Run RO-Crate / BioCompute Object) before designing custom provenance serializers. BioWorkflows should principally add missing scientific decision and cross-run linkage semantics.

Failure recovery must preserve earlier attempts and accepted plan revisions; a new method, reference, parameter, or upstream pipeline revision creates a new auditable revision/attempt relationship rather than being silently overwritten.

## 9. Relationship to PepXplorer and related specialist projects

BioWorkflows and PepXplorer may share principles (reproducibility, validated capability invocation, provenance, agent-neutral interface) **without forcing either project to depend on the other's runtime**.

- BioWorkflows: broad omics analysis and reusable generic bioinformatics capabilities.
- PepXplorer: specialist noncanonical peptide analysis, evidence, and discovery infrastructure.
- Cross-project handoffs should be explicitly declared by the consuming/integrating project when scientifically needed; no assumption of automatic code-level integration.
- Future interoperability may use common versioned manifests, shared conventions or adapters, but should not be prematurely baked into BioWorkflows.

The model is independent platforms with optional evidence-bearing interoperability, **not one expanding monorepo or a BioWorkflows ownership claim over specialized PepXplorer capabilities**.

## 10. Validation and scientific trust

A trust/readiness model should distinguish at least:
1. Descriptor/documentation inspected.
2. Contract and static checks passed.
3. Module/fixture and container execution passed.
4. Integrated workflow on representative inputs passed.
5. Real biological data run with scientific QC/interpretation reviewed.
6. Cross-platform/agent-driven perturbation and replay demonstrated.

These are **discussion-level distinctions, not ratified status labels**. A capability can be usable for a bounded purpose without being validated for all organisms/assays/designs. Evidence should identify exactly what was tested.

Test both positive and negative behavior: corrupt FASTQ, wrong genome build, missing index, mismatched samples, invalid contrast, missing biological replication, contaminated directory, incompatible upstream output layout, and over-ambitious agent suggestions. Rejecting a misleading plan can be the correct result.

Run-success, scientific-adequacy, and biological-inference confidence must remain separate.

## 11. Evidence-informed external ecosystem (checked 2026-10-31)

| Reference | Relevance / limits |
|---|---|
| Nextflow 26.10 migration / new features (https://docs.seqera.io/nextflow/migrations/26-10) | Stable static typing; workflow modules and pipeline composition; experimental agents. Not automatic compatibility with existing code. |
| Typed workflows (https://docs.seqera.io/nextflow/workflow-typed) | Official pipeline-composition docs mark feature experimental; trial behind controlled qualification. |
| nf-core component guide (https://nf-co.re/docs/developing/components/creating-components) | Reuse nf-core modules/subworkflows, component metadata and testing practices rather than duplicating infrastructure. |
| nf-prov (https://github.com/nextflow-io/nf-prov) | Workflow Run RO-Crate and BioCompute Object exports; version compatibility must be tested. |
| Workflow Run RO-Crate (https://pubmed.ncbi.nlm.nih.gov/39255315/) | Interoperable workflow-run provenance; does not replace scientific interpretation. |
| BioMaster, Patterns 2026 (https://pubmed.ncbi.nlm.nih.gov/42630761/) | Planning, execution, error recovery and tool knowledge; its reported task-completion benchmark does not itself establish biological truth. |
| BioAgent Bench, ICML 2026 (https://proceedings.mlr.press/v306/fa26a.html) | End-to-end evaluation and robustness to corrupted inputs, decoys and prompt perturbation. |
| THREAD-Bio discussion, Biology 2026 (https://pmc.ncbi.nlm.nih.gov/articles/PMC13565537/) | Traceable links between scientific question, data, decisions, approvals, validation and claims; reporting framework, not a mandatory runtime schema. |
| Reproducible LLM-assisted scientific analyses (https://pubmed.ncbi.nlm.nih.gov/42746527/) | Preserve generated assets, record models/versions, and prefer off-data-path LLM assistance where possible. |

External facts are reference inputs to architecture discussion, **not immutable dependencies or guaranteed feature support**. Recheck their versions and licenses before development.

## 12. Initial architecture questions requiring subsequent deliberation

**Q1 — Capability scope and granularity:** Scientific capability as central semantic unit; default realization a pipeline/subworkflow, with module-level freedom. How to define capability families without hardcoded assay silos or overly abstract ontologies?

**Q2 — Native interfaces and common contracts:** What is the minimum cross-run/agent-ready envelope while preserving upstream samplesheets, params and output semantics? Are multiple adapter contract classes needed?

**Q3 — Composition compiler:** Should the default engine compose whole pipeline handoffs, native DSL2 subworkflows, typed Nextflow 26.10 pipelines, or choose per case? What validation must occur before emitting executable DSL2?

**Q4 — Scientific knowledge/decision authority:** Who curates method recommendations, parameter bounds, alternative implementations, and method eligibility? How are evidence freshness and contradictions handled?

**Q5 — Accepted plan vs adaptive agent:** Which changes are allowed within the same attempt (retry/restart), and which must create a new plan revision? What approval and stopping rules apply?

**Q6 — Evidence and reproducibility:** Which minimal identities, hash rules, execution snapshots, result lineage, retention settings and provenance export formats are needed in a first useful version?

**Q7 — Trust and validation:** How to grade science-specific verification versus generic execution verification; what real biological test set should qualify first capabilities?

**Q8 — User-facing integration:** Does a small Python CLI with optional MCP/skills suffice? When, if ever, is a service/database needed? Avoid service-first design.

**Q9 — Security/HPC operations:** User-controlled profiles, Apptainer/PBS, quotas, network restrictions, data privacy and human authorization. How to ensure no unintended queue submissions while preflighting?

**Q10 — Cross-project interoperability:** Which conventions may be shared with PepXplorer without introducing build/runtime coupling?

Questions Q1–Q6 are foundational; other decisions can be deferred if they do not affect initial compatibility contracts.

## 13. Proposed discussion-to-implementation sequence (not a delivery commitment)

| Phase | Nature | Exit evidence |
|---|---|---|
| D0 — Unified exploration baseline | **This document** | Scope, existing assets, assumptions, contradictions, open questions documented |
| D1 — Repository and ecosystem reality review | Discussion/audit | Full current code/validation inventory; dependency map; exact migration baseline conflict handling |
| D2 — Scientific Capability & Compatibility Model | Architectural decision | Examples across RNA/ATAC/ChIP/lncRNA; ability to represent valid and invalid handoffs |
| D3 — Composition & Nextflow feasibility | Architectural experiments | Demonstrated alternatives on controlled fixtures; explicit 26.04.6 vs 26.10 compatibility evidence |
| D4 — Plan, Execution, Provenance & Agent Interface | Architectural decision + thin prototype | Plan/revision/run/evidence state transitions, lifecycle permissions and minimal interfaces |
| D5 — Integrated architecture closeout | **Future freeze decision** | Cohesive responsibility map, selected contracts, unresolved risks, incremental migration strategy |
| I0+ — Staged implementation and qualification | Implementation **after** D5 review | Reuse existing pipelines, local/manual tests, then representative real biological data |

Proposed first end-to-end scientific demonstrator: **RNA-seq -> differential analysis -> optional enrichment / lncRNA handoff**, with explicit sample/reference contracts, locked versions, captured QC and provenance. Add **ATAC-seq -> differential accessibility / regulatory association** only once corresponding capabilities and science contracts are qualified.

Why this order: existing adopted and custom routes already exercise full-pipeline adoption, bounded extension and multi-stage handoff. It offers a realistic test of the architectural thesis without a speculative framework rewrite.

## 14. Working principles / current level of agreement

**Directionally supported (not formal freeze):**
- Transform BioWorkflows from personal workflow collection toward a reusable, general-purpose agent-ready scientific infrastructure.
- Nextflow DSL2 + nf-core/community ecosystem are the execution/reuse foundation.
- Allow flexible workflow topology assembled from different levels of reusable scientific capability.
- Scientific correctness, compatibility, reproducibility, version/parameter control, and inspectable deterministic run records are required.
- Humans and general agents access the same validated capabilities.
- Do not prematurely rewrite existing workflows; leverage them for architecture qualification.
- Continue substantive discussion to an integrated architecture before implementation.

**Candidate preferences, to be evaluated rather than silently ratified:**
- Pipeline/subworkflow-capability composition by default; module-level escape hatch.
- Deterministic-first contract binding and validation; LLM for bounded semantic assistance.
- Locally deployed, Python-facing control interface with optional skills/MCP adapters.
- Apptainer/HPC and manual local validation first.
- Provenance via native Nextflow plus standards-aware extensions.

**Not agreed/frozen:**
- Exact object schemas, naming and IDs.
- Exact package layout, registry storage, server/database architecture.
- The final Nextflow target version and pipeline-composition mechanism.
- The scientific grading ladder/agent authorization model.
- Formal replacement or amendment of the October migration baseline.
- Detailed implementation batches or code changes.

## 15. Immediate next discussion agenda

1. Clarify the **Scientific Capability Model** with 3–4 concrete examples and counterexamples: RNA-seq quantification, differential expression with/without replicates, ATAC peak calling, and lncRNA discovery.
2. Define the **Workflow Composition Model** with at least two competing design alternatives and their feasibility/maintenance costs.
3. Resolve what belongs to deterministic contract validation and what can be advised by an agent.
4. Only then refine plan/run/provenance models and define the minimal prototype.
5. Once major decisions converge, produce a coherent project architecture, independent critical review and explicit specification freeze; implementation follows, using the current repository as migration source.

**No code changes, runtime migration, or auto-CI addition are authorized by this baseline.**
