# BioWorkflows — D2 Unified Scientific Capability & Workflow Composition Discussion

- Date: 2026-10-31
- Status: **WORKING / DISCUSSION / NOT FROZEN**
- Repository: Chengy257/BioWorkflows
- Branch: docs/agent-ready-infrastructure-discussion-20261031
- Source snapshot: main at 8839e877a451b51eca460f937a07b34dbe65b107
- Parent: [Agent-Ready Infrastructure Unified Discussion Baseline](BIOWORKFLOWS_AGENT_READY_INFRASTRUCTURE_ARCHITECTURE_DISCUSSION_BASELINE_2026-10-31.md)
- Supporting notes: [D2-A capability examples and counterexamples](BIOWORKFLOWS_D2A_SCIENTIFIC_CAPABILITY_MODEL_WORKING_DISCUSSION_2026-10-31.md) — retained for history, **not** a separate approval gate.
- Intent: One **unified D2 discussion document** for capability semantics and workflow composition. Future findings should amend this document, rather than create D2-B/D2-C review fragments.
- This is not a specification freeze, runtime acceptance, API contract, authorization to refactor code, or permission to add automated CI.

## 1. Core position and question

BioWorkflows should let humans and general-purpose AI agents discover and compose scientifically justified analyses using the existing Nextflow/nf-core ecosystem. It should not prescribe one rigid pipeline, create a competing execution engine, or allow an LLM-generated command alone to define scientific validity.

**Working hypothesis:** An agent builds a *scientific capability plan*, while a deterministic and reviewable binding step maps it to suitable pinned Nextflow pipelines, subworkflows, processes, adapters and data. Actual execution and scientific adequacy remain separately evidenced.

We should distinguish flexibility of **scientific analysis topology** from mutability of an **accepted version-locked run plan**.

## 2. Grounding in the existing code

Source/document review against main@8839e877 (not new runtime qualification):

| Current code/route | Compositional implication |
|---|---|
| Adopted nf-core/rnaseq plus rnaseq-ext | Native upstream pipeline with local extensions; authoritative Salmon vs featureCounts output selection and particular published BAM/matrix paths |
| prepare_da_inputs.py and adopted differentialabundance | Explicit count-matrix/sample/contrast adapter; contrast orientation, replication and model input requirements matter; full DA validation gap remains |
| Adopted chipseq plus chipseq-ext | Revision-sensitive output directory/BAM/peak paths; replicate and control rules; recorded real-container fixture execution |
| Adopted atacseq plus atacseq-ext | Local organelle QC, spike-in and footprint routes; several real-container/full-run validations remain deferred |
| lncrna-nf local DSL2 pipeline | Uses published RNA BAMs, reference FASTA/FAI and GTF, protein-coding annotation, strandedness, optional species-specific CPAT model and auxiliary databases |
| nf-core module meta.yml versus local meta.yml | Technical metadata shape and completeness vary; must inspect rather than automatically treat every metadata record as a qualified scientific capability |

Existing Nextflow DSL2 code, scientific decision notes and validation receipts are assets. Do not rewrite validated routes to satisfy a speculative universal framework. Older fixture/smoke results do not establish publication-ready scientific validity for all contexts.

## 3. One unified conceptual model

| Concept | Responsibility |
|---|---|
| **Scientific Capability** | A bounded scientific analytical function: purpose, meaningful outputs, scope, assumptions, exclusions, scientific QC |
| **Implementation Realization** | Exact versioned pipeline/subworkflow/module/adapter or composition that realizes the capability for stated circumstances |
| **Semantic Artifact Role** | Scientific meaning of the data passed between analyses: measurement unit, reference, biological samples and producer evidence |
| **Binding** | Selected realization(s), exact versions, source-output to consumer-input mappings, adapters, data/reference and parameter choices |
| **Analysis Plan** | Request-specific scientific capability graph, alternatives, rationale and expected results; accepted revision becomes immutable |
| **Execution and Validation Evidence** | Concrete runs/attempts, effective parameters, outputs, provenance, QC and qualified scientific interpretation |

A scientific capability does **not** equal one Nextflow process or one named program. One full pipeline may satisfy several scientific functions; a capability may have multiple method alternatives. Multiple realizations must **not** be presumed scientifically equivalent.

Native nf-core interfaces (schemas, samplesheets, process IO, output publication, container conventions) remain source-specific authority. BioWorkflows adds only minimal semantic annotations, compatibility checks and explicit mappings that solve observed real needs.

No database, universal samplesheet, giant ontology, fixed schema names, CLI signatures or package hierarchy is frozen here.

## 4. Scientific capability and implementation contracts

Proposed small common descriptor vocabulary:

- Scientific intention, admissible analytical questions, exclusions and output interpretation.
- Technical IO and *semantic roles* (not just file suffix); required indexes and supplementary artifacts.
- Biological context: species, genome assembly, annotation/feature-ID namespace, assay, library and strandedness when relevant.
- Experimental/statistical context: biological replicate relationships, sample pairing, group/batch, controls, contrasts and method assumptions.
- Realization source: upstream project, exact commit/tag, Nextflow entry, parameters and dependencies; version-specific published output contracts.
- Trust and evidence: exact scope of static/fixture/container/real-data qualification and known limitations.
- Scientific method and operational parameters distinguished for change control.

A realization may expose only a subset of possible scientific use cases. Versioned upstream migration can invalidate adapters or evidence. A file parser or generic Path type does not establish scientifically compatible results.

**Candidate semantic artifact roles** include AlignedReadSet, GeneCountMatrix, ExpressionEstimate, PeakSet, DifferentialResult, ReferenceBundle and QCReport. These are explanatory working examples; they do not replace actual BAM/CRAM, GTF/GFF, matrices or upstream native artifacts.

## 5. Compatibility is more than typed IO

Checks should consider technical format/index/cardinality; biological reference and assay; experimental sample/control/replicate design; statistical model and units; runtime/security eligibility; and qualification evidence for the intended use.

Working decision vocabulary: **ACCEPT**, **ADAPT** (only a documented valid operation), **REVIEW** (missing evidence or legitimate judgement), **REJECT** (known contradiction). These are candidate meanings, not a frozen enumeration. Missing information differs from known scientific incompatibility. A human approval cannot turn an invalid inference into a sound scientific claim.

Important negative cases:

1. TPM-only matrix flowing into DESeq2 raw-count analysis: a TSV may parse while the inference is invalid.
2. No replicated biological groups or fully confounded batch: a naive differential inference claim is not justified even when sample metadata parse.
3. RNA BAM aligned to assembly A with incompatible lncRNA annotation for assembly B: require a qualified reanalysis or validated conversion, not simple relabeling.
4. Human-trained CPAT model applied without justification to rice: incompatible species-conditioned realization.
5. ATAC peaks represented as ChIP peaks because both are BED-like intervals: wrong scientific semantics.
6. Differential expression plus nearby ATAC peaks: a useful association hypothesis, **not** automatic proof of direct regulation.
7. Missing published upstream BAM: the downstream handoff is not executable; replan upstream instead of fabricating a file.

For actual scientific alternatives, an agent can propose methods, identify missing information or suggest a revised plan. It cannot bypass deterministic scientific constraints or unilaterally change an accepted contrast/reference.

## 6. Three distinct but related graphs

**Scientific intent graph:** a DAG of meaningful capabilities and the biological questions they answer.

**Binding graph:** chosen realization(s), exact native IO, compatible semantic artifacts, adapters, parameter and reference locks, rationale and unresolved conditions.

**Execution graph:** one or more concrete Nextflow runs, their subworkflows/tasks, container/environment resources, retries, outputs and provenance.

There is no obligation to map one scientific capability to one Nextflow process, or one entire analysis plan to one Nextflow run. Scientific flexibility is preserved *above* the concrete execution topology.

## 7. Composition strategies, compared within one D2 topic

| Strategy | Benefits | Difficulties | Working use |
|---|---|---|---|
| **Native standalone upstream pipelines with artifact handoffs** | Maximum adoption of maintained nf-core routes; upgrade isolation | Published-output version drift; materialized data; cross-run lineage/coordination | Sensible default candidate when mature upstream already solves the main assay |
| **One DSL2 graph composed of subworkflows/processes** | Efficient dataflow and flexible custom analysis | More interface coupling; validation burdens and potential duplication of upstream implementation | Suitable for bounded new analyses, custom scientific extensions and already reusable local modules |
| **Typed workflow modules under Nextflow 26.10** | Structured module reuse and type checks | Migration and qualification of existing routes required | Optional future compatible path |
| **Pipeline inclusion/composition in Nextflow 26.10** | Could combine full pipelines under one Nextflow graph | Feature is experimental; typed enablement, assets/config and native upstream compatibility require tests | Isolated feasibility investigation, **not a prerequisite** |
| **Arbitrary AI-generated code each run** | Maximum superficial freedom | Scientific validity, security, provenance and reproducibility risks | Controlled development escape hatch, not a trusted default |

**Provisional candidate: Hybrid composition chosen at binding time**. Prefer native nf-core pipeline execution with checked handoffs when it suffices; use DSL2 subworkflows/modules where they provide genuine flexibility and value. Do not require one giant multi-assay Nextflow pipeline.

Official 2026-10-31 documentation checked:
- [Nextflow 26.10 migration](https://docs.seqera.io/nextflow/migrations/26-10): non-preview static typing, workflow modules, pipeline composition, preview agent primitive.
- [Nextflow typed workflows](https://docs.seqera.io/nextflow/workflow-typed): pipeline composition explicitly **experimental**; typed caller and included pipeline requirements.
- [nf-core component specifications](https://nf-co.re/docs/specifications/components/overview): technical module/subworkflow metadata and testing conventions.

BioWorkflows documented Nextflow baseline remains **26.04.6**. Do not treat the new typing/pipeline inclusion features as validated compatibility with existing adopted nf-core revisions. Nextflow lint and separate scientific validation are necessary; no upgrade without a standalone compatibility/real-container check. The Nextflow agent primitive is not a dependency for a general external agent interface.

## 8. Explicit adapter and data-handoff model

Every cross-boundary connection should declare:

- exact producer and consumer realization identities and versions;
- produced artifact identity/semantic role, actual path and required sidecars/indexes;
- reference/feature/sample/experimental context and cardinality;
- consumer expectations and whether actual data satisfies them;
- explicit adapter operation, exact version, scientific consequences and validation;
- verification timing (pre-run, after upstream artifacts, before downstream admission);
- outcomes, including rejection and requested re-planning.

A lossless format/schema adapter might remain an implementation detail; a transformation altering scientific inference (identifier ambiguity, aggregation, normalization, association modeling) should be represented as a qualified scientific operation with inspectable parameters.

### Case A: RNA-seq -> differential expression -> enrichment

Use existing adopted rnaseq and local quantitative route switch, adapted differentialabundance inputs, and optionally enrichment. Check raw count meaning, group/contrast direction, sample matching, design estimability and gene identifiers before executing dependent analysis. The existing differentialabundance validation gap remains a gate before any broad real-data readiness claim.

### Case B: RNA-seq -> lncRNA identification

Adopted rnaseq publishes BAMs and index when configured appropriately; local lncrna-nf uses them with species/reference/model conditions. Preserve native upstream pipeline and local composite DSL2 workflow; qualify exact handoff.

### Case C: RNA-seq + ATAC-seq -> regulatory association

Potential future composition: parallel mature upstream pipelines, then a separately described integration capability. Confirm assembly, annotation, sample/condition and contrast relationships. Peak-to-gene proximity must be labeled as association, not direct causation.

### Case D: New method absent from nf-core

Allow an agent to propose a module/subworkflow in a bounded development workspace. It is **unqualified code** until review, source/interface/static tests, container execution and appropriate scientific tests. Direct trust promotion or shared-HPC unrestricted execution is not warranted.

## 9. Plan lock, execution and evidence consequences

The scientific composition should end in a reviewable **accepted Plan Revision** with a resolved source and parameter lock before authorized execution. Subsequent plan-method/reference changes must be new revisions rather than invisible retries.

A plan may launch several separate Nextflow runs. Preserve cross-run producer/consumer lineage, task outputs, QC, signed or hashed input/reference identities, exact source revisions and containers. Use native Nextflow trace/report capabilities; evaluate [nf-prov](https://github.com/nextflow-io/nf-prov) for per-run RO-Crate/BioCompute provenance, adding only the cross-run scientific decision layer missing from native provenance.

A process exit code, software compatibility, QC pass and justified biological conclusion are separate states. Data-level numerical reproducibility may have scientific tolerances even where the workflow and input manifest are exact.

Do not design a new scheduler; keep Nextflow as the execution engine and a thin, probably Python-facing local interface as a candidate. Do not introduce automatic GitHub Actions. Avoid accidental HPC submission during validation; earlier repository notes warn that Nextflow preview under machine-local PBS settings can submit real jobs. Keep local/manual gates.

## 10. Design risks and future qualification evidence

- **Over-abstract capability model:** test against actual repo cases, including invalid statistical designs and incompatible reference versions.
- **Duplicated upstream schemas:** use nf-core/native metadata as source of technical truth, add only necessary semantics.
- **Experimental runtime dependency:** keep independent-pipeline handoff viable as a durable alternative.
- **Silent scientific drift:** accepted plan locks, visible adapters and change receipts.
- **Unbounded code generation:** separate code-proposal authority from execution and scientific qualification.
- **Too much framework for a personal repo:** small Git-versioned capability catalog/local index first; no mandatory centralized service/database until real need.
- **Weak end-to-end validation:** actual biological data plus negative admission tests are essential, not merely synthetic fixture PASS.

Useful future technical proof—not being executed in this discussion—should compare a pinned native RNA->DA handoff against a small composed local subworkflow, test version-drift failure, and demonstrate replay/inspection of a cross-run plan. These experiments belong to implementation feasibility **after** architecture priorities are agreed.

## 11. Only three owner-level choices worth raising

**Decision 1: Default composition boundary (primary).**
- **A (provisionally favored):** use mature nf-core pipelines as independent, native execution units; compose subworkflows/modules when justified; permit hybrid plans.
- **B:** prioritize one integrated DSL2 DAG with deeper extraction of upstream internals.
- **C:** no default, let each agent dynamically choose without a baseline preference.

This choice meaningfully affects maintenance cost, reuse fidelity and provenance architecture. The discussion leans A with hybrid flexibility but **does not freeze owner approval**.

**Decision 2: Agent-created methods.** Permit unqualified code proposals, with controlled environment, approval and validation before any trusted use; avoid assuming every agent is allowed to run novel code on the HPC. Authorization details can be decided with the later Agent Interface discussion.

**Decision 3: Catalog authority.** Prefer a small Git-versioned capability manifest/index referencing native sources, not a central database or duplicate ontology from the start. Reassess only if real scale or multi-user access demands it.

Minor field naming, directory layout, ID generation, service packaging and exact schema shape should **not** trigger separate review conversations.

## 12. D2 record status and transition

**Completed for discussion:** unified capability concept, alternative composition modes, version-sensitive boundaries, scientific negative cases, explicit adapter principles, preliminary default preference, and consequential decisions.

**Not completed:** formal selection/freeze, a comprehensive scientific rulebook, real-data qualification, compatibility performance runs, implementation contract or code.

**Document discipline:** maintain this document as the single editable **D2 working record**. Keep D2-A as a historical detailed input; no more D2 micro-stages or micro-closeouts. Revisit D2 only if subsequent full-project discussions expose a real contradiction.

The next broad architecture topic should be **Plan / Run / Evidence / Agent Control** as one integrated discussion, followed by one project-wide consolidated critical review/freeze before incremental implementation.
