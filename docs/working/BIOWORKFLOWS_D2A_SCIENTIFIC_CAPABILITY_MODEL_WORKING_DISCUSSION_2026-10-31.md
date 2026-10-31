# BioWorkflows D2-A — Scientific Capability Model: Meaning, Realizations, and Contracts

- Date: 2026-10-31
- Status: **WORKING DISCUSSION — NOT FROZEN**
- Parent: [Unified Agent-Ready Architecture Discussion Baseline](BIOWORKFLOWS_AGENT_READY_INFRASTRUCTURE_ARCHITECTURE_DISCUSSION_BASELINE_2026-10-31.md)
- Repository reference: main at 8839e877a451b51eca460f937a07b34dbe65b107
- Scope: Scientific capability semantics, realization boundaries, candidate compatibility model; not an implementation spec.
- Governance: Alternatives and open questions remain live. No owner ratification or code changes are implied.

## 1. Central problem and working hypothesis

A general agent must choose analytical functions based on a scientific question, not simply match tool names, file suffixes, or Nextflow process identifiers. A workflow module is a technical implementation unit; a scientific capability is a description of **what a valid analysis can establish, for which data and scientific context, under which assumptions and with which evidence**.

**Working conceptual split:**

1. **Scientific Capability** — analytical objective and bounded scientific meaning (for example, gene-level RNA expression quantification).
2. **Implementation Realization** — a concrete versioned nf-core pipeline, subworkflow, process/module, bounded adapter, or local workflow that can realize the capability within a stated applicability domain.
3. **Analysis Plan** — a request-specific selection of capabilities, realizations, samples, references, scientific method, and parameters.
4. **Execution & Validation Evidence** — observed attempts, results, failures, QC, scientific limitations, and version-specific qualification.

These are conceptual responsibilities, **not a proposed four-table schema**. One capability can have many realizations; a full pipeline can realize several capabilities; a realization can cover only a subset of possible use cases.

## 2. Evidence grounded in this repository

Inspection of current source, not a new runtime acceptance:

| Existing surface | Actual semantics relevant to a new model |
|---|---|
| nextflow/pipelines/rnaseq-ext/README.md and workflows/ext.nf | Explicit Salmon vs featureCounts authoritative matrix route; local optional stages depend on adopted nf-core/rnaseq published outputs |
| nextflow/pipelines/rnaseq-ext/bin/prepare_da_inputs.py | Adapts count matrix + sample/group/batch metadata into differentialabundance inputs; contrast orientation matters; not merely file renaming |
| nextflow/pipelines/rnaseq-ext/modules/nf-core/subread/featurecounts/meta.yml | nf-core technical inputs/outputs and format information, but not sufficient constraints for every downstream statistical interpretation |
| nextflow/pipelines/chipseq-ext/README.md | ChIP/CUT&Tag extension stages depend on precise upstream BAM/peak publish paths and control/replicate semantics |
| nextflow/pipelines/atacseq-ext/README.md | Organellar read fractions, optional spike-in/footprint stages and conditional multi-group requirements |
| nextflow/pipelines/lncrna-nf/README.md and workflows/identify.nf | BAM, genome+FAI, GTF, coding-gene annotation, strandedness, optionally species-specific CPAT model and auxiliary databases |
| nextflow/pipelines/atacseq-ext/modules/local/qc_gates/meta.yml | Local metadata currently less structured than richer nf-core module meta descriptions: extraction must be audited rather than assumed uniform |

Current scientific decision records and validation notes already distinguish tool capabilities and tests, but there is no repo-wide validated capability registry or consistent semantic handoff contract. A documented smoke/fixture pass is not a blanket real-data validation.

## 3. Three candidate models; keep alternatives visible

| Candidate | Strength | Weakness | Working disposition |
|---|---|---|---|
| A — software/module catalog | Easy to derive from existing nf-core metadata and process code | Tool-centric, cannot judge biological or statistical eligibility | Keep as **technical inventory**, not user-facing scientific authority |
| B — typed technical DAG | Better structural checkability and pipeline assembly | A TSV or Path type does not encode assay, unit, design, reference, or claim | Use as **technical composition mechanism**, never full scientific validation |
| C — semantic scientific capability with versioned realizations | Aligns agent intent with science; supports alternative methods and composition levels | Risks excess ontology/curation | **Provisional preference**, test against cases before freezing |

Avoid creating a giant formal ontology or proprietary tool registry. A small stable capability descriptor with scoped domain annotations is a more plausible first milestone.

## 4. Tentative contract responsibilities — NO schema freeze

**Core capability concept (candidate fields):**
- stable human-meaningful purpose, scientific output and explicit exclusions;
- semantic input/output roles and required supporting evidence;
- admissible biological/experimental contexts and method constraints;
- realization choices and the limits under which each applies;
- expected QC and scientific sufficiency statements;
- provenance reference and scoped verification status.

**Realization binding (candidate fields):**
- source location, exact revision, Nextflow entry point/granularity, native parameter/interface requirements;
- containers, software revisions, executor/resource constraints;
- output publication contract and required adapters;
- scientific applicability subset, known limitations and validation receipts.

**Plan and run are distinct**: the plan binds context and scientific decisions; the run captures what actually executed and its evidence. Stable capability identity must not silently imply one canonical implementation or one parameter choice.

**Native authority remains native**:
- adopted nf-core samplesheets, pipeline parameter JSON schemas, module meta.yml, published-output behavior and version pins remain the source of truth for those specific implementations;
- BioWorkflows should annotate and validate the scientific *binding* of such interfaces, not fork a universal replacement samplesheet;
- adapters are explicit versioned scientific/technical operations with test evidence, not implicit parsing conveniences.

The pipeline compatibility boundary is keyed to the exact upstream revision. A revised upstream published output may invalidate an existing adapter until independently requalified.

## 5. Semantic artifacts as the narrow cross-workflow bridge

A semantic artifact role describes what a produced file or collection **means**, without replacing its actual native file format.

| Candidate role | Physical form | Important semantics |
|---|---|---|
| AlignedReadSet | BAM/CRAM plus index | sample/read processing, coordinate reference, alignment policies, pairedness, strand |
| GeneCountMatrix | TSV/HDF5 or equivalent | integer read counts vs estimated counts, feature identity, sample identity, quantification route |
| ExpressionEstimate | matrix | normalized expression units (TPM/FPKM etc.), not interchangeable with raw counts |
| PeakSet | BED/narrowPeak/broadPeak | assay, reference assembly, calling conditions, score definitions |
| DifferentialResult | table | reference/target contrast orientation, effect definition, statistical model, multiple testing |
| GenomeReference / Annotation | FASTA/GTF/GFF plus index | exact assembly/build and annotation lineage, contig/gene-ID compatibility |
| QCReport | JSON/TSV/HTML | metric applicability, thresholds, observed and inferred limits |

It is a **candidate vocabulary**, not a canonical ontology. Machine-readable metadata should capture only the distinctions required for safe selection and composition. Reuse EDAM or another maintained ontology when useful for discovery; an ontology identifier alone is not a complete scientific compatibility rule.

## 6. Stress-test examples and negative cases

### R1. RNA-seq expression quantification

nf-core/rnaseq native quantification and the local featureCounts stage can realize related gene-quantification aims, but produce measurements with method-specific semantics. An agent must not switch the authoritative matrix simply because both outputs have a TSV extension. Track quantification definition, annotation, feature/gene namespace, strandedness and relevant library conditions. Multiple realizations are alternatives, not universal numerical equivalents.

### R2. Differential expression

The downstream DESeq2 route is conditional on suitable counts, experiment metadata, contrasts, and an identifiable design. TPM-only matrices must not be silently fed into a raw-count DESeq2 inference; no-replicate group comparisons cannot be advertised as standard replicate-based differential inference; a batch completely confounded with condition is not cured by having valid-looking CSV headers. Inverting reference and target inverts scientific interpretation. When a study cannot support a requested inference, offer a bounded descriptive alternative instead of false success.

### R3. ATAC-seq and ChIP-seq peaks

An ATAC PeakSet and a ChIP PeakSet may both be valid genomic intervals, but they measure different phenomena. Cross-analysis integration requires the assay, biological question, reference and peak-calling conditions. Plant organelle read fractions can be valuable ATAC QC, but thresholds and interpretation should remain assay/context-specific. TOBIAS motif/group prerequisites should be checked before composition.

### R4. lncRNA identification

The local lncRNA pipeline depends on consistent BAM/reference/annotation plus coding-potential and homology settings. A human CPAT model must not be auto-bound to rice. Disabling Pfam or DIAMOND filtering changes the evidential interpretation; "candidate noncoding transcript" must not be promoted to "experimentally validated lncRNA" merely because processes exit 0.

### X1. RNA-seq + ATAC integration

Linking gene-expression changes and accessible regions requires reference, gene-ID, study-design and contrast compatibility. A peak-to-gene proximity association is not proof of causal regulation. Cross-assay joining may need a distinct scientific integration capability rather than a generic BED-to-TSV join.

## 7. Working compatibility semantics

Check compatibility at separable layers:
1. **Technical**: files, indexes, schemas, channels, cardinality, expected outputs.
2. **Biological**: species, genome/annotation, sample/assay, library and feature definitions.
3. **Experimental**: pairedness, replicates, grouping, controls, treatment, batch.
4. **Statistical/method**: units, assumptions, contrast, model applicability.
5. **Operational**: execution, environment, resource/security constraints.
6. **Evidence**: validated scope appropriate to the intended use (exploratory vs publication).

Candidate validation outcomes (not frozen):
- **ACCEPT**: all necessary verifiable conditions met for the declared purpose.
- **ADAPT**: explicit scientifically justified, tested adapter is needed.
- **REVIEW**: required evidence is incomplete or requires documented domain judgment.
- **REJECT**: demonstrably incompatible for the proposed route or claim.

**Unknown must not silently become acceptable**, but an unknown is not always an irreparable scientific contradiction. A review state is not itself execution authorization; an accepted step is not a guarantee of strong biological inference.

Example counterexamples:
- TPM-only matrix -> DESeq2 raw-count input: technical table may parse; scientific use **REJECT**.
- BAM from assembly A -> annotation from different build B: **REJECT** until an explicit correct reanalysis/qualified transformation is established.
- Missing published BAM for rnaseq -> lncRNA: **REPLAN** upstream/alternative realization, not fabricate a path.
- Cross-omics interval proximity -> categorical direct regulation claim: association could be computed, but the unsupported **causal claim must be rejected or qualified**.

Parameter changes should also be classified: scientific/method, implementation, input/reference, and operational. All effective values are recorded; changing method/reference cannot be disguised as a simple resume/retry.

## 8. Role of deterministic validation versus model-assisted judgment

Prefer code/rules for file integrity, sample identity, reference binding, design-rank checks, units when declared and verified, published-output contracts, exact versions and authorization. Agent assistance can interpret intent, propose candidate capabilities, weigh method alternatives, explain missing data and draft revisions. The agent's assertion of validity is never itself binding scientific evidence.

Do not demand an all-encompassing static biological checker; capability-specific decision logic and recorded expert review will still be necessary.

## 9. Design risks

- Semantic layer grows into an unmaintainable ontology.
- Capability records duplicate upstream schemas and drift after nf-core upgrades.
- False confidence from static type checking or fixture-only validation.
- A capability catalog implies arbitrary tool combinations are safe.
- Over-generalized "gene count"/"BAM"/"peak" roles erase scientific details.
- Artificial one-to-one capability-to-process mapping prevents scientific alternatives.
- Agent autonomy silently changes statistical design, reference or accepted parameters.

Mitigation direction: begin with 3–4 supported capabilities in this repository, real negative compatibility cases, revision-specific binding and explicit abstention.

## 10. Working decisions and unresolved questions

**Provisional preferences, not ratified**: semantic capability as interface; implementations separately versioned; native nf-core contracts retained; minimal semantic artifact roles; two-layer generic + domain checks; deterministic execution admission; LLM assistance only within bounded tasks.

**Open questions for D2-B/C and integrated review:**
- Does one realization expose multiple capabilities, and how should that be advertised without duplicates?
- What is the smallest cross-domain required metadata envelope?
- Which semantic rules are core checks versus assay-specific packs?
- How to express conditionally valid scientific branches and method alternatives?
- Should registry entries be Git-versioned text initially, with no database?
- When should a scientific operation be a capability versus just a technical adapter?
- Does the same run always need a fully resolved realization at planning time, or can alternatives remain until preparation?
- What is the required evidence grade for an agent to propose/prepare/execute a capability?

## 11. Next discussion sequence

- **D2-B**: composition mechanisms, binding boundaries, adapters, and Nextflow feasibility; compare alternative architectures before selecting.
- **D2-C**: adversarial capability/compatibility review against real assay cases, refinements and unresolved science.
- **D3**: execution/plan/provenance and deterministic version-lock semantics (subject to later schedule).
- **Final integrated review**: decide and freeze only after alternatives, constraints and evidence have been consolidated.

No execution, schema migration or scientific qualification is claimed by this document.
