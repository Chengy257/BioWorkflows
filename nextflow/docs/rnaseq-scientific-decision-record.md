# nf-core/rnaseq scientific rebaseline record

- **Date**: 2026-10-11
- **Status**: RATIFIED. The owner resolved all decision points D1-D8 in two
  review rounds on 2026-10-11 (decision log, Section 8). Implementation
  specifications are written against this frozen version.
- **Mandate**: the per-workflow scientific rebaseline record required by the
  migration design doc (Section 9), for the Stage-2 RNA-seq route. Drafted in
  `plans/2026-10-11-stage2-rnaseq-migration-plan.md` and frozen here per the
  Stage-1 order (record precedes code).
- **Pipeline under review**: nf-core/rnaseq 3.27.0 (released 2026-09-23;
  latest stable). The tag still must pass the Nextflow 26.04.6 parse gate at
  implementation (WP0); the fallback is the dev branch per the Stage-1 policy
  (recorded precedent: chipseq 2.1.0 does not parse under Nextflow 26).

---

## 1. Biological use case (Q1)

Bulk RNA-seq of rice (*Oryza sativa*) and other species the user configures
locally: expression quantification, QC reporting, differential expression for
treatment-vs-control designs, and downstream GO/KEGG/GSEA enrichment.
Transcript assembly is NOT part of this route (it moves to the future lncRNA
workflow). Reference genomes are never fixed in the pipeline; species presets
live in machine-local/project configuration (design doc D4 and Section 15
amendment 3).

## 2. Inputs and scientifically meaningful outputs (Q2)

- **In**: raw FASTQ (PE/SE), reference FASTA + GTF + aligner index (or
  `--save_reference`), sample sheet with group/batch metadata; optional UMI
  parameters per project.
- **Out**: aligned BAMs (published for the lncRNA handoff), gene-level
  quantification matrices (authoritative route selected per project by the
  D1 switch: native Salmon summary or local featureCounts matrix), MultiQC +
  QC plots, differential-expression tables with contrasts and thresholds,
  GO/KEGG/GSEA enrichment reports.
- Between-group DEG set comparison and transcript assembly are DELIBERATELY
  DROPPED from this route (D4; D5).

## 3. What the ecosystem provides (Q3)

Parameter pages verified 2026-10-11 (docs.nf-co.re 3.27.0 / 2.0.0); per the
Stage-1 discipline, every fact used for wiring must be re-verified against
the pinned revision's source (`workflows/*.nf`, `conf/modules.config`)
before it becomes a consumed-path contract.

**nf-core/rnaseq 3.27.0**:

- `--aligner` accepts `star_salmon` (default), `star_rsem`, `hisat2`,
  `bowtie2_salmon`; `--pseudo_aligner` adds salmon/kallisto alongside.
- **No STAR+featureCounts quantification route exists**: within STAR routes
  the native quantifiers are Salmon and RSEM; featureCounts appears only as
  a QC tool (biotype QC; `--featurecounts_group_type`, `--featurecounts_feature_type`).
- UMI handling is opt-in: `--with_umi`, `--umi_dedup_tool umitools|umicollapse`,
  `--umitools_extract_method`, `--umitools_bc_pattern[-2]`,
  `--umitools_grouping_method`, `--skip_umi_extract`.
- BAM publication via `--save_align_intermeds`; index reuse/build via
  `--save_reference`; unaligned reads via `--save_unaligned`.

**nf-core/differentialabundance 2.0.0**:

- `--differential_method` accepts `deseq2` (DEFAULT), `limma` (+voom),
  `dream`, `propd`; edgeR is not offered (the legacy workflow is DESeq2-only,
  so nothing is lost). DESeq2 is the ratified backend (D2).
- Contrasts come from a `--contrasts` CSV/TSV/YAML file; blocking factors
  come from samplesheet columns. Residual WP2 verification: DESeq2
  design-covariate (batch) wiring, contrast semantics, matrix input
  contract, container set - gated before adoption is finalized.
- Functional analysis: `--functional_method gsea|gprofiler2|decoupler|grea`
  with GMT gene-set files. OrgDb-based GO/KEGG over-representation (the
  legacy clusterProfiler chain) is NOT provided - D3 stays a local stage.
- Inputs: `--matrix` (TSV/CSV) + samplesheet; DESeq2 knobs exposed (vst/rlog,
  LFC shrink, test/fit type, pre-filtering).

**Module availability** (nf-core/modules): `featurecounts/featurecounts` and
`rseqc/inferexperiment` exist (used by the local featureCounts stage and its
optional strand inference); no clusterProfiler module exists (enrichment is a
local module).

## 4. Preferred current methods (Q4)

- Alignment: native STAR (two-pass) on the adopted route; BAMs published for
  the extension stages and the lncRNA handoff.
- Quantification: **switchable per project** (D1) - the native Salmon gene
  matrix OR a local featureCounts stage on the published STAR BAMs (legacy
  counting semantics). The switch selects the authoritative matrix for the
  differential route, comparisons, and the lncRNA handoff.
- Differential expression: DESeq2 via nf-core/differentialabundance 2.0.0
  (its default backend), gated by the WP2 verification above; fallback is a
  local DESeq2 extension (chipseq DiffBind pattern).
- Enrichment: local clusterProfiler stage wrapping the validated legacy R
  chain; OrgDb supplied by configuration (package name or installable path
  staged as a process input) - the pipeline never binds a species.
- UMI: native `--with_umi` (umitools extract + dedup) when a project opts in;
  default off (most projects have no UMI). Replaces the legacy clip-only
  semantics (8 bp random UMI + T anchor at position 9 configured per project).

## 5. Disposition of the legacy capabilities (Q5)

All dispositions ratified by the owner (rounds 1+2, 2026-10-11).

| # | Legacy capability | Disposition (RATIFIED) | Rationale |
|---|---|---|---|
| 1 | QC + trimming + MultiQC | **replaced** (native) | adopted pipeline covers; parameter/reference configuration only |
| 2 | STAR two-pass alignment | **adopted** (native) | `--save_align_intermeds` publishes the BAMs the ext and the lncRNA handoff consume |
| 3 | featureCounts quantification | **retained as local ext stage, switchable** | no native STAR+featureCounts route in 3.27.0; ext wires the ecosystem module on published BAMs; D1 switch picks the authoritative matrix |
| 4 | TPM/FPKM matrices | **replaced** (TPM from the native route; FPKM dropped) | scientific difference documented, not eliminated |
| 5 | Strandedness inference | **verify native; else ext helper or per-project param** | the featureCounts stage needs a strand setting: RSeQC inference (ecosystem module) or explicit config - never hardcoded (D7) |
| 6 | DESeq2 DE with batch blocking | **adopted** (pending WP2 gate) | DESeq2 is differentialabundance's default backend (D2); fallback = local DESeq2 ext |
| 7 | clusterProfiler GO/KEGG + GSEA | **retained as local ext stage** | OrgDb is a config option; no species binding (D3) |
| 8 | DEGgroupCompare (between-group comparison) | **dropped** | owner round 2; recorded as a deliberately dropped capability for legacy-retirement documentation |
| 9 | StringTie merge + gffcompare + isoform quant | **deferred into the lncRNA workflow** | owner round 1; Stage 2 publishes the BAMs only (D5) |
| 10 | UMI clip-in-STAR (no dedup) | **replaced by native `--with_umi`** | real dedup replaces clip-only; opt-in per project (D6) |

## 6. Local code still necessary (Q6)

Ratified local surface (all default-off, mirroring the legacy flag
discipline):

1. `rnaseq-ext` pipeline (`nextflow/pipelines/rnaseq-ext/`) with stages:
   - `featurecounts_quant` - ecosystem featureCounts module on the published
     STAR BAMs; strand `auto` (RSeQC infer_experiment, legacy mapping) or an
     explicit per-project setting; consumes the D1 route switch.
   - `enrichment` - local clusterProfiler module on the differential results;
     OrgDb/KEGG organism keys are params; offline KEGG behavior documented.
   - `deseq2` - ONLY on the D2 fallback (WP2 gate failure).
   - An input-preparation script for differentialabundance (matrix +
     samplesheet + contrasts from the switch-selected matrix) - a `bin/`
     helper, not a stage.
2. Thin adoption entries: `nextflow/adopted/rnaseq/` and
   `nextflow/adopted/differentialabundance/` (the latter lands only if the
   WP2 gate passes; otherwise the fallback stage above replaces it).

NOT built: `deg_compare` (D4), `stringtie_merge` (D5).

## 7. Validation demonstrating scientific usability (Q7)

- Baseline: test-profile smoke on the HPC (throttled local executor,
  apptainer profile), then a PBS executor run - the Stage-1 validation shape.
- Real-data round (WP7; owner inputs pending): the adopted pipeline + ext
  stages on the owner-named dataset(s), then a targeted legacy comparison
  per design Section 10.3 - DEG sign/rank agreement on the same contrast and
  top-enriched-term agreement. Differences understood and documented, not
  eliminated.
- The D1 switch is exercised in both settings (salmon route and
  featureCounts route) so both matrices are validated.

## 8. Decision log

Ratified by the owner in two review rounds on 2026-10-11. Standing
principles ratified with the items below: **the pipeline never binds a
species** (references, OrgDb, KEGG organism codes are all configuration), and
**nothing is silently hardcoded** - every analysis knob is an overridable
parameter with a documented default (echoes design D4).

| ID | Topic | Resolution (RATIFIED) | Round |
|---|---|---|---|
| D1 | Execution route | Adopted run uses `--aligner star_salmon` with `--save_align_intermeds` (STAR BAMs published). Quantification is switchable per project: native Salmon matrix OR the local featureCounts stage on the published BAMs; the switch selects the authoritative matrix for DA, comparisons, and the lncRNA handoff. Both matrices may exist on disk; the README documents which is authoritative per setting. | 2 |
| D2 | DE statistics route | DESeq2 required; differentialabundance 2.0.0 has `--differential_method deseq2` as the DEFAULT backend (no edgeR; legacy is DESeq2-only). Adoption proceeds if the WP2 gate passes (DESeq2 design covariates, contrast semantics, matrix contract, containers); fallback = local DESeq2 ext (chipseq DiffBind pattern). | 1 |
| D3 | Enrichment route | Local clusterProfiler stage; OrgDb is a config option (package name or installable path staged as input) - no species binding. Offline KEGG behavior documented. | 1 |
| D4 | Between-group DEG comparison | Dropped - no local stage; capability retired with this route. | 2 |
| D5 | Transcript assembly | Moves into the lncRNA workflow; deferred. Stage 2 publishes BAMs only. | 1 |
| D6 | UMI handling | Native `--with_umi` route (umitools extract + dedup); opt-in, default-off per project; patterns configured from the verified data basis when enabled. | 1 |
| D7 | Parameter policy | Native defaults as documented baseline; every parameter exposed and overridable per project; no silent absolute configuration. | 1 |
| D8 | Tracking of the plan document | Stays local under `plans/` (untracked, repo convention). | 1 |

## 9. Version pin

- Pipeline: nf-core/rnaseq 3.27.0, gated by the WP0 Nextflow 26.04.6 parse
  check (dev-HEAD fallback per Stage-1 policy); the final pin and its
  validation date are recorded in `nextflow/versions.yml` only after a
  validated run.
- Containers: per-step modular images (design doc D6) seeded per the
  validated acquisition recipe (implementation guidance Section 3.1); the
  local enrichment module declares its own container at implementation time.
- Module installs record their gitSha in `modules.json` at implementation
  time; no versions are frozen in the mainline (design doc D7).
