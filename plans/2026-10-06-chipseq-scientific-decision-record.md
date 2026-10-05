# nf-core/chipseq scientific rebaseline record

- **Date**: 2026-10-06
- **Status**: RATIFIED. The owner resolved all four decision points on
  2026-10-06: every retained legacy capability is delivered by upgrading the
  existing local implementation as a local extension on the adopted route; no
  replacement alternatives are taken (decision log, Section 8). Implementation
  specifications are written against this frozen version.
- **Mandate**: the per-workflow scientific rebaseline record required by the
  migration design doc, Section 9. This record precedes any implementation
  specification for the chipseq Stage-1 route.
- **Pipeline under review**: nf-core/chipseq dev HEAD
  `1b3078d71c4c8126a9b976d3b0248e9b9b525608` (the 2.1.0 release is incompatible
  with Nextflow 26; see implementation guidance Section 3.1).

---

## 1. Biological use case (Q1)

ChIP-seq and CUT&Tag analyses of rice (*Oryza sativa*) and other species the
user configures locally: histone-mark and transcription-factor binding profiles,
replicate-aware peak sets, signal tracks, QC reporting, and differential
binding between conditions. Reference genomes are never fixed in the pipeline;
species presets live in machine-local/project config (design doc, D4 and
Section 15 amendment 3).

## 2. Inputs and scientifically meaningful outputs (Q2)

- **In**: raw FASTQ (PE/SE), reference FASTA + GTF + BWA/Bowtie2 indices,
  optional control (IgG/input) per sample, optional spike-in genome FASTA,
  sample sheet with group/antibody/control/batch metadata.
- **Out**: per-sample and consensus peak BED/narrowPeak, annotated peak tables,
  per-sample normalized bigWig signal tracks, QC report (MultiQC + NSC/RSC +
  fingerprint + correlation), replicate-support annotations, differential
  binding tables (contrasts with FDR/fold-change), motif enrichment reports.

## 3. What the ecosystem provides (Q3)

Verified against the pinned dev HEAD checkout on 2026-10-05/06 and a live
`nf-core modules list remote` query on 2026-10-06 (https://nf-co.re/modules).

**Native dev HEAD model** (read from `workflows/chipseq.nf` and the local
subworkflows, not from docs):

- Per-sample alignment (BWA / Bowtie2 / Chromap / STAR), picard markduplicates,
  bamtools filtering (MAPQ + blacklist).
- **Per-sample MACS3 peak calling** — each sample (with its own control when
  declared) is called individually; there is no pooled-per-group call and no
  separate per-replicate-vs-pooled two-tier scheme.
- Cross-replicate handling: `MACS3_CONSENSUS` (custom local script) merges peak
  sets per antibody into a consensus BED + SAF, and
  `ANNOTATE_BOOLEAN_PEAKS` adds a per-sample presence boolean to every
  consensus peak — replicate support is carried as annotation, not as an IDR
  statistic. No IDR anywhere in the pipeline.
- Consensus quantification: featureCounts over all samples on the consensus
  SAF, then `DESEQ2_QC` — PCA and sample-distance QC plots from the count
  matrix. **This is QC only; no contrast testing, no differential binding.**
- Signal tracks: per-sample bigWigs via bedtools genomecov (CPM scaling) +
  UCSC bedGraphToBigWig; deepTools computematrix/plotProfile/plotHeatmap/
  plotFingerprint suite; IGV sessions.
- QC: FastQC, phantompeakqualtools NSC/RSC, preseq complexity, picard metrics,
  HOMER annotatepeaks (peak annotation only — no motif discovery),
  MultiQC aggregation.

**Module availability in nf-core/modules (live-verified)**:

| Gap tool | Module present? |
|---|---|
| SEACR | yes — `seacr/callpeak` |
| TOBIAS | yes — `tobias/atacorrect`, `tobias/scorebigwig`, `tobias/bindetect` |
| IDR | yes — `idr` |
| MACS2 | yes — `macs2/callpeak` |
| HOMER motif discovery | **no** — only `homer/{annotatepeaks,findpeaks,maketagdirectory,makeucscfile,pos2bed}` |
| DiffBind (R/Bioconductor) | **no** |
| csaw | **no** |

**Related pipelines**:

- nf-core/cutandrun natively supports spike-in normalization, controls, and
  SEACR/MACS2 calling for CUT&RUN/CUT&Tag data (design doc, Section 5.6).
- nf-core/atacseq ships TOBIAS footprinting natively as an optional
  subworkflow (`--with_tobias`: ATACorrect -> scoreBigwig -> BINDetect),
  verified 2026-10-06.
- No nf-core ChIP-seq pipeline performs differential binding contrasts.

## 4. Preferred current methods (Q4)

- Peak calling: MACS3 (upstream default) over legacy MACS2 — maintained,
  same statistical model, actively versioned.
- Reproducibility: native consensus (with per-sample support annotation) as
  the base, plus a local IDR stage for narrow peaks (ratified, D3) — consensus
  is upstream-maintained, IDR restores the ENCODE-style reproducibility
  statistic the legacy pipeline used for narrow peaks.
- Signal tracks: bedtools+UCSC CPM bigWigs (native) over deeptools bamCoverage
  RPGC (legacy) — fewer moving parts; fold-enrichment style tracks remain
  available through MACS3's own outputs if a project needs them (scientific
  note: normalization semantics differ from the legacy default).
- Differential binding: DiffBind remains the scientifically preferred method
  (summit-window counting, DESeq2/edgeR, batch blocking) and must be locally
  provided (decision D1).
- Annotation: HOMER annotatepeaks (native) replaces the legacy
  ChIPseeker/annot step; motif discovery has no ecosystem module and is
  provided by a local wrapper (ratified, D4).

## 5. Disposition of the nine legacy functions (Q5)

All dispositions ratified by the owner on 2026-10-06.

| # | Legacy capability | Disposition (RATIFIED) | Rationale |
|---|---|---|---|
| 1 | Replicate-aware peak calling (v0.5) | **replaced** | dev HEAD calls peaks per sample natively; replicate support surfaces in the consensus boolean matrix |
| 2 | IDR / overlap consensus | **expanded** — native consensus + local IDR extension (D3) | native consensus is retained (upstream-maintained); a local IDR stage restores ENCODE-style reproducibility for narrow peaks (legacy v0.5 parity) |
| 3 | Spike-in normalization | **retained as local extension** (D2) | no native support on the chipseq route; port legacy co-alignment + scaling logic |
| 4 | SEACR caller | **retained as local extension** — wire the ecosystem `seacr/callpeak` module as the alternative caller for CUT&Tag samples (D2) | CUT&Tag data routes through the chipseq route per the owner's no-alternatives decision; the module exists, so this is wiring plus config, not new science; SEACR stays irrelevant for ChIP-seq |
| 5 | DiffBind differential binding | **retained as local module** (D1) | ecosystem vacuum; legacy R logic is validated and scientifically preferred |
| 6 | QC gate table (v0.6) | **retained as local extension** (D4) | no native equivalent; cheap to port (PASS/WARN/FAIL from existing QC outputs) |
| 7 | TOBIAS footprinting | **replaced** (moves to ATAC route) | `tobias/*` modules verified available and nf-core/atacseq ships TOBIAS natively (`--with_tobias`); install in the atacseq route, not here |
| 8 | HOMER motif discovery | **retained as local wrapper** (D4) | no `findmotifs` module in the ecosystem; native annotatepeaks covers annotation only |
| 9 | Per-sample bigWigs | **replaced** | native per-sample bigWigs (CPM) + full deepTools suite; normalization semantics differ from legacy FE/RPGC (documented, accepted) |

Additional local capabilities outside the nine:

- **Organelle QC** (chloroplast/mitochondrial read fraction): **retained as
  local extension** — rice-critical, no native equivalent, trivial to port.
- TSS enrichment: covered by native deepTools computematrix over a gene BED;
  no separate port needed.

## 6. Local code still necessary (Q6)

Ratified local surface for the chipseq route (all default-off, mirroring the
legacy flag discipline):

1. `diffbind` — local nf-core module wrapping the legacy R analysis
   (DESeq2/edgeR, summit+/-flank windows, contrasts, optional batch blocking),
   own container built per the Section 3.1 recipe.
2. `spike-in` — local extension: spike-in genome co-alignment, scaling-factor
   computation, optional bigWig rescaling; reuses native aligner containers.
3. `idr` — local IDR stage for narrow peaks: install the ecosystem `idr`
   module, feed it the native per-replicate narrowPeak calls, configurable
   threshold (legacy default 0.05); broad peaks stay on native consensus.
4. `seacr` — wire the ecosystem `seacr/callpeak` module as the alternative
   peak caller for CUT&Tag samples (caller selection by config, mirroring the
   legacy `peak.caller` semantics).
5. `qc-gates` — local script emitting the PASS/WARN/FAIL gate table from
   flagstat/picard/MultiQC JSON outputs.
6. `homer-findmotifs` — local module wrapping `findMotifsGenome.pl`
   (genome tag from config).
7. `organelle-qc` — local script (idxstats fraction against configurable
   organelle patterns).

None of these blocks the pure-adoption baseline; each lands as its own
default-off stage.

## 7. Validation demonstrating scientific usability (Q7)

- Baseline (already done): test-profile dry-run + throttled smoke run on the
  HPC head node, 29 containers seeded, outputs verified (BASELINE.yml).
- PBS executor end-to-end run before real-data work (workq was full during
  seeding).
- Per retained extension, one real-data check against the corresponding legacy
  result on the same input: DiffBind contrasts (sign direction + rank
  agreement), spike-in scaling (expected ChIP/input ratio recovery), QC gate
  table (row agreement with legacy gates on a legacy project), HOMER motifs
  (top motifs agree on the same peak set).
- Legacy comparison stays targeted per design doc Section 10.3: no parity
  acceptance on filenames/paths/defaults.

## 8. Decision log

Ratified by the owner on 2026-10-06. Guiding principle for every item:
**upgrade the existing local implementation as a local extension; no
replacement alternatives.**

| ID | Topic | Resolution (RATIFIED) | Rejected alternatives |
|---|---|---|---|
| D1 | Differential binding | Port DiffBind as a local nf-core module + container (Section 6 item 1) | DESeq2 on native consensus counts; defer to post-adoption |
| D2 | Spike-in + CUT&Tag routing | Local spike-in extension on the chipseq route, used by ChIP and CUT&Tag data; implies local SEACR caller wiring (Section 6 item 4) | Route CUT&Tag data to nf-core/cutandrun; drop spike-in until needed |
| D3 | IDR vs consensus | Consensus + IDR dual track: native consensus plus a local IDR stage for narrow peaks via the ecosystem `idr` module (Section 6 item 3) | Native consensus only |
| D4 | Local extras (HOMER motifs, QC gates, organelle QC) | Keep all three as local extensions at v1 (Section 6 items 5-7) | Drop HOMER motifs; pure adoption with no extras at v1 |

The record is frozen as of the ratification date; the Stage 1 implementation
specification is written against this version only.

## 9. Version pin

- Pipeline: nf-core/chipseq dev HEAD `1b3078d71c4c8126a9b976d3b0248e9b9b525608`.
- Containers: 29-image runtime set already seeded machine-locally (2026-10-05,
  BASELINE.yml); local extensions add their own containers at implementation
  time via the validated acquisition recipe (design doc Section 15 amendment 6,
  guidance doc Section 3.1).
- Module installs (`nf-core modules install`) record their gitSha in
  `modules.json` at implementation time; no versions are frozen in the
  mainline (design doc D7).
