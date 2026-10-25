# nf-core/atacseq scientific rebaseline record

- **Date**: 2026-10-12
- **Status**: RATIFIED. The owner instructed implementation per the Stage-3
  plan (`plans/2026-10-12-stage3-atac-faire-migration-plan.md`, machine-local);
  its D1-D8 recommendations stand as ratified, with one evidence-driven
  correction recorded in the decision log (D2: nf-core/atacseq ships NO
  TOBIAS - the design amendment 7 `--with_tobias` premise is wrong; the plan's
  own D2 override path applies: a default-off ext port via the central
  nf-core tobias modules). Real-data rounds are skipped by owner decision
  (design Section 15 amendment 8); execution validation (test-profile smoke,
  real-container fixture execution) is deferred with the stage (WP6).
- **Mandate**: the per-workflow scientific rebaseline record required by the
  migration design doc, Section 9. This record precedes the Stage-3
  implementation commits.
- **Pipeline under review**: nf-core/atacseq dev HEAD
  `6a9307709662bad3e7bcc10aa873f7e81fdd1646` (v2.2.0dev, 2026-09-30,
  `nextflowVersion = '!>=25.10.4'`). Pin gate 2026-10-12 under Nextflow
  26.04.6: tag 2.1.2 FAILS to parse (`Unable to parse config file` - same
  old-template failure class as chipseq 2.1.0); dev HEAD parses
  (`nextflow config -flat -profile test` exit 0).

---

## 1. Biological use case (Q1)

ATAC-seq and FAIRE-seq analyses of rice (*Oryza sativa*) and other locally
configured species: open-chromatin accessibility profiles, replicate-aware
peak sets, signal tracks, insertion/fragment QC, organelle contamination
diagnostics for plant libraries, Tn5-bias-corrected footprinting, and
differential-accessibility-ready peak quantification. References are never
fixed in the pipeline (design D4; amendment 3).

## 2. Inputs and scientifically meaningful outputs (Q2)

- **In**: raw FASTQ (PE only; the legacy workflow enforces PE and this route
  inherits that limit - SE is documented as dropped), reference FASTA + GTF +
  BWA/Bowtie2/Chromap indices (igenomes or local paths), optional control per
  group, blacklist file, optional combined spike-in reference.
- **Out (adopted route)**: filtered markdup BAMs `{id}.mLb.clN.sorted.bam`
  (published unconditionally), Tn5-shifted BAMs, MACS3 narrow/broad peaks +
  boolean-support consensus peaks, consensus-site quantification with native
  DESeq2 QC, bigWigs (FE + shifted), FRiP scores, ataqv/insertion metrics,
  deepTools profiles/fingerprints (TSS enrichment among them), MultiQC.
- **Out (ext stages, all default-off)**: organelle fraction table, QC gate
  table, spike-in scale factors + rescaled bigWigs, TOBIAS footprinting chain
  outputs (corrected tracks, footprint scores, BINDetect results).

## 3. What the ecosystem provides (Q3)

nf-core/atacseq dev HEAD (verified against the pinned source 2026-10-12, not
docs): BWA/Bowtie2/Chromap alignment with `--save_unaligned` support, Picard
MarkDuplicates + filter chain, native Tn5 shifting (BAM_SHIFT_READS via
deepTools alignmentsieve), MACS3 per-sample + `MACS3_CONSENSUS` boolean
support matrix with `--min_replicates`, FRiP, ataqv/mkarv, deepTools
computeMatrix/plotProfile/plotHeatmap/plotFingerprint, HOMER annotatePeaks,
GFFREAD/GTF2BED, preseq, khmer unique kmers, MultiQC, and a native
consensus-site DESeq2 QC route (BED_CONSENSUS_QUANTIFY...DESEQ2).

Publish contract consumed by the ext (pinned `conf/modules.config`):

| Artifact | Published location |
|---|---|
| Filtered BAM (+bai) | `{outdir}/{aligner}/merged_library/{id}.mLb.clN.sorted.bam(+.bai)` - unconditional |
| Stats/flagstat/idxstats | `{outdir}/{aligner}/merged_library/samtools_stats/{id}.mLb.clN.sorted.bam.{stats,flagstat,idxstats}` (clN; a parallel mkD-prefixed set publishes from the markdup BAM) |
| Picard metrics | `{outdir}/{aligner}/merged_library/picard_metrics/{id}.mLb.mkD.sorted*` |
| MACS3 peaks | `{outdir}/{aligner}/merged_library/macs3/{narrow_peak|broad_peak}/{id}_peaks.narrowPeak|broadPeak` |
| Consensus peaks | `{outdir}/{aligner}/merged_library/macs3/{type}/consensus/consensus_peaks.mLb.clN*` |
| bigWig | `{outdir}/{aligner}/merged_library/bigwig/{id}.mLb.clN*.bigWig` |
| Unmapped FASTQs | `{outdir}/{aligner}/library/unmapped/*.fastq.gz` - only with `--save_unaligned` (default false) |

Param defaults that shape routing: `aligner = 'bwa'` (bowtie2 available),
`narrow_peak = false` (broad peaks are the ATAC default), `min_reps_consensus
= 1`, `save_unaligned = false`.

NOT provided by the ecosystem route (verified by source inspection + a
full-source grep on both 2.1.2 and dev): **TOBIAS footprinting** - there is
no `--with_tobias` flag and no tobias module anywhere in the pipeline (the
design amendment 7 claim was wrong for every released/dev revision). Central
nf-core modules DO carry tobias/atacorrect, tobias/scorebigwig, and
tobias/bindetect (container `biocontainers/tobias:0.17.5--py310h3479294_0`).

FAIRE-seq: no dedicated maintained community/nf-core workflow exists (design
5.10 finding re-checked 2026-10-12, unchanged). The legacy BioWorkflows
workflow treats FAIRE identically to ATAC (same rule, same dedup default,
narrowPeak); the community norm is open-chromatin peak calling without
Tn5-specific processing.

## 4. Preferred current methods (Q4)

- **Route**: adopt nf-core/atacseq dev HEAD (thin entry, no fork) with
  Bowtie2 as the documented plant-route aligner (end-to-end very-sensitive
  PE, matching the legacy alignment policy; BWA remains the upstream default
  and stays selectable). Broad peaks by default (upstream default; the
  legacy default narrow for ATAC is superseded by the ENCODE-oriented
  broad-for-accessibility convention of the adopted route - a per-project
  switch, documented in the adopted params example).
- **Peak reproducibility**: native MACS3_CONSENSUS boolean support matrix
  replaces the legacy classic-IDR route for ATAC/FAIRE (recorded
  replacement; the classic idr module remains available on the chipseq
  route where ChIP narrow peaks are the norm). FAIRE runs the same entry.
- **Footprinting**: TOBIAS via the central nf-core modules, ported as a
  default-off ext stage (D2 override; see Section 8).
- **Peak annotation**: native HOMER annotatePeaks replaces legacy ChIPseeker
  annotation for the ATAC route (recorded replacement, design 10.3: tools
  changing is not a failure).

## 5. Disposition of the legacy ATAC/FAIRE capabilities (Q5)

| # | Legacy capability (chip_cuttag_atac_faire) | Disposition |
|---|---|---|
| 1 | trim_galore -q 25 trimming | adopted (native trimgalore) |
| 2 | bowtie2 end-to-end very-sensitive PE, mapq >= 30 | adopted (native bowtie2; tuning documented in adopted params) |
| 3 | picard dedup for atac/faire | adopted (native markduplicates + filter chain) |
| 4 | MACS2 ATAC calling: bampe or shifted (Tn5 -100/200) | adopted (native MACS3 + native shift; bampe/shifted map to upstream defaults) |
| 5 | FE + per-sample bigWigs | adopted (native bedgraph/bigwig chain) |
| 6 | FRiP (per group, sample) | adopted (native FRIP_SCORE) |
| 7 | deepTools QC incl. fragment sizes | adopted (native plotFingerprint/computeMatrix/plots; ataqv covers insertion metrics) |
| 8 | TSS enrichment score | adopted (native reference-point TSS matrix + profile) |
| 9 | ChIPseeker peak annotation | replaced (native HOMER annotatePeaks) |
| 10 | Replicate IDR / consensus (atac/faire as narrow) | replaced (native MACS3_CONSENSUS support matrix) |
| 11 | QC gate summary table (v0.6, default off) | ext port (atacseq-ext qc_gates) |
| 12 | Organelle fraction QC (v0.6, default off) | ext port (atacseq-ext organelle_qc) |
| 13 | Spike-in scaling via unmapped-pair realignment (v0.6, default off) | ext port, combined-reference model (atacseq-ext spike_in; identical to the chipseq-ext mechanism - the adopted run is configured with a combined host+spike reference and the ext recomputes idxstats from the published BAMs and reads spike contigs from them) |
| 14 | TOBIAS ATACorrect/ScoreBigwig/BINDetect (v0.6, default off) | ext port via central nf-core tobias modules (atacseq-ext footprint; D2 override) |
| 15 | SEACR alternative caller | not ported (CUT&Tag/RUN-oriented; lives on the chipseq route) |
| 16 | DiffBind differential binding | not ported in this stage (recorded as a future default-off candidate against the native consensus-peak quantification) |
| 17 | SE input support | dropped (legacy never supported it; documented at retirement) |
| 18 | HMMRATAC / FSeq2 callers | dropped (existed only in the archived v0.1.0 snapshot) |

## 6. Local code still necessary (Q6)

`nextflow/pipelines/atacseq-ext/` (thin, default-off stages only; cloned from
the chipseq-ext conventions):

- `organelle_qc_ext` - samtools idxstats + organelle fraction summary
  (plant ATAC/mito diagnostic; ported as-is from chipseq-ext).
- `qc_gates_ext` - QC gate table over published flagstats + picard metrics +
  organelle summary (ported as-is; thresholds are config, not code).
- `spikein_scaling` - spike scale factors + optional rescaled bigWigs
  (ported as-is; combined-reference model).
- `footprint_ext` - per-group pooled-BAM TOBIAS chain: group merge ->
  TOBIAS_ATACCORRECT -> TOBIAS_SCOREBIGWIG -> TOBIAS_BINDETECT (new port of
  the legacy v0.6 stage, using the central nf-core tobias modules; default
  off; known legacy trap carried as a note: BINDetect can exit on sklearn
  NaN with degenerate synthetic signal - re-validate on real research data).

Everything else is adopted upstream code - no forking, no re-implementation.

## 7. Validation demonstrating scientific usability (Q7)

Per owner decision (design Section 15 amendment 8, 2026-10-12): real-data
rounds are skipped; validation evidence is

- the WP0 pin gate (this record, Section header);
- synthetic fixtures mirroring the pinned publish contract + dry-run
  baselines (all-off = 0 tasks; all-on = N processes) recorded in
  `nextflow/docs/atacseq-ext-validation-notes.md`;
- pytest unit suites for every ext helper script;
- DEFERRED execution validation (recorded in the validation notes when run):
  adopted test-profile smoke on the HPC + real-container fixture execution
  of the ext stages (Stage-1 model). versions.yml rows follow the ledger
  rule - only after those runs.

## 8. Decision log

- **D1 pin** - 2.1.2 fails the Nextflow 26 parse gate; dev HEAD
  `6a9307709662bad3e7bcc10aa873f7e81fdd1646` passes. Pin dev HEAD; re-run
  the consumed-path contract audit on any re-pin (Stage-1 watch-item rule).
- **D2 footprinting route** - CORRECTION: the design amendment 7 premise
  ("nf-core/atacseq covers ATAC-seq natively, including TOBIAS footprinting
  `--with_tobias`") is false for every released and dev revision (verified
  2026-10-12 by full-source grep). Resolution: the Stage-3 plan's own D2
  override path applies - the legacy TOBIAS stage is ported as a default-off
  atacseq-ext stage using the central nf-core tobias modules. No silent
  capability drop.
- **D3 ext scope** - four default-off stages: organelle_qc, qc_gates,
  spike_in, footprint (the first three are chipseq-ext clones; the fourth is
  the D2 override port). SEACR and DiffBind are NOT ported (15/16 above).
- **D4 FAIRE routing** - FAIRE runs through the SAME adopted atacseq entry
  with a documented FAIRE settings block in the adopted params example
  (control optional; TSS interpreted as a sonication/insertion-bias check
  for FAIRE; footprinting OFF for FAIRE). No separate pipeline.
- **D5 UMI** - CORRECTION: the pinned workflow hardcodes UMI extraction OFF
  (the shared trimgalore subworkflow supports it, but the main workflow
  passes `with_umi = false` and exposes no opt-in flag). UMI-based dedup is
  also scientifically marginal for ATAC, where PCR-duplicate removal is the
  standard. Disposition: not part of the adopted surface; the params example
  documents `clip_r1/clip_r2` for trimming-only use. No ext port (the legacy
  workflow never had ATAC UMI support either).
- **D6 retirement mechanics** (chip_cuttag_atac_faire, design Stage 4) -
  DEFERRED with owner sign-off pending: tag the final legacy state, write
  the dropped-capabilities note, archive the directory outside the repo,
  update root README/AGENTS. Runs only after the owner confirms no active
  legacy project.
- **D7 assay routing table** - chip/cuttag -> chipseq route (+chipseq-ext);
  atac/faire -> atacseq route (+atacseq-ext). Recorded here and in the root
  README migration section.
- **D8 plan-document tracking** - the Stage-3 plan stays local under
  `plans/` (untracked, repo convention).

## 9. Version pin

- nf-core/atacseq: dev HEAD `6a9307709662bad3e7bcc10aa873f7e81fdd1646`
  (v2.2.0dev, 2026-09-30). `nextflowVersion !>=25.10.4` - compatible with
  the repo runtime Nextflow 26.04.6. Gate evidence: `nextflow config -flat
  -profile test` exit 0 on 2026-10-12 (2.1.2 fails: "Unable to parse config
  file", the chipseq-2.1.0 failure class).
- versions.yml row: added only after the adopted test-profile smoke
  (deferred; ledger rule).
