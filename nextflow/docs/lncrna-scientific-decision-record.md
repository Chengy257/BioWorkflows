# BioWorkflows lncRNA identification workflow scientific rebaseline record

- **Date**: 2026-10-12
- **Status**: RATIFIED. The owner instructed implementation per the Stage-4
  plan (`plans/2026-10-12-stage4-lncrna-migration-plan.md`, machine-local);
  its D1-D9 recommendations stand as ratified. Real-data rounds are skipped
  by owner decision (design Section 15 amendment 8); execution validation
  (real-container fixture run) is deferred with the stage (WP6).
- **Mandate**: the per-workflow scientific rebaseline record required by the
  migration design doc, Section 9. This record precedes the Stage-4
  implementation commits.
- **Route category**: design-doc Stage 3 (local DSL2 build) - the review
  below proves no adoptable upstream exists.

---

## 1. Biological use case (Q1)

De novo long-noncoding-RNA identification from RNA-seq alignments of rice
(*Oryza sativa*) and other locally configured species: reference-guided
transcript assembly, novel (intergenic, classcode-u) transcript extraction
against a protein-coding reference, coding-potential + length filtering,
protein-domain (Pfam) and homology (NR/DIAMOND) filters, the final lncRNA
set, and its expression matrix against the reference-plus-novel GTF.

## 2. Inputs and scientifically meaningful outputs (Q2)

- **In** (the Stage-2 handoff, ratified in the rnaseq decision record D1/D5):
  published markdup STAR BAMs
  `{rnaseq_outdir}/{aligner}/{sample}.markdup.sorted.bam(+.bai)`
  (`bam_suffix` switchable to the pre-dedup `{sample}.sorted.bam` when the
  adopted run passed `--save_align_intermeds`), the reference genome FASTA,
  the annotation GTF, and the protein-coding CDS GTF (`gtf_pcg`; for rice
  the user-built MSU+RAP combined CDS GTF - build recipe documented, D7).
- **Out**: per-sample assemblies, the merged novel-transcript set
  (classcode u), the final lncRNA FASTA/GTF, `ref_withLncRNA.gtf`, and the
  lncRNA-level count/TPM matrices (`lncrna.count.matrix.tsv`,
  `lncrna.tpm.matrix.tsv` - named distinctly from the Stage-2 gene-level
  matrices).

## 3. What the ecosystem provides (Q3)

- **nf-core/lncpipe**: still "UNDER DEVELOPMENT" (dockstore v0.2.4; no
  stable release; absent from recent nf-core release newsletters as of
  2026-10-12). NOT adoptable; watch item (re-review if it releases).
- **Central nf-core modules verified present** (API check 2026-10-12):
  `stringtie` (assembly; container stringtie:3.0.3), `stringtie/merge`,
  `gffcompare`, `gffread`, `cpc2`, `diamond/blastx` (+ `diamond/makedb`),
  `subread/featurecounts`. NOT present: `pfamscan` (legacy pfam_scan.pl),
  `cpat`.
- **Containers verified on quay.io/biocontainers** (2026-10-12):
  `cpat:3.0.5--py310h20b60a1_4`, `pfam_scan:1.6--hdfd78af_5`.
- **Coding-potential tool status**: CNCI is a python2-only, effectively
  unmaintained repository (its ML successor is CNIT); CPC2 is
  maintenance-frozen but packaged and callable (central module exists);
  CPAT 2.x/3.x is python3, bioconda-packaged, actively maintained.

## 4. Preferred current methods (Q4)

- **Execution model**: local DSL2 pipeline `nextflow/pipelines/lncrna-nf/`
  consuming the Stage-2 handoff; modules-first (central modules vendored
  verbatim; local modules in nf-core module format for the true gaps).
- **Assembly**: StringTie per-sample `--conservative -j 5` (legacy
  parameters, validated on rice projects) + `stringtie --merge` against
  `gtf_pcg`, then `gffcompare` classcode-u novel extraction and `gffread`
  FASTA extraction - unchanged scientifically, containerized.
- **Coding potential** (the scientific upgrade point, plan D3): CPC2 AND
  CPAT 3.0.5 intersection + length >= 200 nt, replacing the legacy
  CPC2 AND CNCI AND length rule (CNCI retired with the evidence recorded
  here; two independent ORF/ML approaches, both container-clean).
- **Homology filters**: pfam_scan (`pfam_scan.pl -translate`, E < 1e-5) and
  DIAMOND blastx vs `nr.dmnd` (any-hit removal) - kept, containerized,
  default-ON with per-project OFF switches.
- **Expression**: subread/featureCounts (vendored from rnaseq-ext's copy)
  against `ref_withLncRNA.gtf`; strand mode from the project config.
- **Species neutrality**: genome/gtf/gtf_pcg/org context are params only;
  rice is the documented example configuration (design D4, amendment 3).

## 5. Disposition of the legacy lncrna-mode capabilities (Q5)

| # | Legacy capability (rna-seq lncrna mode) | Disposition |
|---|---|---|
| 1 | Per-sample StringTie `--conservative -j 5` assembly | kept (central stringtie module) |
| 2 | `stringtie --merge -G gtf_PcGs -F 0.1 -T 0.1` + gffcompare classcode-u + gffread FASTA | kept (central stringtie/merge, gffcompare, gffread modules + a local classcode-u helper) |
| 3 | CPC2 noncoding prediction | kept (central cpc2 module) |
| 4 | CNCI noncoding prediction (python2) | replaced (CPAT 3.0.5; CNCI retired, python2 container not maintainable - Q3 evidence) |
| 5 | length >= 200 nt + intersection of predictions | kept (local helper; now CPC2 AND CPAT AND length) |
| 6 | pfam_scan.pl -translate, E < 1e-5 | kept (local module, container pfam_scan:1.6) |
| 7 | DIAMOND blastx vs nr.dmnd, any-hit removal | kept (central diamond/blastx module) |
| 8 | final set (Pfam U NR removed) + ref_withLncRNA.gtf | kept (local helper) |
| 9 | featureCounts vs ref_withLncRNA + TPM/FPKM matrices | kept (vendored subread/featurecounts + the rnaseq-ext merge helper) |
| 10 | `as`-mode isoform expression (`stringtie -e`) | deferred (plan D8 - no current consumer; future default-off stage) |
| 11 | assembly-optimized STAR args (`STAR_ARGS_ASSEMBLY`) | dropped (input BAMs come from the adopted rnaseq route; no local alignment) |

## 6. Local code still necessary (Q6)

`nextflow/pipelines/lncrna-nf/`: workflows/identify.nf wiring the chain;
vendored central modules (stringtie, stringtiemerge, gffcompare, gffread,
cpc2, diamond/blastx, subread/featurecounts); local modules for CPAT and
pfam_scan and the final-set assembly; bin/ helpers (classcode-u parsing,
coding-potential intersection, pfam/diamond hit parsing, final-set build,
count merge). No alignment code - the route consumes the adopted rnaseq
run's published BAMs.

## 7. Validation demonstrating scientific usability (Q7)

Per owner decision (design Section 15 amendment 8): no real-data round.
Evidence at this stage: synthetic fixtures with planted novel transcripts
(some coding - must be filtered out; some noncoding - must be recovered),
mini Pfam HMM + mini DIAMOND databases so the homology stages run at
fixture scale, dry-run baselines, pytest unit suites for every helper, and
the DEFERRED real-container fixture run (sensitivity: planted noncoding
recovered; specificity: planted coding excluded). versions.yml tool rows
follow the ledger rule - only after the validated fixture run.

## 8. Decision log

- **D1 execution model** - local DSL2 `pipelines/lncrna-nf/` (lncpipe
  immature; no adoption). 
- **D2 assembly** - StringTie --conservative + --merge vs gtf_pcg kept.
- **D3 coding potential** - CPC2 AND CPAT 3.0.5 AND length >= 200
  intersection; CNCI retired with this evidence.
- **D4 homology filters** - pfam_scan + DIAMOND kept, containerized,
  default-ON with per-project OFF switches.
- **D5 final set + expression** - final FASTA/GTF + ref_withLncRNA.gtf +
  featureCounts; matrices named `lncrna.*` to stay distinct from the
  Stage-2 gene-level matrices.
- **D6 species neutrality** - params only; rice documented as the example.
- **D7 external databases + gtf_pcg** - Pfam dir + nr.dmnd stay HPC path
  params; the legacy `gtf_PcGs` placeholder is replaced by a documented
  build recipe (docs section in the pipeline README).
- **D8 isoform quantification** - deferred (no consumer; future default-off
  stage).
- **D9 plan-document tracking** - the Stage-4 plan stays local under
  `plans/`.

## 9. Version pin

- Pipeline-local pins (recorded in versions.yml only after the validated
  fixture run): stringtie 3.0.3 (central module), gffcompare/gffread per
  central modules, cpc2 per central module, cpat 3.0.5, pfam_scan 1.6,
  diamond per central module.
- Watch item: nf-core/lncpipe release (re-review if it ever ships stable).
