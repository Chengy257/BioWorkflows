# pipelines/lncrna-nf - BioWorkflows lncRNA identification workflow

Local DSL2 pipeline (the design-doc local-foundation route): consumes the
adopted nf-core/rnaseq route's published BAMs and runs the ratified
identification chain. Scientific record:
`nextflow/docs/lncrna-scientific-decision-record.md`. This is NOT an
upstream adoption - nf-core/lncpipe remains unreleased (watch item).

## Input contract (the Stage-2 rnaseq handoff)

| Input | Param | Notes |
|---|---|---|
| Published BAMs | `bam_dir` | `{rnaseq_outdir}/{aligner}`; default suffix `.markdup.sorted.bam` (publish unconditionally on rnaseq 3.27.0); `bam_suffix = '.sorted.bam'` consumes the pre-dedup STAR BAMs when the adopted run passed `--save_align_intermeds` |
| Sample sheet | `sample_sheet` | CSV with a `sample_id` column (the Stage-2 samplesheet works as-is) |
| Genome FASTA | `genome` | sibling `.fai` must exist |
| Annotation GTF | `gtf` | |
| Protein-coding GTF | `gtf_pcg` | see the build recipe below |
| Strandedness | `strand` | `rf` / `fr` / empty - from the Stage-2 strandedness inference |

## Chain (decision-record dispositions)

StringTie per-sample `--conservative -j 5` assembly -> `stringtie --merge`
vs `gtf_pcg` (`-F 0.1 -T 0.1`) -> gffcompare vs `gtf_pcg` ->
**classcode-u novel extraction** -> candidate FASTA -> coding potential
**CPC2 AND CPAT 3.0.5 AND length >= 200** (CNCI retired: python2-only,
unmaintained) -> noncoding FASTA -> **pfam_scan** (`-translate`, E < 1e-5)
and **DIAMOND blastx vs nr.dmnd** homology filters (default ON, per-project
OFF switches) -> **final lncRNA set** + `ref_withLncRNA.gtf` ->
featureCounts -> `lncrna.count/tpm/fpkm.matrix.tsv` (named distinctly from
the Stage-2 gene-level matrices). Central nf-core modules are vendored
verbatim; CPAT, pfam_scan and the set-building steps are local modules.

## The gtf_pcg build recipe (rice example)

The legacy `lncrna.gtf_PcGs` placeholder is a user-built artifact: the
protein-coding CDS entries of the MSU (IRGSP-1.0) and RAP-DB annotations
merged into one GTF. Build it once per reference release:

```bash
# 1. Download both GFF3/GTF releases for the same genome build.
# 2. Keep only protein-coding CDS/exon features (gene biotype filter).
# 3. Concatenate and coordinate-sort:
cat msu_pcg.gtf rap_pcg.gtf | sort -k1,1 -k4,4n -k5,5n -u > osa_pcg.gtf
```

## CPAT species models

CPAT requires a species logit model + hexamer table (`cpat.model`,
`cpat.hexamer`). The human files ship with CPAT; for rice, build the model
with `cpat.py make_logitModel` against a curated rice coding/noncoding
training set, or use a published rice model. The pipeline never bundles a
species model (species neutrality, design D4).

## Run

```bash
nextflow run <this pipeline directory> \
    -profile apptainer \
    -params-file my-project.params.yaml \
    -ansi-log false
```

Dry-run baselines (fixtures): `bash tests/make_fixtures.py` then
`nextflow run . -profile test` / `-profile test_full` (see the validation
notes for the recorded counts). `tests/unit/` pins every helper's
semantics via pytest; run everything through the repo-level
`nextflow/tests/static-checks.sh` (no hosted CI).

## Deferred validation

Real-container fixture execution (sensitivity: planted noncoding
transcripts recovered; specificity: planted coding transcripts excluded)
is deferred with the owner's real-data skip decision of 2026-10-12 (design
amendment 8). versions.yml tool rows follow the ledger rule - after the
validated fixture run.
