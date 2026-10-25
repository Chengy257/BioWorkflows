# Legacy Snakemake suite retirement record

Date: 2026-10-26. Owner decision, recorded per the frozen migration mainline
design doc Section 4 (Stage 4 - Legacy retirement) and D10 (Git-native legacy
preservation). This is the "document deliberately dropped capabilities" item
of the retirement checklist; the "confirm use cases are covered" item is the
owner's own research-need call, given explicitly with the retirement
authorization.

## What was retired

All five Snakemake workflows and their shared layer, removed from the active
tree in one retirement commit:

- `rna-seq/` (v0.8.x)
- `chip_cuttag_atac_faire/` (v0.6.x)
- `seclip-seq/` (v0.2.x)
- `srna-seq/` (v0.2.x)
- `bs-seq/` (v0.2.x)
- `shared/` (the Snakemake-only runtime layer: `python/bioworkflows_runtime.py`
  the WorkflowSpec framework, `python/bioworkflows_versions.py`, and
  `lib/launcher.sh` - consumed in place by the five workflows, nothing else)

## Preservation (both copies verified at retirement time)

1. Annotated git tag `legacy-snakemake-suite` at commit `c1ff802` (the last
   commit carrying the full suite). Recover any directory with:
   `git checkout legacy-snakemake-suite -- <dir>`
2. Local archive copy (HPC machine-local, outside the repository):
   `/home/chengyu/BioWorkflows_archive/2026-10-26_legacy_snakemake_suite/`
   containing all six directories (2.1 MB, `cp -a`).

## Replacement status per workflow

| Legacy workflow | Nextflow replacement | Status |
|---|---|---|
| rna-seq (alignment, quantification, DEG, enrichment, lncRNA chain) | `nextflow/adopted/rnaseq` (3.27.0) + `pipelines/rnaseq-ext` + `adopted/differentialabundance` (2.0.0) + `pipelines/lncrna-nf` | COVERED (star_salmon route; featureCounts route + enrichment in rnaseq-ext; DE via DA; lncRNA chain as lncrna-nf) |
| chip_cuttag_atac_faire - ChIP-seq & CUT&Tag branches | `nextflow/adopted/chipseq` + `pipelines/chipseq-ext` (7 default-off extension stages, real-container validated) | COVERED |
| chip_cuttag_atac_faire - ATAC-seq & FAIRE-seq branches | `nextflow/adopted/atacseq` + `pipelines/atacseq-ext` | COVERED (execution validation of the ext stages still deferred) |
| seclip-seq | none | NOT REPLACED - route not started (design Section 5.8: fresh CLIP/eCLIP ecosystem review before any port; the old pipeline is NOT a porting target) |
| srna-seq | none | NOT REPLACED - route not started (local-foundation route planned: absorb nf-core/smrnaseq capabilities) |
| bs-seq | none | NOT REPLACED - route not started (nf-core/methylseq adopt-first) |

## Deliberately dropped capabilities (no replacement anywhere)

Recorded so nothing is rediscovered as a surprise later:

- The unified Snakemake launcher model (`run.sh` preflight/scheduler
  auto-detection, `config.local.yaml` layering, species presets merge,
  per-rule resources.yaml) - the Nextflow entries use native nf-core
  parameters, machine-local executor config, and per-process resources.
- rna-seq legacy options not carried into the rnaseq route: DEGgroupCompare
  (dropped at Stage 2, owner round 2), the legacy numerical DEG/enrichment
  output layout (D9: scientific acceptance, not parity), trim_galore-era
  adapter-trim artifact naming.
- chip SE-mode support for ATAC/FAIRE (dropped at Stage 3), the mixed-assay
  single-sheet seqtype switch (D3 retires the combined-workflow design).
- The Snakemake 7.32.4 cluster profiles (SGE/SLURM/PBS) - the Nextflow
  entries run local or PBS via machine-local config; per-scheduler porting
  happens when a real project needs it.
- Small-RNA cascade filter, novel-miRNA discovery (miRDeep-P2), Bismark
  methylation calling + methylKit DMR, PureCLIP/eCLIP chain: live ONLY in
  the archive/tag until their routes are started.

## Retirement checklist (design doc Stage 4)

- [x] Dropped capabilities documented (this file)
- [x] Final state preserved by annotated tag (`legacy-snakemake-suite`)
- [x] Local archive copy (BioWorkflows_archive)
- [x] Currently required use cases: owner confirmation given with the
      retirement authorization (2026-10-26); small-RNA/BS-seq/CLIP assays
      remain archive-only until their routes start
- [x] Root README rewritten (assay-route table points at the Nextflow
      entries; retired section added)
