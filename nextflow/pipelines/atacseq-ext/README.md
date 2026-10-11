# pipelines/atacseq-ext - local extension stages for the adopted ATAC/FAIRE route

Local DSL2 extension stages for the adopted nf-core/atacseq entry
(`nextflow/adopted/atacseq/`, pin `6a9307709662bad3e7bcc10aa873f7e81fdd1646`).
Serves ATAC-seq AND FAIRE-seq projects (the FAIRE route decision, D4: same
adopted entry, footprint stays off). This is the chipseq-ext model: consume
the adopted run's published results, add only the ratified local capabilities
as default-off stages. The scientific record is
`nextflow/docs/atacseq-scientific-decision-record.md`.

## Consumed published contract (pinned dev HEAD, from `conf/modules.config`)

| Artifact | Published location (under the adopted run's `outdir`) |
|---|---|
| Filtered BAM (+ BAI) | `{aligner}/merged_library/{id}.mLb.clN.sorted.bam(+.bai)` - unconditional |
| flagstat (+ stats/idxstats) | `{aligner}/merged_library/samtools_stats/{id}.mLb.clN.sorted.bam.flagstat` |
| Picard metrics | `{aligner}/merged_library/picard_metrics/{id}.mLb.mkD.sorted.metrics.txt` |
| Consensus peaks | `{aligner}/merged_library/macs3/{peak_type}/consensus/consensus_peaks.mLb.clN.bed` |

`ext_results_dir` points at the adopted run's `outdir`; the `ext_aligner` /
`ext_peak_type` params must match the adopted run (`bowtie2` / `broad_peak`
defaults here = the documented plant route; upstream defaults are `bwa` /
broad). Every default has an explicit override param (`ext_bam_dir`,
`ext_flagstat_dir`, `ext_consensus_bed`, ...).

Re-run the consumed-path contract audit whenever the adopted pin changes
(Stage-1 watch-item rule).

## Sample sheet

CSV with header `sample_id,group,condition,control,batch` (only
`sample_id,group` are required). A group's `control` column names its control
sample (empty when none); one control may serve several groups. Controls are
processed but excluded from footprint group pooling.

## Stages (all default-off)

| Stage | Params gate | What it does |
|---|---|---|
| `organelle_qc` | `organelle_qc.enabled` | Plant ATAC/mito diagnostic: per-sample `samtools idxstats` -> organelle read-fraction table (chloroplast/mito patterns configurable). |
| `qc_gates` | `qc_gates.enabled` | QC gate table joining published flagstats (mapping rate) + picard metrics (dup rate) + the organelle summary; thresholds are config, gates are informational. |
| `spike_in` | `spike_in.enabled` | Combined-reference runs (host+spike concat index): scale factors from spike-contig idxstats, optional rescaled bigWigs (`spike_in.rescale_bigwigs`). |
| `footprint` | `footprint.enabled` | Per-group TOBIAS chain via the central nf-core tobias modules: pool group treats (`samtools merge`) -> `TOBIAS ATACorrect` (genome FASTA + sibling `.fai` + consensus peaks) -> `TOBIAS ScoreBigwig` -> one `TOBIAS BINDetect` over all groups (needs >= 2 groups; JASPAR motifs; xlsx export skipped). Outputs under `footprint/{group}/{atacorrect,scorebigwig}` and `footprint/bindetect`. |

Known trap (carried from the legacy chip v0.6 stage): BINDetect can exit on a
sklearn NaN with degenerate/toy signal - re-validate on real research data.

## Run

```bash
nextflow run <this pipeline directory> \
    -profile apptainer \
    -params-file my-project.ext-params.yaml \
    -ansi-log false
```

with a params file carrying `ext_results_dir`, `ext_sample_sheet`, and the
enabled stages, e.g.:

```yaml
ext_results_dir: "/path/to/adopted-run/results"
ext_sample_sheet: "/path/to/sheet.csv"
organelle_qc:
  enabled: true
qc_gates:
  enabled: true
```

Dry-run baselines (fixtures): `bash tests/make_fixtures.py` then

```bash
nextflow run . -profile test        # all-off: exit 0, 0 tasks
nextflow run . -profile test_full   # all-on: 22 tasks over 9 processes
```

## Testing

- `tests/make_fixtures.py` - fixed-seed synthetic fixtures mirroring the
  pinned publish layout (git-ignored `tests/fixtures/`).
- `tests/unit/` - pytest suite pinning the helper-script semantics
  (organelle/spike contig pattern convention, QC-gates file join); runs via
  the repo-level `nextflow/tests/static-checks.sh` (no hosted CI).
- Real-container fixture execution (Stage-1 model) and a
  real-BAM fixture helper are deferred with the stage's execution-validation
  round (owner decision 2026-10-12, design amendment 8).
