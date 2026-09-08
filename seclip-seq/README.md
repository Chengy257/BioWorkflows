# seclip-seq

A **Snakemake** workflow for single-end enhanced CLIP (eCLIP-style) protein-RNA binding maps. One command takes raw SE FASTQ files with a 10-base random UMI at the start of each read through: UMI extraction, two-pass 3' adapter trimming, read sorting, an optional sncRNA/repeats pre-filter, unique genome alignment, UMI-based deduplication, and peak calling with PureCLIP (CLIPper when configured), ending in a combined MultiQC report. Two optional v0.2 stages (both default off) add cross-sample reproducible peaks per condition (bedtools multiinter consensus) and GTF-based peak annotation (nearest gene / distance / biotype per peak set); the consensus stage can additionally peak-call `role: input` control samples and flag (optionally filter) the consensus against their per-condition background union.

> **Status (v0.2.0)**: first-class subproject aligned with the `rna-seq` / `chip_cuttag_atac_faire` engineering model — unified `run.sh` launcher (four scheduler profiles + auto detection + preflight + resource overrides + unlock), one main software environment resolved via `config/software.yaml` (never auto-created by Snakemake), per-rule cluster resources in `config/resources.yaml`, species presets, three-layer config stacking, parse-time config/sample-table validation, software-version provenance, and a synthetic-data dry-run regression test (plus a pytest suite for the annotation scripts). The default regression DAG is **23 jobs** (both optional stages off); the `--consensus` scenario adds 5 (28 total); the `--input-control` scenario (2 ip + 2 input samples, background flagging/filtering on) dry-runs 45 jobs; for end-to-end validation see [Development and testing](#development-and-testing).

## Workflow overview

| Stage | Rule(s) | What happens |
|---|---|---|
| `0.index` | `star_index_genome` | STAR genome index from the genome FASTA + GTF (`--sjdbOverhang`, `--genomeSAindexNbases`) |
| | `star_index_repeats` | STAR index over the sncRNA/repeats FASTA (small reference, relaxed RAM setting) |
| `2.cleandata` | `umi_extract` | `umi_tools extract` with the configured UMI pattern (default 10 random N bases) |
| | `cutadapt_trim` | two-pass 3' adapter trimming against the 20 shifted adapter variants (relaxed `-O 1` sweep, then strict `-O 5`) |
| | `fastq_sort` | `seqkit sort -n` + bgzip, deterministic order before alignment |
| | `fastqc` | FastQC on the sorted trimmed reads |
| `3.align` | `star_filter_repeats` | when `filter_repeats: true`: align to the repeats index first, keep only unmapped reads |
| | `star_align` | unique genome alignment (end-to-end, `--outFilterMultimapNmax 1` default) |
| `4.rmdup` | `umi_dedup` | `umi_tools dedup --method unique` + samtools sort/index |
| | `read_count` | `samtools view -c -F 4` mapped-read count of the deduplicated BAM |
| `5.callpeak` | `callpeak_pureclip` | PureCLIP crosslink-site peaks on the deduplicated BAM |
| | `callpeak_clipper` | CLIPper peak clusters, only when a CLIPper executable is configured |
| `6.reproducible_peaks` (optional) | `consensus_peaks` | when `reproducible_peaks.enabled`: bedtools multiinter consensus of one condition's ip-role PureCLIP beds; sites present in >= `min_replicates` beds are kept, column 4 of the consensus BED reports the support |
| | `input_background` | with `reproducible_peaks.input_control`: union of the condition's input-control PureCLIP beds (per-condition background) |
| | `flag_input_background` / `filter_input_background` | with `input_control`: append the binary `in_input_background` column to the consensus BED; with `filter_by_input` additionally write the consensus minus background sites |
| `6.annotation` (optional) | `gtf_gene_regions` | when `annotate_peaks.enabled`: GTF -> gene/exon BEDs + gene attribute table (stdlib parser) |
| | `annotate_sample_peaks` / `annotate_consensus_peaks` | bedtools intersect (exon/gene overlap) + closest (nearest gene + distance) merged into one TSV per peak set |
| `5.QC` | `multiqc` | combined QC report over all per-sample modules |
| | `software_versions` | record of the tool versions actually resolved for the run |

A rendered DAG of the regression-test project is in [docs/dag_test.svg](docs/dag_test.svg).

## Environment setup

Three ways to get an environment (pick one; details in [docs/user-guide.md](docs/user-guide.md) §1):

1. **Fresh server**: `mamba env create -f workflow/environment.yaml` (all-in-one environment `seclip-seq`, pinning snakemake-minimal 7.32.4 / STAR 2.7.10b / samtools 1.17 / bedtools 2.31.0 / cutadapt 4.6 / umi-tools 1.1.5 (pip, python 3.9) / seqkit 2.13.0 / fastqc 0.11.9 / multiqc 1.21 / pureclip 1.3.1);
2. **Reuse an existing conda environment**: in `config/software.yaml` set `environment.type: conda` + `conda_prefix` (recommended) or `conda_name`; `run.sh` injects the prefix's `bin` into PATH automatically, no activate needed;
3. **System mode**: with `environment.type: system`, all tools come from PATH.

Environments are created explicitly by the user; Snakemake never deploys them automatically. Before launching, preflight with `bash run.sh -P <workdir> --check-software` (required executables + configured CLIPper + databases). CLIPper is deliberately not packaged: point `paths.clipper` in `software.yaml` at an existing absolute install path, or leave it empty and the workflow auto-skips CLIPper (PureCLIP only) with a warning.

Snakemake version matrix: **7.32.4** is the reference version (pinned in `workflow/environment.yaml`; the four cluster profiles target the 7.x classic `--cluster` interface). Snakemake 8.x moved cluster submission to the executor plugin system — parsing and dry-runs work, but test before real cluster runs.

## Quick start

The full walkthrough is in [docs/user-guide.md](docs/user-guide.md) §2; a ready-to-edit project lives in [example/](example/README.md). Summary:

```bash
# 1) Working directory and data (SE fastq: {sample}.fastq.gz or {sample}_R1.fq.gz)
mkdir -p ~/work/fbl/1.rawdata
cp FC_rep1.fastq.gz FC_rep2.fastq.gz ~/work/fbl/1.rawdata/

# 2) Sample table (sample_id column, optional condition/role columns) + project config
cp example/samples.csv ~/work/fbl/samples.csv
cp example/config.yaml  ~/work/fbl/config.yaml       # fill in the reference paths

# 3) Preflight -> dry-run -> run (from the repository root / this subproject)
bash run.sh -P ~/work/fbl --check-software
bash run.sh -P ~/work/fbl -n
bash run.sh -P ~/work/fbl -j 10
# positional form:  bash run.sh ~/work/fbl ~/work/fbl/config.yaml 10
# cluster example:  bash run.sh -P ~/work/fbl --profile pbs --queue workq --memory 16G --runtime 600
```

Key points:

- Raw inputs live in `1.rawdata/` at the working-directory root; every derived artifact goes under `results/` (rename via `results_dir` in config) in numbered stage dirs;
- Config chain: repository `config/config.yaml` defaults → project config (`-c` or positional; auto-detected `<workdir>/config.yaml`) → `config.local.yaml` in the working directory (auto-layered; later files win);
- `species: "hsa"` selects the GRCh38.p14 / GENCODE v44 reference preset from `config/species.yaml`; explicit `genome`/`gtf`/`repeats_fa` keys in the project config win over the preset;
- Cluster jobs are submitted with the per-rule resources declared in `config/resources.yaml` (`threads/mem_mb/runtime_min`); copy that file into the working directory as `resources.yaml` to override (run.sh auto-detects it), or override globally via `--memory`/`--runtime`. Full options in `run.sh --help`.

## Directory structure

```
seclip-seq/
├── run.sh                    # unified launcher CLI (four profiles / auto detection / preflight / --unlock; see run.sh --help)
├── workflow/
│   ├── Snakefile             # single entry point (species presets + config layering + target aggregation)
│   ├── environment.yaml      # all-in-one conda environment template (pinned; created explicitly by the user)
│   ├── rules/                # common / index / upstream / align / callpeak / consensus / meta
│   ├── scripts/              # runtime_config.py (software.yaml resolver), collect_versions.py,
│   │                         # gtf_to_gene_regions.py + annotate_peaks.py (annotation stage)
│   ├── profile/              # default / pbs / sge / slurm profiles + README (cluster commands and pinned params)
│   └── multiqc_config.yaml
├── config/
│   ├── config.yaml           # repository defaults (hsa preset paths; outputs rooted at results/)
│   ├── config.template.yaml  # annotated project template (copy into the working directory)
│   ├── species.yaml          # species presets (hsa: genome / gtf / repeats_fa)
│   ├── resources.yaml        # per-rule scheduler resources (threads/mem_mb/runtime_min)
│   ├── software.yaml         # unified software runtime (conda_prefix / system; paths.clipper external path)
│   └── samples.csv           # sample table template (single sample_id column)
├── tests/                    # run_test.sh (dry-run + --real-run regression) / lint.sh / make_testdata.py
├── example/                  # example project (legacy FBL-CLIP sample table + filled config + start guide)
├── docs/                     # user guide + TODO backlog + test DAG
├── Makefile                  # make check / lint / test
├── CHANGELOG.md
└── LICENSE                   # Apache-2.0 (CI lives at the repository root: .github/workflows/ci.yml)
```

A project working directory looks like this:

```
workdir/
├── config.local.yaml        # optional config overlay (auto-layered, usually kept out of version control)
├── samples.csv              # sample table (filename free; pointed to by SampleListFile)
├── snakemake.logs.txt       # launcher/Snakemake main log (change with --log FILE)
├── 1.rawdata/               # raw SE fastq inputs (the only data directory outside results/)
└── results/                 # ALL derived artifacts (rename via results_dir)
    ├── 0.index/             # STAR genome / repeats indices
    ├── 2.cleandata/         # UMI-extracted, trimmed, sorted fastq + FastQC
    ├── 3.align/             # repeats-unmapped reads + unique genome alignments
    ├── 4.rmdup/             # deduplicated BAMs + read counts
    ├── 5.callpeak/          # PureCLIP / CLIPper peak BEDs
    ├── 6.reproducible_peaks/  # optional: per-condition consensus BEDs (support in column 4; input flag with input_control)
    ├── 6.annotation/        # optional: {set}.annotation.tsv + _ref/ gene/exon BEDs
    ├── 5.QC/                # multiqc_report.html + software_versions.yaml
    └── logs/                # per-rule logs
```

## Results path quick reference

All derived artifacts live under `results/` in the working directory (rename via `results_dir` in config); raw inputs `1.rawdata/` stay at the working-directory root. `{sample}` is a `sample_id` from the sample table.

| Result | Path |
|---|---|
| STAR genome index | `results/0.index/genome_STARindex/` (checkpoint: `.../SA`) |
| STAR repeats/sncRNA index (`filter_repeats: true`) | `results/0.index/repeats_STARindex/` (checkpoint: `.../SA`) |
| UMI-extracted reads | `results/2.cleandata/{sample}_umi.fq.gz` |
| umi_tools extract metrics | `results/2.cleandata/logs/{sample}_umi_extract.metrics` |
| Trimmed reads (after both cutadapt passes) | `results/2.cleandata/{sample}_clean.fqTrTr.fq.gz` |
| cutadapt report (second pass) | `results/2.cleandata/logs/{sample}_cutadapt.metrics` |
| Sorted trimmed reads (alignment input) | `results/2.cleandata/{sample}_clean.fqTrTr.sorted.fq.gz` |
| FastQC report | `results/2.cleandata/fastqc/{sample}_fastqc.html` (+ `.zip`) |
| Repeats-unmapped reads (`filter_repeats: true`) | `results/3.align/repeats/{sample}_Unmapped.out.mate1` |
| Genome alignment BAM | `results/3.align/genome/{sample}_Aligned.out.bam` |
| STAR alignment statistics | `results/3.align/genome/{sample}_Log.final.out` |
| Deduplicated BAM (+ index) | `results/4.rmdup/{sample}.rmDupSo.bam` (+ `.bam.bai`) |
| umi_tools dedup statistics | `results/4.rmdup/{sample}_stats/{sample}_edit_distance.tsv` |
| Mapped-read count | `results/4.rmdup/{sample}_readnum.txt` |
| PureCLIP peaks | `results/5.callpeak/{sample}.pureclip.bed` — 7 columns: BED6 (chromosome, start (0-based), end, site name, crosslink-site score, strand) plus a trailing score-attributes field (`[score_CL=...;score_E=...]`) |
| CLIPper peaks (when configured) | `results/5.callpeak/{sample}.clipper.peakClusters.bed` |
| Per-condition consensus peaks (`reproducible_peaks.enabled`) | `results/6.reproducible_peaks/{condition}.consensus.bed` (with `input_control`, conditions with inputs carry the `in_input_background` flag in column 5) |
| Input background (`reproducible_peaks.input_control`) | `results/6.reproducible_peaks/{condition}.input_background.bed` |
| Filtered consensus (`reproducible_peaks.filter_by_input`) | `results/6.reproducible_peaks/{condition}.consensus.filtered.bed` |
| Peak-set annotation (`annotate_peaks.enabled`) | `results/6.annotation/{set}.annotation.tsv` (`{set}` = each sample id and each `{condition}.consensus`) |
| GTF-derived annotation reference (`annotate_peaks.enabled`) | `results/6.annotation/_ref/genes.bed` (+ `exons.bed`, `genes.tsv`) |
| Combined QC report | `results/5.QC/multiqc/multiqc_report.html` |
| Software version record | `results/5.QC/software_versions.yaml` |
| Per-rule logs | `results/logs/` |

STAR also writes sibling files next to the two declared alignment outputs (`{sample}_Log.out`, `{sample}_SJ.out.tab`, ...); they are informational only.

Optional v0.2 outputs:

- `results/6.reproducible_peaks/{condition}.consensus.bed` — one file per condition (the `condition` of the ip-role samples in the sample table): bedtools multiinter consensus of that condition's per-sample PureCLIP beds, filtered to sites present in >= `min_replicates` beds; BED4 with the support count (number of replicates carrying the site) in column 4. With `reproducible_peaks.input_control: true`, conditions that declare input samples get the binary `in_input_background` flag appended as column 5 (1 = the site also occurs in the condition's input background);
- `results/6.reproducible_peaks/{condition}.input_background.bed` — with `input_control`: union of the condition's input-control PureCLIP beds (BED4, column 4 = number of inputs covering the feature);
- `results/6.reproducible_peaks/{condition}.consensus.filtered.bed` — with `filter_by_input: true` (needs `input_control`): the consensus minus the background-flagged sites (BED4 again);
- `results/6.annotation/{set}.annotation.tsv` — one file per peak set (`{set}` is each peak-called sample id plus each `{condition}.consensus`), with the header `chrom, start, end, score, nearest_gene, nearest_gene_id, distance, feature_class, gene_biotype` (`score` = PureCLIP crosslink-site score, or the consensus support for consensus sets; `feature_class` = `exon` when the peak overlaps an exon, `gene` when it falls in a gene body, `intergenic` otherwise; `distance` = signed bedtools closest distance to the nearest gene, 0 = overlapping, `NA` when the contig has no gene); with `filter_by_input` the consensus annotation consumes the filtered BED, otherwise the flagged BED (its flag column is dropped per the column contract);
- `results/6.annotation/_ref/` — the GTF-derived machinery behind the table above (`genes.bed`, `exons.bed`, `genes.tsv`).

Both stages are off by default; see [docs/user-guide.md](docs/user-guide.md) §5.5 for the config switches and the sample-table columns that drive them.

## QC notes

The MultiQC report aggregates, per sample:

- **umi_tools extract** (`*_umi_extract.metrics`): how many reads carried the expected UMI pattern and entered the pipeline;
- **cutadapt** (`*_cutadapt.metrics`, the strict second pass): adapter-contaminated read fraction and length distribution after trimming;
- **FastQC** (`*_fastqc.zip`): per-base quality, GC content, and other module checks on the trimmed reads;
- **STAR** (`*_Log.final.out`): input reads, uniquely mapped reads and percentage, splices, and mismatch rates for the end-to-end genome alignment;
- **umi_tools dedup** (`*_edit_distance.tsv`): deduplication statistics of the UMI-collapse step.

There is no Bismark or RNA-seq-style strandedness module here; pipeline health is judged from the UMI extraction rate → trimmed yield → unique mapping rate → deduplicated mapped count (`{sample}_readnum.txt`) chain, plus the peak BED sizes in `results/5.callpeak/`. Peak files are per sample (single-sample calling); with `reproducible_peaks.enabled` the per-condition consensus (`results/6.reproducible_peaks/`) and its support counts add the cross-replicate view, and `annotate_peaks.enabled` ships every peak set with a nearest-gene/biotype table (see [docs/TODO.md](docs/TODO.md) for the remaining backlog).

## Differences from the legacy pipeline

The workflow is a refactor of the legacy `seCLIP.smk` Snakemake draft. Four deliberate deviations:

1. **`seCLIP_nofilter.smk` variant replaced by a config switch.** The legacy repo carried a second Snakefile that skipped the sncRNA/repeats filter; here the same behavior is `filter_repeats: false` in config (the repeats index and filter rule drop out of the DAG).
2. **Redundant name-sort step dropped.** The source sorted reads unsorted → by name → by coordinate before deduplication; the intermediate name sort was immediately overwritten by the coordinate sort, so the workflow goes straight from the trimmed reads to the coordinate-sorted dedup input.
3. **FastQC added on trimmed reads.** The legacy pipeline had no read-level QC; the `fastqc` rule now feeds per-base quality plots into MultiQC.
4. **CLIPper de-hardcoded and auto-skipped.** The legacy Snakefile hardcoded the CLIPper executable path and species; the executable is now resolved from `software.yaml` (`paths.clipper`, exported as `SECLIP_TOOL_CLIPPER`) with `--species` from `callpeak.clipper_species`, and when no executable is configured CLIPper is skipped with a warning while PureCLIP still runs.

## Development and testing

```bash
make check    # bash -n syntax checks (no snakemake needed)
make lint     # static check suite tests/lint.sh (missing optional tools are skipped)
make unit     # pytest suite for the annotation scripts (tests/test_gtf_regions.py)
make test     # check + lint + unit
```

Regression tests (need snakemake):

```bash
bash tests/run_test.sh               # synthetic-data dry-run: generate -> assemble working directory -> validate DAG integrity
bash tests/run_test.sh --consensus   # dry-run with the optional v0.2 stages enabled (condition/role table, both stages on)
bash tests/run_test.sh --input-control  # dry-run with the input-control consensus scenario (2 ip + 2 input samples, flagging/filtering on)
bash tests/run_test.sh --real-run    # end-to-end run + output assertions (server validation; needs the full analysis environment)
```

Test data is generated by `tests/make_testdata.py` with a fixed seed (2 x 20 kb chromosomes, 3 snRNA-like repeats, 2 samples of 50 bp SE UMI reads; `--with-inputs` adds the input controls FC_in1/FC_in2) and is never committed. The default dry-run DAG is **23 jobs** (2 samples, `filter_repeats: true`, PureCLIP on, CLIPper absent/auto-skipped, both optional stages off) — compare against this count after refactors (recorded in [CHANGELOG.md](CHANGELOG.md)); the `--consensus` scenario dry-runs 28 jobs and asserts the consensus/annotation rules join the DAG; the `--input-control` scenario dry-runs 45 jobs and asserts the background/flagging/filtering rules join the DAG. A lightweight CI job (lint + `--reads 2000` dry-run regression) runs at the repository root (`.github/workflows/ci.yml`, added 2026-09-07); locally, use `make test` + `bash tests/run_test.sh --reads 2000`. Real-run validation uses the script's `--reads 50000` default — PureCLIP's parameter learning needs tens of thousands of reads.

## Known limitations and TODOs

Mirrors [docs/TODO.md](docs/TODO.md):

1. ~~**Cross-sample reproducible peaks**~~ — resolved in v0.2.0 (optional `reproducible_peaks` stage); IDR-style ranking remains a possible future refinement.
2. ~~**Peak annotation**~~ — resolved in v0.2.0 (optional `annotate_peaks` stage); repeat-feature annotation beyond the GTF is future work.
3. ~~**IP vs input control**~~ — resolved in v0.2.0 (`reproducible_peaks.input_control` / `filter_by_input`: condition-scoped input background, flagged/filtered consensus); IDR-style contrast modes remain possible future work.

## License

[Apache-2.0](LICENSE) © 2026 ChengYu

## Versions

v0.2.0 (2026-09-08): optional cross-sample reproducible peaks (per-condition bedtools multiinter consensus), GTF-based peak annotation stages, and the IP-vs-input extension (`input_control` / `filter_by_input`: condition-scoped input background, flagged/filtered consensus), driven by the optional condition/role sample-table columns; bedtools 2.31.0 pinned; pytest unit suite and the `--consensus` / `--input-control` regression scenarios (see [CHANGELOG.md](CHANGELOG.md)).

v0.1.0 (2026-09-05): initial release — config layer, STAR index / UMI / trim / align / dedup / peak-calling rules, unified launcher, regression tests, and docs (see [CHANGELOG.md](CHANGELOG.md)).
