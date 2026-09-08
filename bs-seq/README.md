# bs-seq

A **Snakemake** workflow for bisulfite sequencing (BS-seq) methylation calling. One command takes raw FASTQ files (paired-end or single-end, auto-detected per sample) through: optional 3' adapter/quality trimming (Trim Galore), bisulfite genome preparation (`bismark_genome_preparation`), Bismark alignment (bowtie2 backend), PCR-duplicate removal (`deduplicate_bismark`), per-sample methylation extraction (`bismark_methylation_extractor`) with an optional merged CpG table (`coverage2cytosine --merge_CpG`) and a per-sample HTML report, ending in a run-level `bismark2summary` and a combined MultiQC report — plus an optional methylKit differential-methylation stage (DMCs/DMRs over the merged CpG tables, `dmr.enabled`, default off) into `6.DMR/`.

> **Status (v0.2.0)**: first-class subproject aligned with the `rna-seq` / `srna-seq` / `chip_cuttag_atac_faire` engineering model — unified `run.sh` launcher (four scheduler profiles + auto detection + preflight + resource overrides + unlock), one main software environment resolved via `config/software.yaml` (never auto-created by Snakemake), per-rule cluster resources in `config/resources.yaml`, species presets (rice `osa` / human `hsa`), three-layer config stacking, parse-time config/sample-table validation, PE/SE auto-detection, software-version provenance, and a synthetic-data dry-run regression test. The default regression baseline DAG is **20 jobs** (`--dmr` scenario: 35); for end-to-end validation see [Development and testing](#development-and-testing).

## Workflow overview

| Stage | Rule(s) | What happens |
|---|---|---|
| `0.index` | `bismark_genome_prep` | copies the configured genome FASTA into `0.index/bismark_genome/` and runs `bismark_genome_preparation` (C→T and G→A bisulfite-converted bowtie2 indexes) |
| | `bam2nuc_genome` | genome-wide nucleotide composition totals (`bam2nuc --genomic_composition_only`) |
| `2.cleandata` | `trim_pe` / `trim_se` | Trim Galore 3' quality + adapter trimming, `--gzip`, with `--fastqc` on the trimmed reads; the rule (PE or SE) is selected per sample by the detected library layout; skipped entirely when `trim.enabled: false` (raw reads are aligned directly) |
| `3.align` | `bismark_align` | Bismark (bowtie2 backend) alignment of the trim-aware reads against the bisulfite genome |
| `4.dedup` | `deduplicate` | `deduplicate_bismark` removes PCR duplicates from the alignment BAM |
| `5.methylation/{sample}` | `bam2nuc_sample` | per-sample nucleotide composition on the deduplicated BAM |
| | `methylation_extractor` | `bismark_methylation_extractor` on the deduplicated BAM: gzipped cytosine report + bedGraph (+ splitting report and M-bias with a genome) |
| | `coverage2cytosine` | `coverage2cytosine --merge_CpG`: merged CpG-context table per sample (only requested when `methylation_extractor.merge_cpg: true`) |
| | `bismark2report` | per-sample HTML report (alignment + dedup + splitting + M-bias + nucleotide stats) |
| `5.QC` | `bismark2summary` | run-level Bismark summary HTML across all samples |
| | `multiqc` | combined QC report over trimming, FastQC, and the Bismark reports |
| | `software_versions` | record of the tool versions actually resolved for the run |
| `6.DMR` | `dmr_methylkit` | optional (dmr.enabled, default off) single methylKit job over every sample's merged CpG table: per-contrast DMC tables (all/hyper/hypo), tiled DMR table, and a summary; contrasts are treat-vs-`control_group` from the sample table's optional `group` column (>= 2 replicates per group) |

## Bismark output naming

Bismark tools derive most output file names from the **input** file name, so the bs-seq result paths follow a derived-name chain rather than a flat per-stage scheme. The workflow declares these contract paths (kept in sync between the rule modules and `TARGETS` in `workflow/rules/common.smk`):

| Step | Input (basename) | Declared contract output (verified against Bismark 0.24.0, 2026-09-05) |
|---|---|---|
| `bismark_align` (`--basename {sample}`) | trim-aware FASTQs | `results/3.align/{sample}.bam` + `results/3.align/{sample}_report.txt`; the tool writes `{sample}_pe.bam`/`_se.bam` and `{sample}_PE_report.txt`/`_SE_report.txt`, which are renamed/copied — the native report copy and a `{sample}_pe.bam` symlink are kept beside the BAM for `bismark2summary` discovery |
| `deduplicate_bismark` | `3.align/{sample}.bam` | `results/4.dedup/{sample}.deduplicated.bam` + `results/4.dedup/{sample}.deduplication_report.txt` |
| `bismark_methylation_extractor` | `4.dedup/{sample}.deduplicated.bam` | `results/5.methylation/{sample}/{sample}.deduplicated.bismark.cov.gz`, `.../{sample}.deduplicated.bedGraph.gz`, `.../{sample}.deduplicated_splitting_report.txt` (underscore), `.../{sample}.deduplicated.M-bias.txt` |
| `coverage2cytosine --merge_CpG` | `{sample}.deduplicated.bismark.cov.gz` | `results/5.methylation/{sample}/{sample}.CpG_merged.CpG_report.merged_CpG_evidence.cov.gz` (the tool writes `...merged_CpG_evidence.cov`; the rule gzips it; the large `{sample}.CpG_merged.CpG_report.txt` cytosine report stays beside it as a side effect) |
| `bam2nuc` (genome) | prepared genome folder | `results/0.index/bismark_genome/genomic_nucleotide_frequencies.txt` |
| `bam2nuc` (sample) | `4.dedup/{sample}.deduplicated.bam` | `results/5.methylation/{sample}/{sample}.deduplicated.nucleotide_stats.txt` |
| `bismark2report --output {sample}.html` | the per-sample Bismark reports | `results/5.methylation/{sample}/{sample}.html` (the `--output` name is used verbatim) |
| `bismark2summary` | alignment BAMs (via `_pe`/`_se` symlinks) | `results/5.QC/bismark2summary.html` (`-o bismark2summary`; the rule `cd`s into `5.QC/` with absolute paths) |

These names were **reconciled against the pinned Bismark 0.24.0 in the 2026-09-05/06 WSL real-run validation** (13 reconciliation rounds) and confirmed by the full 20/20 real-run on 2026-09-07 (which also verified the `bam2nuc_sample` genome-folder path fix and the `bismark2report` output name). Remaining caveats — per-sample `--parallel` alignment deferred (Bismark 0.24 rejects `--basename` + `--multicore` together; parallelism comes from Snakemake scheduling samples, like the legacy ParaFly model), and `bismark2summary` skipping dedup/splitting stats unless reports sit beside the BAM — are tracked in [docs/TODO.md](docs/TODO.md) §1.

## Environment setup

Three ways to get an environment (pick one; details in [docs/user-guide.md](docs/user-guide.md) §1):

1. **Fresh server**: `mamba env create -f workflow/environment.yaml` (all-in-one environment `bs-seq`, pinning snakemake-minimal 7.32.4 / bismark 0.24.0 / bowtie2 2.5.2 / samtools 1.17 / trim-galore 0.6.10 / fastqc 0.11.9 / multiqc 1.21, plus the R stack for the optional DMR stage: r-base 4.3 / bioconductor-methylkit / r-getopt);
2. **Reuse an existing conda environment**: in `config/software.yaml` set `environment.type: conda` + `conda_prefix` (recommended) or `conda_name`; `run.sh` injects the prefix's `bin` into PATH automatically, no activate needed;
3. **System mode**: with `environment.type: system`, all tools come from PATH.

Environments are created explicitly by the user; Snakemake never deploys them automatically. Before launching, preflight with `bash run.sh -P <workdir> --check-software` (required executables: bismark / bismark_genome_preparation / deduplicate_bismark / bismark_methylation_extractor / coverage2cytosine / bam2nuc / bismark2report / bismark2summary / bowtie2 / samtools / trim_galore / fastqc / multiqc / python3; Rscript is resolved through the `r:` section for the optional DMR stage but is never demanded while `dmr` is disabled).

Snakemake version matrix: **7.32.4** is the reference version (pinned in `workflow/environment.yaml`; the four cluster profiles target the 7.x classic `--cluster` interface). Snakemake 8.x moved cluster submission to the executor plugin system — parsing and dry-runs work, but test before real cluster runs.

## Quick start

The full walkthrough is in [docs/user-guide.md](docs/user-guide.md) §2; a ready-to-edit project lives in [example/](example/README.md). Summary:

```bash
# 1) Working directory and data (PE: {sample}_1.fastq.gz + {sample}_2.fastq.gz; SE: {sample}.fastq.gz)
mkdir -p ~/work/rice/1.rawdata
cp s1_1.fastq.gz s1_2.fastq.gz s2_1.fastq.gz s2_2.fastq.gz ~/work/rice/1.rawdata/

# 2) Sample table (single sample_id column) + project config
cp example/samples.csv ~/work/rice/samples.csv
cp example/config.yaml  ~/work/rice/config.yaml       # fill in the genome path

# 3) Preflight -> dry-run -> run (from this subproject directory)
bash run.sh -P ~/work/rice --check-software
bash run.sh -P ~/work/rice -n
bash run.sh -P ~/work/rice -j 10
# positional form:  bash run.sh ~/work/rice ~/work/rice/config.yaml 10
# cluster example:  bash run.sh -P ~/work/rice --profile pbs --queue workq --memory 32G --runtime 720
```

Key points:

- Raw inputs live in `1.rawdata/` at the working-directory root; every derived artifact goes under `results/` (rename via `results_dir` in config) in numbered stage dirs (`0.index` / `2.cleandata` / `3.align` / `4.dedup` / `5.methylation` / `5.QC` / optional `6.DMR`);
- Config chain: repository `config/config.yaml` defaults → project config (`-c` or positional; auto-detected `<workdir>/config.yaml`) → `config.local.yaml` in the working directory (auto-layered; later files win);
- `species: "osa"` or `"hsa"` fills an unset `genome` path from the `config/species.yaml` preset (explicit non-empty keys win; `species: "none"` disables the fallback entirely);
- Library layout is auto-detected per sample from the `1.rawdata/` file names (PE first: `{sample}_1.fastq.gz` + `{sample}_2.fastq.gz`; then SE: `{sample}.fastq.gz`) — PE and SE samples can coexist in one run;
- Cluster jobs are submitted with the per-rule resources declared in `config/resources.yaml` (13 entries: `threads/mem_mb/runtime_min`); copy that file into the working directory as `resources.yaml` to override (run.sh auto-detects it), or override globally via `--memory`/`--runtime`. Full options in `run.sh --help`.

## Directory structure

```
bs-seq/
├── run.sh                    # unified launcher CLI (four profiles / auto detection / preflight / --unlock; see run.sh --help)
├── workflow/
│   ├── Snakefile             # single entry point (species presets + config layering + target aggregation)
│   ├── environment.yaml      # all-in-one conda environment template (pinned; created explicitly by the user)
│   ├── rules/                # common / index / upstream / align / methylation / dmr / meta
│   ├── scripts/              # runtime_config.py, collect_versions.py
│   ├── profile/              # default / pbs / sge / slurm profiles + README (cluster commands and pinned params)
│   └── multiqc_config.yaml
├── config/
│   ├── config.yaml           # repository defaults (osa preset; outputs rooted at results/)
│   ├── config.template.yaml  # annotated project template (copy into the working directory)
│   ├── species.yaml          # species presets (osa: IRGSP-1.0, hsa: GRCh38.p14 genome)
│   ├── resources.yaml        # per-rule scheduler resources (threads/mem_mb/runtime_min; 13 rules)
│   ├── software.yaml         # unified software runtime (conda_prefix / system)
│   └── samples.csv           # sample table template (single sample_id column)
├── tests/                    # run_test.sh (dry-run + --real-run regression) / lint.sh / make_testdata.py
├── example/                  # example project (two-sample table + filled config + start guide)
├── docs/                     # user guide + TODO backlog
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
├── 1.rawdata/               # raw FASTQ inputs (the only data directory outside results/)
└── results/                 # ALL derived artifacts (rename via results_dir)
    ├── 0.index/             # bisulfite genome (bismark_genome/) + genome-wide nucleotide totals
    ├── 2.cleandata/         # trimmed fastq + trimming reports + FastQC
    ├── 3.align/             # Bismark alignment BAMs + reports
    ├── 4.dedup/             # deduplicated BAMs + dedup reports
    ├── 5.methylation/       # per-sample cytosine/bedGraph/M-bias/splitting outputs + reports
    ├── 5.QC/                # bismark2summary + multiqc_report.html + software_versions.yaml
    ├── 6.DMR/               # optional (dmr.enabled): methylKit DMC/DMR tables + summary
    └── logs/                # per-rule logs
```

## Results path quick reference

All derived artifacts live under `results/` in the working directory (rename via `results_dir` in config); raw inputs `1.rawdata/` stay at the working-directory root. `{sample}` is a `sample_id` from the sample table. See [Bismark output naming](#bismark-output-naming) for how the derived names come about.

| Result | Path |
|---|---|
| Bisulfite genome index | `results/0.index/bismark_genome/Bisulfite_Genome/{CT,GA}_conversion/` |
| Genome nucleotide totals | `results/0.index/bismark_genome/genomic_nucleotide_frequencies.txt` |
| Trimmed reads (PE) | `results/2.cleandata/{sample}_1_val_1.fq.gz` (+ `{sample}_2_val_2.fq.gz`) |
| Trimmed reads (SE) | `results/2.cleandata/{sample}_trimmed.fq.gz` |
| Trim Galore reports | `results/2.cleandata/{sample}_trimming_report.txt` (PE: one per mate, keyed by the mate name) |
| FastQC on trimmed reads | `results/2.cleandata/fastqc/{sample}_trimmed_fastqc.zip` (+ `.html`; PE: per mate) |
| Alignment BAM | `results/3.align/{sample}.bam` |
| Bismark alignment report | `results/3.align/{sample}_report.txt` |
| Deduplicated BAM | `results/4.dedup/{sample}.deduplicated.bam` |
| Deduplication report | `results/4.dedup/{sample}.deduplication_report.txt` |
| Nucleotide stats (per sample) | `results/5.methylation/{sample}/{sample}.deduplicated.nucleotide_stats.txt` |
| Cytosine report | `results/5.methylation/{sample}/{sample}.deduplicated.bismark.cov.gz` |
| Methylation bedGraph | `results/5.methylation/{sample}/{sample}.deduplicated.bedGraph.gz` |
| Splitting report / M-bias | `results/5.methylation/{sample}/{sample}.deduplicated_splitting_report.txt` / `{sample}.deduplicated.M-bias.txt` |
| Merged CpG table | `results/5.methylation/{sample}/{sample}.CpG_merged.CpG_report.merged_CpG_evidence.cov.gz` (when `merge_cpg: true`) |
| Per-sample HTML report | `results/5.methylation/{sample}/{sample}.html` |
| DMC/DMR tables (optional) | `results/6.DMR/{treat}_vs_{control}_DMC_{all,hyper,hypo}.tsv` + `_DMR_tiles.tsv` + `DMR_summary.tsv` (when `dmr.enabled: true`) |
| Run-level Bismark summary | `results/5.QC/bismark2summary.html` |
| Combined QC report | `results/5.QC/multiqc/multiqc_report.html` |
| Software version record | `results/5.QC/software_versions.yaml` |
| Per-rule logs | `results/logs/` |

Reading the outputs: the **cytosine report** (`.bismark.cov.gz`) has one line per covered cytosine — chromosome, start, end, methylation percentage, methylated and unmethylated read counts; the **bedGraph** is the same context-level methylation percentage without the counts; the **splitting report** summarises methylation per sequence context (CpG / CHG / CHH) and strand; the **M-bias file** shows per-read-position methylation (the basis for `--mbias` trimming decisions); the **merged CpG table** is the coverage2cytosine `--merge_CpG` product for CpG sites merged across strands. The per-sample and run-level HTML reports aggregate all of the above.

## QC notes

The MultiQC report aggregates, per sample:

- **Trim Galore** (`*_trimming_report.txt`): adapter-contaminated fraction and the length distribution after trimming;
- **FastQC** (`*_trimmed_fastqc.zip`): per-base quality, GC content, and other module checks on the trimmed reads;
- **Bismark** (`3.align/{sample}_report.txt`, `4.dedup/{sample}.deduplication_report.txt`, `5.methylation/{sample}/{sample}.deduplicated_splitting_report.txt`): alignment efficiency, deduplication rate, and per-context methylation percentages.

`bismark2summary.html` is the run-level view of the same numbers (one row per sample). Pipeline health is read as: raw → trimmed yield → alignment efficiency → duplication rate → CpG methylation level sane for the organism (rice ~80-90% CG context methylation; a near-zero CpG level usually means a wrong genome or a non-bisulfite library). Differential methylation across conditions is out of scope for v0.1 (see [docs/TODO.md](docs/TODO.md)).

## Differences from the legacy pipeline

The workflow reproduces the legacy Bismark BS-seq pipeline (Trim Galore → bismark_genome_preparation → bismark → deduplicate_bismark → bismark_methylation_extractor → coverage2cytosine → reports). The tool invocations and parameters are kept faithful; the deliberate deviations:

1. **Trim stage added, default ON.** The workflow adds explicit Trim Galore adapter/quality trimming before alignment (the legacy flow aligned raw reads). `trim.enabled: false` restores the legacy raw-read alignment exactly.
2. **FastQC / MultiQC added.** FastQC runs on the trimmed reads inside the trim jobs (`trim_galore --fastqc`), and a whole-run MultiQC report aggregates trimming, FastQC, and the Bismark reports; the legacy run left per-step console logs only.
3. **CX / cytosine_report behind a flag.** The full CX-context cytosine report (`--CX --cytosine_report`) is off by default (`methylation_extractor.cx_report: false`) because it is large; the standard contexts + bedGraph + CpG merge cover the legacy outputs.
4. **Orphan `buffer_CpG` line dropped.** The legacy script set a `buffer_CpG` value that nothing consumed; it is gone.
5. **`--buffer_size` derived from resources.** Instead of a hard-coded extractor buffer, the workflow derives `--buffer_size` from the rule's scheduler memory: `mem_mb / methylation_extractor.buffer_frac` (in GB; the cluster memory request and the tool buffer can no longer drift apart).
6. **PE/SE auto-detect.** The legacy script was hard-wired for one library layout; the workflow detects the layout of every declared sample from its `1.rawdata/` file names and picks the matching trim/align/extractor flags, so PE and SE samples can share one run.
7. **ParaFly replaced by Snakemake.** The legacy script ran per-sample command lists through ParaFly with hand-managed logs; Snakemake provides per-sample/per-rule job scheduling, cluster profiles (pbs/sge/slurm), restart safety, and per-rule logs.

## Development and testing

```bash
make check    # bash -n syntax checks (no snakemake needed)
make lint     # static check suite tests/lint.sh (missing optional tools are skipped)
make test     # CI-equivalent full check (= check + lint)
```

Regression tests (need snakemake):

```bash
bash tests/run_test.sh               # synthetic-data dry-run: generate -> assemble working directory -> validate DAG integrity
bash tests/run_test.sh --dmr         # differential-methylation scenario: 4 samples (2x2 group design), asserts dmr_methylkit in the DAG
bash tests/run_test.sh --real-run    # end-to-end run + output assertions (server validation; needs the full analysis environment)
```

Test data is generated by `tests/make_testdata.py` with a fixed seed (2 x 20 kb chromosomes, paired-end 100 bp reads simulated post-bisulfite: C→T on read 1, G→A on read 2, 1% errors, 30% adapter-tailed pairs) and is never committed. The dry-run baseline DAG is **20 jobs** (rule all + 2 x trim_pe + bismark_genome_prep + 2 x bismark_align + 2 x deduplicate + bam2nuc_genome + 2 x bam2nuc_sample + 2 x methylation_extractor + 2 x coverage2cytosine + 2 x bismark2report + bismark2summary + multiqc + software_versions, 2 samples) — compare against this count after refactors (recorded in [CHANGELOG.md](CHANGELOG.md)); the `--dmr` scenario (4 samples, `dmr.enabled: true`) dry-runs **35 jobs** including the single `dmr_methylkit` job.

The regression config runs with `species: "none"` and a relative miniature `ref/genome.fa`, exercising the full repository config schema. A handful of Bismark derived output names were reconciled against the pinned Bismark 0.24.0 in WSL `--real-run` validation (see the naming table above and [docs/TODO.md](docs/TODO.md) §1).

A lightweight CI job (lint + `--reads 2000` dry-run regression) runs at the repository root (`.github/workflows/ci.yml`, added 2026-09-07); locally, use `make test` + `bash tests/run_test.sh --reads 2000`.

## Known limitations and TODOs

Mirrors [docs/TODO.md](docs/TODO.md):

1. **Bismark derived output names**: reconciled against Bismark 0.24.0 in WSL real-run validation (naming table above, [docs/TODO.md](docs/TODO.md) §1). Remaining caveats: `bismark2summary` skips dedup/splitting stats unless reports sit beside the BAM under the tool's own names, and per-sample `--parallel` alignment is deferred (Bismark 0.24 rejects `--basename` + `--multicore`).
2. **Differential methylation / DMR**: done in v0.2.0 — methylKit DMC/DMR calling over the merged CpG tables behind the default-off `dmr:` section (see the workflow table above and [docs/user-guide.md](docs/user-guide.md) §4.2).
3. **No unit-test suite**: the regression is a synthetic-data dry-run; pytest coverage for the pure-Python pieces and the `common.smk` validation is deferred.

## License

[Apache-2.0](LICENSE) © 2026 ChengYu

## Versions

v0.2.0 (2026-09-08): optional methylKit differential-methylation stage (`dmr:` section, `6.DMR/` outputs, optional group/batch sample-table columns, R stack in the environment, `--dmr` test scenario) — see [CHANGELOG.md](CHANGELOG.md).

v0.1.0 (2026-09-05): initial release — config layer, bisulfite-genome / trim / Bismark alignment / dedup / methylation-extraction / reporting rules, unified launcher, pinned environment, regression tests, and docs (see [CHANGELOG.md](CHANGELOG.md)).
