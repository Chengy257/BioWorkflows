# srna-seq

A **Snakemake** workflow for small-RNA sequencing (plant / plant+microbiome). One command takes raw single-end FASTQ files through: 3' adapter trimming (Trim Galore), a sequential bowtie1 cascade filter over ordered sncRNA classes (rRNA → snoRNA → snRNA → tRNA → miRNA → mRNA → optional exogenous classes), a bowtie1 genome alignment of the trimmed reads, per-class count and RPM matrices, ending in a cascade read-fate summary and a combined MultiQC report — plus two optional, default-off stages: differential expression (DESeq2 over the per-class count matrices, `6.DEG/`) and novel miRNA discovery (miRDeep-P2, `6.novel_mirna/`).

> **Status (v0.2.0)**: first-class subproject aligned with the `rna-seq` / `chip_cuttag_atac_faire` engineering model — unified `run.sh` launcher (four scheduler profiles + auto detection + preflight + resource overrides + unlock), one main software environment resolved via `config/software.yaml` (never auto-created by Snakemake), per-rule cluster resources in `config/resources.yaml`, species presets, three-layer config stacking, parse-time config/sample-table validation, software-version provenance, and a synthetic-data dry-run regression test. v0.2.0 adds the optional DESeq2 differential-expression stage (default off; see [Differential expression](#differential-expression)) and the optional miRDeep-P2 novel-miRNA discovery stage (default off; see [Novel miRNA discovery](#novel-mirna-discovery-optional-default-off)). The regression baseline DAG is **27 jobs**; for end-to-end validation see [Development and testing](#development-and-testing).

## Workflow overview

| Stage | Rule(s) | What happens |
|---|---|---|
| `0.index` | `bowtie_index` | bowtie1 index per configured cascade class (`0.index/{class}/`) plus the genome index (`0.index/genome/`) — built automatically from the configured FASTAs |
| `2.cleandata` | `trim` | Trim Galore (SE): 3' quality + adapter trimming, `--gzip`, with `--fastqc` on the trimmed reads |
| `3.align/filter` | `cascade_stage` | one bowtie1 stage per class, chained via `--un`: mapped reads are assigned to the class, unmapped reads pass to the next stage |
| `3.align/genome` | `genome_align` | bowtie1 alignment of the **trimmed reads** (not the cascade remainder) against the reference genome |
| `4.expression` | `count_stage` | per-sample per-class counts: primary alignments per reference feature in each cascade SAM |
| | `raw_count` | raw read count per sample (feeds the cascade summary) |
| | `merge_counts` | per-class count/RPM matrices + the combined `all_classes_counts.tsv` |
| `5.QC` | `cascade_summary` | per-sample read fate across the cascade (wide TSV + MultiQC custom bargraph) |
| | `multiqc` | combined QC report over all per-sample modules + the cascade bargraph |
| | `software_versions` | record of the tool versions actually resolved for the run |
| `6.DEG` (optional) | `deg_deseq2` | DESeq2 differential expression per analyzed class (default off; `deg.enabled: true` + a `group` column in the sample table) |
| `6.novel_mirna` (optional) | `novel_mirna` | miRDeep-P2 novel miRNA prediction per sample from the trimmed reads against the genome (default off; `novel_mirna.enabled: true` + a mature-miRNA fasta + a configured genome) |

## Cascade design

The cascade is the heart of the workflow; its semantics are fixed by the ordered `cascade:` config list:

1. **Ordered, chained stages.** Classes run strictly top to bottom. Stage *k* aligns stage *k-1*'s unmapped reads (`3.align/filter/{class}/{sample}_unmapped.fq`) against the class index; the first stage starts from the trimmed reads. A read is assigned to the **first** class whose index it maps to, and leaves the flow; everything else continues (`bowtie ... --un` chaining).
2. **Data-driven class list.** Any number of classes with arbitrary (legal) names — not just the shipped rRNA/snoRNA/snRNA/tRNA/miRNA/mRNA/rhizo set. Names become directory names and bowtie arguments.
3. **Empty fasta = skipped.** An entry whose `fasta` is empty after the species-preset fallback is dropped from the DAG entirely: no index, no stage, no counts, and no columns in the matrices or the cascade summary. This is how the optional `mRNA` and `rhizo` stages are switched off.
4. **mRNA from cDNA references.** The `mRNA` class is meant to align against mRNA/cDNA sequences (e.g. spliced transcripts from the genome annotation), catching reads derived from mature transcripts before the genome step.
5. **rhizo = exogenous index.** `rhizo` carries the AM-fungus / other exogenous reference in the rice+microbiome design; it is optional like every other class.
6. **Genome alignment is independent.** `genome_align` uses the **trimmed reads** directly, not the cascade remainder — faithful to the legacy script. The cascade SAMs classify reads against known classes; the genome SAM is a complete trimmed-read alignment for genome-based analyses; both are produced.

## Environment setup

Three ways to get an environment (pick one; details in [docs/user-guide.md](docs/user-guide.md) §1):

1. **Fresh server**: `mamba env create -f workflow/environment.yaml` (all-in-one environment `srna-seq`, pinning snakemake-minimal 7.32.4 / bowtie 1.3.1 / trim-galore 0.6.10 / fastqc 0.11.9 / multiqc 1.21 / mirdeep-p2 1.1.4, plus the R/DESeq2 stack used by the optional differential-expression stage);
2. **Reuse an existing conda environment**: in `config/software.yaml` set `environment.type: conda` + `conda_prefix` (recommended) or `conda_name`; `run.sh` injects the prefix's `bin` into PATH automatically, no activate needed;
3. **System mode**: with `environment.type: system`, all tools come from PATH.

Environments are created explicitly by the user; Snakemake never deploys them automatically. Before launching, preflight with `bash run.sh -P <workdir> --check-software` (required executables: bowtie / bowtie-build / trim_galore / fastqc / multiqc / python3).

Snakemake version matrix: **7.32.4** is the reference version (pinned in `workflow/environment.yaml`; the four cluster profiles target the 7.x classic `--cluster` interface). Snakemake 8.x moved cluster submission to the executor plugin system — parsing and dry-runs work, but test before real cluster runs.

## Quick start

The full walkthrough is in [docs/user-guide.md](docs/user-guide.md) §2; a ready-to-edit project lives in [example/](example/README.md). Summary:

```bash
# 1) Working directory and data (SE fastq: {sample}.fastq.gz or {sample}_R1.fq.gz)
mkdir -p ~/work/rice/1.rawdata
cp root_rep1.fq.gz root_rep2.fq.gz leaf_rep1.fq.gz leaf_rep2.fq.gz ~/work/rice/1.rawdata/

# 2) Sample table (sample_id column; add a group column for the DE stage) + project config
cp example/samples.csv ~/work/rice/samples.csv
cp example/config.yaml  ~/work/rice/config.yaml       # fill in the reference paths

# 3) Preflight -> dry-run -> run (from the repository root / this subproject)
bash run.sh -P ~/work/rice --check-software
bash run.sh -P ~/work/rice -n
bash run.sh -P ~/work/rice -j 10
# positional form:  bash run.sh ~/work/rice ~/work/rice/config.yaml 10
# cluster example:  bash run.sh -P ~/work/rice --profile pbs --queue workq --memory 16G --runtime 600
```

Key points:

- Raw inputs live in `1.rawdata/` at the working-directory root; every derived artifact goes under `results/` (rename via `results_dir` in config) in numbered stage dirs;
- Config chain: repository `config/config.yaml` defaults → project config (`-c` or positional; auto-detected `<workdir>/config.yaml`) → `config.local.yaml` in the working directory (auto-layered; later files win);
- `species: "osa"` selects the rice (IRGSP-1.0) preset from `config/species.yaml`: a cascade entry with an empty `fasta` inherits the preset path for that class, and the same fallback applies to `genome.fasta` (explicit non-empty keys win; `species: "none"` disables the fallback entirely);
- Cluster jobs are submitted with the per-rule resources declared in `config/resources.yaml` (`threads/mem_mb/runtime_min`); copy that file into the working directory as `resources.yaml` to override (run.sh auto-detects it), or override globally via `--memory`/`--runtime`. Full options in `run.sh --help`.

## Directory structure

```
srna-seq/
├── run.sh                    # unified launcher CLI (four profiles / auto detection / preflight / --unlock; see run.sh --help)
├── workflow/
│   ├── Snakefile             # single entry point (species presets + config layering + target aggregation)
│   ├── environment.yaml      # all-in-one conda environment template (pinned; created explicitly by the user)
│   ├── rules/                # common / upstream / cascade / quant / deg / novel_mirna / meta
│   ├── scripts/              # runtime_config.py, collect_versions.py, count_features.py, merge_counts.py, cascade_summary.py, run_deseq2.R
│   ├── profile/              # default / pbs / sge / slurm profiles + README (cluster commands and pinned params)
│   └── multiqc_config.yaml
├── config/
│   ├── config.yaml           # repository defaults (osa preset paths; outputs rooted at results/)
│   ├── config.template.yaml  # annotated project template (copy into the working directory)
│   ├── species.yaml          # species presets (osa: per-class sncRNA fastas + genome)
│   ├── resources.yaml        # per-rule scheduler resources (threads/mem_mb/runtime_min)
│   ├── software.yaml         # unified software runtime (conda_prefix / system + r: section for the DE stage)
│   └── samples.csv           # sample table template (single sample_id column)
├── tests/                    # run_test.sh (dry-run / --real-run / --deg regression) / lint.sh / make_testdata.py / test_common.py
├── example/                  # example project (legacy rice sample table + filled config + start guide)
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
├── 1.rawdata/               # raw SE fastq inputs (the only data directory outside results/)
└── results/                 # ALL derived artifacts (rename via results_dir)
    ├── 0.index/             # per-class bowtie indices + genome index
    ├── 2.cleandata/         # trimmed fastq + trimming reports + FastQC
    ├── 3.align/             # filter/{class}/ cascade stages + genome/ alignment
    ├── 4.expression/        # per-class counts/RPM matrices + all_classes_counts.tsv
    ├── 5.QC/                # cascade summary + multiqc_report.html + software_versions.yaml
    ├── 6.DEG/               # optional DE results per analyzed class (deg stage)
    ├── 6.novel_mirna/       # optional novel-miRNA predictions per sample (novel_mirna stage)
    └── logs/                # per-rule logs
```

## Results path quick reference

All derived artifacts live under `results/` in the working directory (rename via `results_dir` in config); raw inputs `1.rawdata/` stay at the working-directory root. `{sample}` is a `sample_id` from the sample table; `{class}` is a configured cascade class.

| Result | Path |
|---|---|
| Class bowtie index | `results/0.index/{class}/{class}.1.ebwt` (+ `.2.ebwt`, `.rev.1.ebwt`, `.rev.2.ebwt` siblings) |
| Genome bowtie index | `results/0.index/genome/genome.1.ebwt` (+ siblings) |
| Trimmed reads | `results/2.cleandata/{sample}_trimmed.fq.gz` |
| Trim Galore report | `results/2.cleandata/{sample}_trimming_report.txt` |
| FastQC on trimmed reads | `results/2.cleandata/fastqc/{sample}_trimmed_fastqc.zip` (+ `.html`) |
| Cascade stage SAM | `results/3.align/filter/{class}/{sample}.sam` |
| Cascade stage unmapped reads | `results/3.align/filter/{class}/{sample}_unmapped.fq` |
| Genome alignment SAM | `results/3.align/genome/{sample}.sam` |
| Genome unmapped reads | `results/3.align/genome/{sample}_unmapped.fq` |
| Per-sample per-class counts | `results/4.expression/{class}/{sample}_counts.txt` |
| Per-class count matrix | `results/4.expression/{class}/{class}_counts.tsv` |
| Per-class RPM matrix | `results/4.expression/{class}/{class}_RPM.tsv` |
| Combined count matrix | `results/4.expression/all_classes_counts.tsv` |
| Raw read counts | `results/5.QC/raw_counts/{sample}.txt` |
| Cascade read-fate summary | `results/5.QC/cascade_summary.tsv` |
| Combined QC report | `results/5.QC/multiqc/multiqc_report.html` |
| Software version record | `results/5.QC/software_versions.yaml` |
| DE contrast table | `results/6.DEG/{class}/{treat}_vs_{control}_DESeq2.output.tsv` (deg stage) |
| DE plots | `results/6.DEG/{class}/*_VolcanoPlot.pdf`, `*_MAPlot.pdf`, `{class}_DESeq2.normalized.vst.PCA_plot.pdf`, `{class}_DESeq2.normalized.vst.Pearson_heatmap.pdf` (deg stage) |
| Novel-miRNA predictions | `results/6.novel_mirna/{sample}/{sample}_filter_P_prediction` (+ `.bed`) (novel_mirna stage) |
| Novel-miRNA tool tree | `results/6.novel_mirna/{sample}/` — precursor fasta/structures, signatures, raw `{sample}_predictions` (novel_mirna stage) |
| Per-rule logs | `results/logs/` |

`cascade_summary.tsv` has one row per sample with columns `sample`, `raw`, `trimmed`, then `{class}_mapped` / `{class}_unmapped` for each configured class, then `genome_unmapped`:

- `{class}_mapped` — reads assigned to the class in that stage (the `__mapped_total` line of `{sample}_counts.txt`, i.e. primary alignments in the stage SAM);
- `{class}_unmapped` — reads that passed through to the next stage (line count / 4 of `{sample}_unmapped.fq`);
- `genome_unmapped` — trimmed reads that did not map to the genome;
- `-1` — the corresponding unmapped FASTQ was missing (an upstream stage failed; the value could not be counted);
- `NA` — `genome_unmapped` only: the genome FASTA is not configured, so no genome alignment ran.

## QC notes

The MultiQC report aggregates, per sample:

- **Trim Galore** (`*_trimming_report.txt`): adapter-contaminated fraction and the length distribution after trimming;
- **FastQC** (`*_trimmed_fastqc.zip`): per-base quality, GC content, and other module checks on the trimmed reads;
- **Cascade read fate** (`cascade_summary_mqc.tsv`): a MultiQC custom-content **bargraph** ("sRNA cascade read fate") plotting `raw → trimmed → {class}_mapped → {class}_unmapped → genome_unmapped` per sample — the at-a-glance view of where every read went.

There is no bowtie-specific MultiQC module: per-stage mapping behaviour is judged from the cascade summary numbers. Pipeline health is read as the chain raw → trimmed yield → per-class assignment → genome-unmapped remainder (§7 of the user guide).

## Differential expression (optional, default off)

The `deg` stage (v0.2.0) runs DESeq2 over the per-class count matrices, ported from the rna-seq workflow. It is off unless you ask for it:

1. Add a `group` column to the sample table (and optionally `batch`); accepted headers are exactly `sample_id`, `sample_id,group`, or `sample_id,group,batch`. Single-column tables keep working with the stage off.
2. Set `deg.enabled: true` in the project config; tune `classes` (cascade classes to analyze), `control_group`, `foldchange`, `padj`, `batch_correction` (`"T"` fits `~ batch + group`), and `pca_ntop`.
3. `results/6.DEG/{class}/` then receives one contrast table per `<treat>_vs_<control>` comparison (first column `feature`), volcano + MA plots per contrast, a vst PCA plot and sample Pearson heatmap, and `sessionInfo.txt`.

Design rules enforced at parse time: the group column must exist when the stage is enabled, `control_group` must be one of the observed groups, at least two distinct groups are required, and a single replicate per group only warns (legal but weak). Details in [docs/user-guide.md](docs/user-guide.md) §3/§4.2/§8.

## Novel miRNA discovery (optional, default off)

The `novel_mirna` stage (v0.2.0) runs miRDeep-P2 (bioconda `mirdeep-p2=1.1.4`) per sample over the trimmed reads: candidate hairpins are excised from the reference genome and scored by mapping signature. It is off unless you ask for it:

1. Configure the genome (`genome.fasta`) — candidate hairpins are excised from it.
2. Set `novel_mirna.enabled: true` and `novel_mirna.mature_fasta` (a known mature-miRNA fasta, e.g. the miRBase mature set of the species) in the project config.
3. `results/6.novel_mirna/{sample}/` then receives the tool-native result tree; the main table is `{sample}_filter_P_prediction` (plant-criteria-filtered, redundancy-removed predictions; may legitimately be empty), with a `.bed` coordinate export, precursor fasta/structures, and raw `{sample}_predictions` scores. A `flag.log` marker tracks completion in the DAG.

Enabling the stage requires `miRDP2-v1.1.4_pipeline.bash` on PATH (pinned in `workflow/environment.yaml`); the preflight does not check it while the stage is off. Two verified packaging facts worth knowing: miRDP2 1.1.4 has no user-facing mature-reference flag (known-miRNA filtering uses the bundled plant mature index — the configured `mature_fasta` is enforced for provenance, not passed to the tool), and the rfam ncRNA filter silently no-ops because the package omits its index. Details in [docs/user-guide.md](docs/user-guide.md) §9; the full verified command surface is recorded in the header of `workflow/rules/novel_mirna.smk`.

## Differences from the legacy pipeline

The workflow is a refactor of the legacy `srna-seq.sh` PBS script (sequential bowtie cascade over rice sncRNA classes). The Trim Galore parameters, the bowtie invocation shape, and the trimmed-reads genome step are kept faithful; five deliberate deviations:

1. **New counting/quantification stage.** The legacy script produced SAM files only (its only counting hint was a commented-out `cut -f3 | sort | uniq -c` line); the workflow adds per-class feature counting, per-class count/RPM matrices, and the combined `all_classes_counts.tsv`.
2. **Per-class bowtie indices auto-built.** The legacy flow required manually pre-building every `0.index/<class>` bowtie index before running; `rule bowtie_index` builds each configured class index (and the genome index) from the configured FASTAs.
3. **Data-driven class list.** The legacy script hard-coded seven sncRNA stages plus the genome step; the ordered `cascade:` config list defines any number of classes with arbitrary names, in any order.
4. **mRNA / rhizo stages optional via empty fasta.** The legacy script unconditionally ran every stage; here an entry whose `fasta` is empty after the species-preset fallback is skipped entirely (no index, no stage) — and this applies to every class, not just `mRNA` and `rhizo`.
5. **MultiQC + cascade summary new.** The legacy run left per-stage console logs only; the workflow adds `cascade_summary.tsv` (per-sample read fate), the MultiQC custom bargraph, and `software_versions.yaml` provenance.

## Development and testing

```bash
make check    # bash -n syntax checks (no snakemake needed)
make lint     # static check suite tests/lint.sh (missing optional tools are skipped)
make test     # CI-equivalent full check (= check + lint)
```

Regression tests (need snakemake):

```bash
bash tests/run_test.sh               # synthetic-data dry-run: generate -> assemble working directory -> validate DAG integrity
bash tests/run_test.sh --deg         # dry-run the optional DE stage (group column + deg enabled; asserts deg_deseq2 in the DAG)
bash tests/run_test.sh --novel-mirna # dry-run the optional novel-miRNA stage (asserts novel_mirna in the DAG)
bash tests/run_test.sh --real-run    # end-to-end run + output assertions (server validation; needs the full analysis environment)
```

A pytest unit suite covers the sample-table loader, config/deg validation, and the species-merge parity: `python -m pytest tests -q` (needs pytest; no snakemake required for the extracted common.smk logic).

Test data is generated by `tests/make_testdata.py` with a fixed seed (2 x 20 kb chromosomes, rRNA/tRNA/miRNA reference fastas, 2 samples of 21-26 nt SE reads) and is never committed. The dry-run baseline DAG is **27 jobs** (2 samples x 9 sample-local jobs (trim, raw count, 3 cascade stages, genome alignment, 3 counts) + 4 index builds (rRNA/tRNA/miRNA/genome) + merge_counts + cascade_summary + multiqc + software_versions + rule all) — compare against this count after refactors (recorded in [CHANGELOG.md](CHANGELOG.md)).

The regression config runs with `species: "none"`: no `config/species.yaml` preset merge happens, and the four unconfigured classes (snoRNA/snRNA/mRNA/rhizo) are dropped through the empty-fasta skip path. This is also the supported way to run a fully custom class set — with a real species preset, an empty fasta inherits the preset path instead of skipping.

A lightweight CI job (lint + `--reads 2000` dry-run regression) is planned as part of Phase D integration (plan Task D3); until it lands, use `make test` + `bash tests/run_test.sh --reads 2000` locally.

## Known limitations and TODOs

Mirrors [docs/TODO.md](docs/TODO.md):

1. ~~Differential expression~~ — done in v0.2.0: the optional, default-off `deg` stage (DESeq2 over the per-class count matrices, group/batch sample-table design).
2. ~~Novel miRNA discovery~~ — done in v0.2.0: the optional, default-off `novel_mirna` stage (miRDeep-P2 per sample over the trimmed reads against the genome).
3. **Unit-test coverage**: `tests/test_common.py` covers the sample-table loader, config/deg validation, and the species-merge parity; the counting scripts (`count_features.py`, `merge_counts.py`, `cascade_summary.py`) are still uncovered.

## License

[Apache-2.0](LICENSE) © 2026 ChengYu

## Versions

v0.2.0 (2026-09-08): optional, default-off differential-expression stage — DESeq2 over the per-class count matrices with an optional group/batch sample-table design (`6.DEG/`), R plumbing in `software.yaml`/`environment.yaml`, `--deg` dry-run scenario, and loader/validation coverage in `tests/test_common.py`; plus the optional, default-off novel-miRNA discovery stage — miRDeep-P2 (`mirdeep-p2=1.1.4`) per sample over the trimmed reads against the genome (`6.novel_mirna/`), `novel_mirna:` config section, and the `--novel-mirna` dry-run scenario (see [CHANGELOG.md](CHANGELOG.md)).

v0.1.0 (2026-09-05): initial release — config layer, trim / bowtie-index / cascade / genome-alignment / counting rules, cascade summary + MultiQC, unified launcher, regression tests, and docs (see [CHANGELOG.md](CHANGELOG.md)).
