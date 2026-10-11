# lncrna-nf validation notes

Validation record for the Stage-4 local lncRNA identification workflow
(`nextflow/pipelines/lncrna-nf/`). Companion to
`lncrna-scientific-decision-record.md` (ratified 2026-10-12). Dates are
absolute; every entry states what ran and what was observed.

## 2026-10-12 - ecosystem evidence + module availability

- nf-core/lncpipe: still "UNDER DEVELOPMENT" (dockstore v0.2.4; no stable
  release; absent from recent release newsletters). Not adoptable; watch
  item.
- Central nf-core modules verified present via the GitHub API:
  `stringtie` (assembly, container stringtie:3.0.3), `stringtie/merge`,
  `gffcompare`, `gffread`, `cpc2`, `diamond/blastx`, `subread/featurecounts`
  (reused from the rnaseq-ext vendored copy). NOT present: pfamscan, cpat.
- Containers verified on quay.io/biocontainers: `cpat:3.0.5--py310h20b60a1_4`,
  `pfam_scan:1.6--hdfd78af_5`.
- Coding-potential evidence for D3: CNCI python2-only and unmaintained;
  CPC2 maintenance-frozen but packaged; CPAT 3.x python3 and maintained.
  CNCI retired, CPAT adopted as the second independent method.

## 2026-10-12 - lncrna-nf local validation

- Config syntax: `nextflow config -flat -profile test` and `-profile
  test_full` exit 0 for `nextflow/pipelines/lncrna-nf/` under Nextflow
  26.04.6.
- Fixtures: `tests/make_fixtures.py` (fixed seed 20261012) mirrors the
  Stage-2 handoff ({aligner}/ with markdup BAM placeholders, 4 samples) and
  the reference set (300 kb chr1 genome(+.fai), 5-gene reference GTF,
  protein-coding-only GTF, CPAT model/hexamer placeholders, mini Pfam HMM
  file, empty nr.dmnd placeholder). Three planted novel loci: novA/novB
  noncoding (must be recovered), novC coding-like (must be filtered).
- Dry-run baselines: `-profile test` (homology off) resolves 12 process
  slots; `-profile test_full` (pfam + diamond on) resolves 14 slots over 14
  processes. Expected task counts at execution: 18 (test) / 20 (test_full)
  - 4x STRINGTIE + 4x FEATURECOUNTS are the per-sample slots; exact counts
  confirmed at the deferred real-container round.
- pytest: `tests/unit/` 8 tests green (ms env) pinning: classcode-u
  extraction (dedup, only-u rows), the D3 intersection (CPC2 AND CPAT AND
  length >= 200, exonic-length arithmetic over multi-exon transcripts),
  pfam_scan column semantics (id col 1, sequence E-value col 13), the
  final-set union + coordinate-sorted ref+novel merge + exon-spliced
  revcomp FASTA, and the count/TPM/FPKM matrix arithmetic from the Length
  column + Assigned reads (with a header-mismatch failure test). Full
  local sweep: `nextflow/tests/static-checks.sh` green under `nf` (config
  x4 pipelines, YAML 33 files, doc links 17 files, py_compile) and
  `--unit-only` under `ms` (67 tests: 27 chipseq + 19 rnaseq + 13 atacseq
  + 8 lncrna).

## Known limitations at this point (deferred, not blockers)

- DEFERRED EXECUTION VALIDATION (owner decision 2026-10-12, design
  amendment 8): real-container fixture execution of the whole chain (needs
  real fixture BAMs, the stringtie/gffcompare/gffread/cpc2/cpat/pfam_scan/
  diamond/subread images, and mini Pfam/DIAMOND databases built for real -
  `diamond makedb` + `hmmbuild` on toy seed alignments). The scientific
  gate at that round: planted noncoding transcripts recovered, planted
  coding transcripts excluded. versions.yml tool rows (stringtie 3.0.3,
  cpat 3.0.5, pfam_scan 1.6, diamond, cpc2, gffcompare, gffread,
  featureCounts) are added only after that validated run (ledger rule).
- CPAT species models: the human model ships with CPAT; a rice logit model
  must be user-provided (species neutrality). The deferred round should
  decide whether to build one with `cpat.py make_logitModel` (documented
  in the pipeline README).
- The DIAMOND NR filter defaults ON with a 24 h resource ceiling (legacy
  resource defaults); projects without a real nr.dmnd must disable it.
- Development-box note: this CentOS 7 machine intermittently served stale
  file contents to freshly spawned processes during this session
  (md5-verified writes invisible to next readers; one orphaned directory
  entry had to be renamed away before deletion worked). All files were
  md5-verified after writing; re-verify the tree on the HPC at the
  deferred execution round. One real finding of the debugging round: an
  include-name typo (TOBIAS_ATACCORRECT vs the upstream TOBIAS_ATACORRECT)
  initially presented as a confusing "process not defined" failure - the
  error message names the missing export exactly; read it literally.
