# seclip-seq backlog (TODO)

> Items identified during v0.1 and deferred to a later release; none are scheduled yet.

## 1. Cross-sample reproducible peaks

- v0.1 calls peaks per sample (PureCLIP, plus CLIPper when configured). A reproducible peak set across replicates (IDR-style ranking or overlap-based merging of `results/5.callpeak/*.bed`) is out of scope for v0.1 and is the first v0.2 candidate.

## 2. Peak annotation

- Annotate called peaks against gene models / repeat features (GTF-overlap based, or a HOMER/ChIPseeker-style step) so every `results/5.callpeak/*.bed` ships with nearest-gene and biotype tables.

## 3. IP vs input control

- v0.1 has no input/background channel concept. Supporting an paired input control (PureCLIP background estimation from a control BAM, or CLIPper contrast modes) needs a sample-table design decision first (how to declare IP/input pairs).

## 4. ea-utils fallback for fastq-sort

- `rule fastq_sort` depends on `fastq-sort` from ea-utils, which has few conda builds (see the note in `workflow/environment.yaml`). If a server cannot solve ea-utils, switch the rule to an equivalent name-based sort (e.g. `seqkit sort -N` or a sort-by-id one-liner) and record the change here and in the CHANGELOG.
