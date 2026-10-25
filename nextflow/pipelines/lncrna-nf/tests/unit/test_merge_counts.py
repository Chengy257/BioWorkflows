"""Unit tests for merge_lncrna_counts.py: the raw featureCounts merge
(gene order anchored to the first sample, header sample-column check
against the BAM-basename labels featureCounts really emits), and the
TPM/FPKM arithmetic from the Length column + Assigned reads."""
import sys

import merge_lncrna_counts as mlc


def write_sample(d, sid, genes, header=None):
    # genes: list of (gene_id, length, count); featureCounts labels the
    # count column with the input BAM basename, not the bare sample id
    # (the S1 ground truth; the bare form stays accepted for hand-made runs)
    if header is None:
        header = f"{sid}.markdup.sorted.bam"
    with open(d / f"{sid}.featureCounts.tsv", "w") as fh:
        fh.write("# Program: featureCounts (fixture)\n")
        fh.write("Geneid\tChr\tStart\tEnd\tStrand\tLength\t" + header + "\n")
        for g, length, c in genes:
            fh.write(f"{g}\tchr1\t100\t{99 + length}\t+\t{length}\t{c}\n")
    assigned = sum(c for (_g, _l, c) in genes) + 10
    with open(d / f"{sid}.featureCounts.tsv.summary", "w") as fh:
        fh.write("Status\tx\n")
        fh.write(f"Assigned\t{assigned}\n")
        fh.write("Unassigned_Unmapped\t0\n")


def run(tmp_path, ids_text):
    ids = tmp_path / "ids.txt"
    ids.write_text(ids_text, encoding="utf-8")
    old = sys.argv
    sys.argv = ["merge_lncrna_counts.py", str(tmp_path), str(ids), "lncrna."]
    try:
        mlc.main()
    finally:
        sys.argv = old


def test_merge_and_matrices(tmp_path):
    write_sample(tmp_path, "w1", [("g1", 100, 10), ("g2", 300, 20)])
    write_sample(tmp_path, "w2", [("g1", 100, 5), ("g2", 300, 0)])
    run(tmp_path, "w1\nw2\n")

    count = (tmp_path / "lncrna.count.matrix.tsv").read_text().splitlines()
    assert count[0].split("\t") == ["gene_id", "w1", "w2"]
    assert count[1].split("\t") == ["g1", "10", "5"]       # order anchored to w1
    assert count[2].split("\t") == ["g2", "20", "0"]

    tpm = (tmp_path / "lncrna.tpm.matrix.tsv").read_text().splitlines()
    # w1: rpk = 10/(100/1000)=100, 20/(300/1000)=66.667; sum=166.667
    # g1 tpm = 100/166.667*1e6 = 600000
    w1_g1 = float(tpm[1].split("\t")[1])
    assert abs(w1_g1 - 600000.0) < 1.0

    fpkm = (tmp_path / "lncrna.fpkm.matrix.tsv").read_text().splitlines()
    # w1 assigned = 40; g1 fpkm = 100/40*1e6 = 2.5e6
    w1_g1 = float(fpkm[1].split("\t")[1])
    assert abs(w1_g1 - 2500000.0) < 1.0

    stats = (tmp_path / "lncrna.count_stats.tsv").read_text().splitlines()
    assert stats[1].split("\t")[0] == "Assigned"
    assert stats[1].split("\t")[1] == "40"


def test_bare_id_header_still_accepted(tmp_path):
    write_sample(tmp_path, "w1", [("g1", 100, 10)], header="w1")
    run(tmp_path, "w1\n")
    assert (tmp_path / "lncrna.count.matrix.tsv").read_text().splitlines()[1] \
        .split("\t") == ["g1", "10"]


def test_header_mismatch_fails(tmp_path):
    write_sample(tmp_path, "w1", [("g1", 100, 10)])
    # tamper: count column carries ANOTHER sample's BAM name
    p = tmp_path / "w1.featureCounts.tsv"
    p.write_text(p.read_text().replace(
        "\tw1.markdup.sorted.bam\n", "\tw2.markdup.sorted.bam\n"),
        encoding="utf-8")
    ids = tmp_path / "ids.txt"
    ids.write_text("w1\n", encoding="utf-8")
    old = sys.argv
    sys.argv = ["merge_lncrna_counts.py", str(tmp_path), str(ids)]
    try:
        try:
            mlc.main()
            assert False, "expected a sample-mismatch failure"
        except SystemExit:
            pass
    finally:
        sys.argv = old
