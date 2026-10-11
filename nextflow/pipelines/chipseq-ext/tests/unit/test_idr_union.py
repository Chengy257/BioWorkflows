"""Unit tests for idr_union.py: pair-name parsing, per-chromosome interval
merging, distinct-replicate support counting and the narrowPeak output
contract (pinned after the real-execution debugging rounds - see
nextflow/docs/chipseq-ext-validation-notes.md)."""
import sys

import pytest

import idr_union


def write_idr(path, rows):
    # idr's *idrValues.txt is narrowPeak-like: comment lines allowed, no
    # column header row (a header line would fail the int() coordinate parse)
    lines = ["# savedContent of idr"]
    for chrom, start, end, score in rows:
        lines.append(f"{chrom}\t{start}\t{end}\t1000\t{score}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_union(monkeypatch, out, files, min_replicates="2"):
    argv = ["idr_union.py", "--group", "g", "--min-replicates", min_replicates,
            "--out", str(out)] + [str(p) for p in files]
    monkeypatch.setattr(sys, "argv", argv)
    idr_union.main()


def test_parse_pair():
    assert idr_union.parse_pair("G__rep1_vs_rep2.idrValues.txt") == ("rep1", "rep2")


def test_parse_pair_rejects_bad_name(tmp_path):
    bad = tmp_path / "somethingelse.txt"
    bad.write_text("")
    with pytest.raises(SystemExit):
        idr_union.parse_pair(str(bad))


def test_merge_overlapping_and_touching():
    merged = idr_union.merge([(chr_, s, e) for chr_, s, e in [
        ("chr1", 100, 200), ("chr1", 150, 300), ("chr1", 300, 400),
        ("chr2", 50, 60), ("chr1", 500, 600)]])
    assert merged == {"chr1": [(100, 400), (500, 600)], "chr2": [(50, 60)]}


def test_overlap_len():
    assert idr_union.overlap_len(0, 100, 50, 150) == 50
    assert idr_union.overlap_len(0, 100, 100, 200) == 0
    assert idr_union.overlap_len(0, 100, 200, 300) == 0


def test_support_counting_two_pairs(tmp_path, monkeypatch):
    # one interval covered by two pairs -> 4 distinct replicates; one interval
    # covered by a single pair -> 2 distinct replicates (documented semantics)
    p1 = tmp_path / "g__r1_vs_r2.idrValues.txt"
    p2 = tmp_path / "g__r3_vs_r4.idrValues.txt"
    write_idr(p1, [("chr1", 1000, 1500, 3.5), ("chr1", 5000, 5500, 2.0)])
    write_idr(p2, [("chr1", 1200, 1700, 4.5)])
    out = tmp_path / "out.narrowPeak"
    run_union(monkeypatch, out, [p1, p2])
    rows = [line.split("\t") for line in out.read_text().splitlines()[1:]]
    assert len(rows) == 2
    # merged union of the two overlapping intervals carries the max IDR score
    assert rows[0][0:3] == ["chr1", "1000", "1700"]
    assert rows[0][4] == "4"          # support of r1, r2, r3, r4
    assert rows[0][6] == "4.5"        # signalValue = max overlapping IDR score
    assert rows[1][0:3] == ["chr1", "5000", "5500"]
    assert rows[1][4] == "2"          # single-pair interval


def test_min_replicates_filters(tmp_path, monkeypatch):
    p1 = tmp_path / "g__r1_vs_r2.idrValues.txt"
    p2 = tmp_path / "g__r3_vs_r4.idrValues.txt"
    write_idr(p1, [("chr1", 1000, 1500, 3.5)])
    write_idr(p2, [("chr1", 1200, 1700, 4.5)])
    out = tmp_path / "out.narrowPeak"
    run_union(monkeypatch, out, [p1, p2], min_replicates="5")
    assert out.read_text().splitlines()[1:] == []
