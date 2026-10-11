"""Pin the merge_featurecounts.py semantics: raw featureCounts inputs, gene
order anchored to the sample-id list, one sample column per file, summary
status aggregation, strict sample-name agreement."""
import pytest

import merge_featurecounts as mf


def write_fc(path, sample, genes):
    rows = ["# Program: featureCounts v2.1.1; Command: test", ""]
    rows.append("Geneid\tChr\tStart\tEnd\tStrand\tLength\t" + sample)
    for gid, n in genes:
        rows.append(f"{gid}\tchr1\t100\t200\t+\t100\t{n}")
    path.write_text("\n".join(rows) + "\n")


def write_summary(path, assigned, unassigned):
    rows = ["Status\t" + "x.sample", f"Assigned\t{assigned}", f"Unassigned_NoFeatures\t{unassigned}"]
    path.write_text("\n".join(rows) + "\n")


def test_merge_produces_matrix_and_stats(tmp_path):
    write_fc(tmp_path / "s1.featureCounts.tsv", "s1", [("g1", "10"), ("g2", "20")])
    write_fc(tmp_path / "s2.featureCounts.tsv", "s2", [("g1", "5"), ("g2", "1")])
    write_summary(tmp_path / "s1.featureCounts.tsv.summary", "30", "0")
    write_summary(tmp_path / "s2.featureCounts.tsv.summary", "6", "0")
    (tmp_path / "sample_ids.txt").write_text("s1\ns2\n")

    import os
    os.chdir(tmp_path)
    import sys
    sys.argv = ["merge_featurecounts.py", ".", "sample_ids.txt"]
    mf.main()

    matrix = (tmp_path / "count.matrix.tsv").read_text().strip().split("\n")
    assert matrix[0] == "gene_id\ts1\ts2"
    assert matrix[1] == "g1\t10\t5"
    assert matrix[2] == "g2\t20\t1"

    stats = (tmp_path / "count_stats.tsv").read_text().strip().split("\n")
    assert stats[0] == "sample\tAssigned\tUnassigned_NoFeatures"
    assert stats[1] == "s1\t30\t0"
    assert stats[2] == "s2\t6\t0"


def test_missing_sample_column_fails(tmp_path):
    write_fc(tmp_path / "s1.featureCounts.tsv", "OTHER", [("g1", "10")])
    (tmp_path / "sample_ids.txt").write_text("s1\n")
    import sys
    sys.argv = ["merge_featurecounts.py", str(tmp_path), str(tmp_path / "sample_ids.txt")]
    with pytest.raises(SystemExit):
        mf.main()


def test_gene_union_fills_na(tmp_path):
    write_fc(tmp_path / "s1.featureCounts.tsv", "s1", [("g1", "10")])
    write_fc(tmp_path / "s2.featureCounts.tsv", "s2", [("g1", "5"), ("g2", "1")])
    write_summary(tmp_path / "s1.featureCounts.tsv.summary", "10", "0")
    write_summary(tmp_path / "s2.featureCounts.tsv.summary", "6", "0")
    (tmp_path / "sample_ids.txt").write_text("s1\ns2\n")
    import sys
    sys.argv = ["merge_featurecounts.py", str(tmp_path), str(tmp_path / "sample_ids.txt")]
    mf.main()
    matrix = (tmp_path / "count.matrix.tsv").read_text().strip().split("\n")
    assert matrix[2] == "g2\tNA\t1"
