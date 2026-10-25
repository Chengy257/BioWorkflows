"""Unit tests for parse_classcode_u.py: the classcode-u extraction keeps a
novel transcript AND its exon rows (gffcompare writes class_code only on
transcript rows - exon rows carry the same transcript_id and no class_code;
the real-run fixture round caught the earlier exon-dropping variant)."""
import sys

import parse_classcode_u

# REAL gffcompare *.annotated.gtf layout: class_code on transcript rows only
GTF = "\n".join([
    "#!genome-build test",
    'chr1\ttest\ttranscript\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; class_code "=";',
    'chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; exon_number "1";',
    'chr1\ttest\ttranscript\t300\t500\t.\t+\t.\tgene_id "g2"; transcript_id "t2"; class_code "u";',
    'chr1\ttest\texon\t300\t400\t.\t+\t.\tgene_id "g2"; transcript_id "t2"; exon_number "1";',
    'chr1\ttest\texon\t400\t500\t.\t+\t.\tgene_id "g2"; transcript_id "t2"; exon_number "2";',
    'chr1\ttest\ttranscript\t600\t700\t.\t-\t.\tgene_id "g3"; transcript_id "t3"; class_code "x";',
    'chr1\ttest\texon\t600\t700\t.\t-\t.\tgene_id "g3"; transcript_id "t3"; exon_number "1";',
]) + "\n"


def test_classcode_u_extraction(tmp_path):
    gtf = tmp_path / "annotated.gtf"
    gtf.write_text(GTF, encoding="utf-8")
    out_gtf = tmp_path / "novel.gtf"
    out_ids = tmp_path / "ids.txt"
    argv = ["parse_classcode_u.py", "--gtf", str(gtf),
            "--out-gtf", str(out_gtf), "--out-ids", str(out_ids)]
    old = sys.argv
    sys.argv = argv
    try:
        parse_classcode_u.main()
    finally:
        sys.argv = old
    ids = out_ids.read_text().split()
    assert ids == ["t2"]  # deduped, only u
    lines = out_gtf.read_text().splitlines()
    # transcript row + BOTH exon rows (no class_code attr on exons)
    assert len(lines) == 3
    assert lines[0].split("\t")[2] == "transcript"
    assert sorted(l.split("\t")[2] for l in lines[1:]) == ["exon", "exon"]
    assert all('transcript_id "t2"' in l for l in lines)
    # the known and classcode-x transcripts are fully excluded
    assert not any('"t1"' in l or '"t3"' in l for l in lines)
