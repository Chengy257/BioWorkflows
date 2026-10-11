"""Unit tests for parse_classcode_u.py: the classcode-u extraction pins the
novel-transcript rule (only class_code='u' rows, transcript ids deduped)."""
import sys

import parse_classcode_u

GTF = "\n".join([
    "#!genome-build test",
    "chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\"; class_code \"=\";",
    "chr1\ttest\texon\t300\t400\t.\t+\t.\tgene_id \"g2\"; transcript_id \"t2\"; class_code \"u\";",
    "chr1\ttest\texon\t400\t500\t.\t+\t.\tgene_id \"g2\"; transcript_id \"t2\"; class_code \"u\";",
    "chr1\ttest\ttranscript\t300\t500\t.\t+\t.\tgene_id \"g2\"; transcript_id \"t2\"; class_code \"u\";",
    "chr1\ttest\texon\t600\t700\t.\t-\t.\tgene_id \"g3\"; transcript_id \"t3\"; class_code \"x\";",
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
    assert len(lines) == 3
    assert all('class_code "u"' in l for l in lines)
