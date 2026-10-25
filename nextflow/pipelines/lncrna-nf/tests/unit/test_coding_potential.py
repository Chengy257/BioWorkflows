"""Unit tests for coding_potential_intersect.py: the D3 intersection rule
(CPC2 noncoding AND CPAT noncoding AND exonic length >= min), the CPAT 3.x
best-ORF table layout (one row per sequence; no_ORF ids are noncoding; a
malformed header is a hard error), and the exonic-length arithmetic across
multi-exon transcripts."""
import sys

import coding_potential_intersect as cpi

CPC2 = "\n".join([
    "# CPC2 result",
    "ID\tlen\tORF\tFitch\tpI\tcharge\tcoding_prob\tcoding",
    "novA\t300\t60\t0.1\t7.0\t0\t0.1\tnoncoding",
    "novB\t400\t90\t0.2\t7.2\t0\t0.2\tnoncoding",
    "novC\t500\t450\t0.9\t9.0\t10\t0.9\tcoding",
    "novD\t100\t30\t0.1\t7.0\t0\t0.1\tnoncoding",
    "novE\t300\t0\t0.1\t7.0\t0\t0.1\tnoncoding",
]) + "\n"

# CPAT 3.x <prefix>.ORF_prob.best.tsv: exactly one row per sequence with an
# ORF; novE has no ORF at all (absent here, listed in cpat.no_ORF.txt).
CPAT = "\n".join([
    "Sequence_ID\tKeytranscript_ID\tORF_ID\tORF_Begin\tORF_End\tORF_Length\t"
    "Fickett_Score\tHexamer_Score\tCoding_Probability",
    "novA\tnovA\tnovA_1\t1\t60\t60\t0.1\t0.1\t0.1",
    "novB\tnovB\tnovB_1\t1\t90\t90\t0.2\t0.2\t0.2",
    "novC\tnovC\tnovC_1\t1\t450\t450\t0.9\t0.9\t0.9",
    "novD\tnovD\tnovD_1\t1\t30\t30\t0.1\t0.1\t0.1",
]) + "\n"

NO_ORF = "novE\n"

GTF = "\n".join([
    "chr1\ttest\texon\t1000\t1400\t.\t+\t.\tgene_id \"novA\"; transcript_id \"novA\";",   # 401 nt
    "chr1\ttest\texon\t2000\t2100\t.\t+\t.\tgene_id \"novB\"; transcript_id \"novB\";",   # 101
    "chr1\ttest\texon\t2200\t2400\t.\t+\t.\tgene_id \"novB\"; transcript_id \"novB\";",   # 201 -> 302 total
    "chr1\ttest\texon\t3000\t3500\t.\t+\t.\tgene_id \"novC\"; transcript_id \"novC\";",   # 501
    "chr1\ttest\texon\t4000\t4050\t.\t+\t.\tgene_id \"novD\"; transcript_id \"novD\";",   # 51
    "chr1\ttest\texon\t4500\t4799\t.\t+\t.\tgene_id \"novE\"; transcript_id \"novE\";",   # 300
    "chr1\ttest\ttranscript\t2000\t2400\t.\t+\t.\tgene_id \"novB\"; transcript_id \"novB\";",
]) + "\n"


def run(tmp_path, min_length, no_orf=True):
    cpc2 = tmp_path / "cpc2.txt"
    cpc2.write_text(CPC2, encoding="utf-8")
    cpat = tmp_path / "cpat.ORF_prob.best.tsv"
    cpat.write_text(CPAT, encoding="utf-8")
    no_orf_file = tmp_path / "cpat.no_ORF.txt"
    if no_orf:
        no_orf_file.write_text(NO_ORF, encoding="utf-8")
    gtf = tmp_path / "novel.gtf"
    gtf.write_text(GTF, encoding="utf-8")
    out_ids = tmp_path / "ids.txt"
    out_gtf = tmp_path / "nc.gtf"
    old = sys.argv
    sys.argv = ["coding_potential_intersect.py", "--cpc2", str(cpc2),
                "--cpat", str(cpat), "--cpat-no-orf", str(no_orf_file),
                "--gtf", str(gtf),
                "--min-length", str(min_length),
                "--out-ids", str(out_ids), "--out-gtf", str(out_gtf)]
    try:
        cpi.main()
    finally:
        sys.argv = old
    return out_ids.read_text().split(), out_gtf.read_text().splitlines()


def test_intersection_and_length(tmp_path):
    ids, lines = run(tmp_path, 200)
    # novA: noncoding both + 401 nt -> kept; novB: noncoding both + 302 nt -> kept;
    # novC: coding -> dropped; novD: noncoding but 51 nt -> dropped;
    # novE: no ORF (CPAT verdict: noncoding) + 300 nt -> kept
    assert ids == ["novA", "novB", "novE"]
    # the input carries novA exon x1 + novB exon x2 + novB transcript x1 + novE exon x1
    assert len(lines) == 5


def test_length_floor_drops_short(tmp_path):
    ids, _ = run(tmp_path, 350)
    assert ids == ["novA"]


def test_no_orf_ids_are_noncoding(tmp_path):
    # without the no_ORF list novE is unknown to CPAT and drops out
    ids, _ = run(tmp_path, 200, no_orf=False)
    assert ids == ["novA", "novB"]


def test_cpat_missing_header_fails(tmp_path):
    bad = tmp_path / "cpat.bad.tsv"
    bad.write_text("seq_id\tfoo\nnovA\t1\n", encoding="utf-8")
    try:
        cpi.read_cpat(str(bad))
        assert False, "expected a hard failure on the missing header"
    except SystemExit:
        pass


def test_exonic_length_math(tmp_path):
    (tmp_path / "novel.gtf").write_text(GTF, encoding="utf-8")
    lengths = cpi.exonic_lengths(str(tmp_path / "novel.gtf"))
    assert lengths["novA"] == 401
    assert lengths["novB"] == 302
    assert lengths["novE"] == 300
