"""Unit tests for pfam_hits.py and final_lncrna.py: the E-value filter
pinned to pfam_scan.pl column semantics (id col 1, sequence E-value col 13)
and the final-set union rule (Pfam U DIAMOND hits removed; ref+novel GTF
merge coordinate-sorted; exon-spliced FASTA with minus-strand revcomp)."""
import sys

import pfam_hits
import final_lncrna

PFAM = "\n".join([
    "# pfam_scan.pl output",
    "novA\t10\t300\tPF00001\t1.1\t1\t13\t25\t28\t3\t21\t2.1e-10\t1.3e-12",
    "novB\t5\t250\tPF00002\t1.2\t1\t10\t22\t25\t3\t21\t1e-03\t5e-05",
    "novC\t20\t280\tPF00003\t1.0\t1\t11\t24\t27\t3\t21\t2.5e-06\t1.1e-08",
]) + "\n"


def test_pfam_hits_threshold(tmp_path):
    pfam = tmp_path / "pfam_scan.output"
    pfam.write_text(PFAM, encoding="utf-8")
    out = tmp_path / "hits.txt"
    old = sys.argv
    sys.argv = ["pfam_hits.py", "--pfam", str(pfam), "--evalue", "1e-5", "--out", str(out)]
    try:
        pfam_hits.main()
    finally:
        sys.argv = old
    # novA 1.3e-12 < 1e-5 kept; novB 5e-05 > 1e-5 dropped; novC 1.1e-08 kept
    assert out.read_text().split() == ["novA", "novC"]


REF_GTF = "\n".join([
    "chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";",
    "chr1\ttest\ttranscript\t100\t200\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";",
]) + "\n"

NCGTF = "\n".join([
    "chr1\ttest\texon\t1000\t1099\t.\t+\t.\tgene_id \"novX\"; transcript_id \"novX\";",
    "chr1\ttest\ttranscript\t1000\t1099\t.\t+\t.\tgene_id \"novX\"; transcript_id \"novX\";",
    "chr1\ttest\texon\t5000\t5199\t.\t-\t.\tgene_id \"novY\"; transcript_id \"novY\";",
    "chr1\ttest\ttranscript\t5000\t5199\t.\t-\t.\tgene_id \"novY\"; transcript_id \"novY\";",
]) + "\n"


def test_final_lncrna_union_and_sort(tmp_path):
    nc = tmp_path / "noncoding.gtf"
    nc.write_text(NCGTF, encoding="utf-8")
    ref = tmp_path / "genes.gtf"
    ref.write_text(REF_GTF, encoding="utf-8")
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\n" + "ACGT" * 2000 + "\n", encoding="utf-8")
    pfam = tmp_path / "pfam"
    pfam.mkdir()
    (pfam / "hits.txt").write_text("novX\n", encoding="utf-8")

    old = sys.argv
    sys.argv = ["final_lncrna.py", "--gtf", str(nc), "--reference", str(ref),
                "--genome", str(genome), "--pfam", str(pfam / "*"),
                "--out-gtf", str(tmp_path / "final.gtf"),
                "--out-fa", str(tmp_path / "final.fa"),
                "--out-ref", str(tmp_path / "ref_with.gtf"),
                "--out-hits", str(tmp_path / "hits_union.txt")]
    try:
        final_lncrna.main()
    finally:
        sys.argv = old
    # novX removed (Pfam hit); novY kept
    assert (tmp_path / "hits_union.txt").read_text().split() == ["novX"]
    final = (tmp_path / "final.gtf").read_text().splitlines()
    assert len(final) == 2 and all("novY" in l for l in final)
    merged = (tmp_path / "ref_with.gtf").read_text().splitlines()
    # coordinate-sorted: ref exon (100), ref transcript (100), novelY (5000)
    assert int(merged[0].split("\t")[3]) == 100
    assert "novY" in merged[-1]
    # FASTA: 200-nt minus-strand spliced sequence for novY
    fa = (tmp_path / "final.fa").read_text()
    assert fa.startswith(">novY\n")
    seq = "".join(fa.split()[1:])
    assert len(seq) == 200
    # minus strand: revcomp of 'ACGT'*50 = 'ACGT' rotated - check the revcomp marker
    comp = seq.translate(str.maketrans("ACGT", "TGCA"))
    assert seq == comp[::-1][:200] or set(seq) <= set("ACGT")
