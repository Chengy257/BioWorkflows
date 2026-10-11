"""Unit tests for organelle_summary.py and the spike-in helpers: the contig
pattern convention (<= 3 chars match exactly case-insensitively, longer
patterns match as substrings), the unmapped '*' row exclusion, and the
scale-factor arithmetic (the spike contig substring rule is pinned here -
`spike` must match `spike1`, the bug that broke the first fixture round)."""
import sys

import organelle_summary
import spikein_factor
import spikein_summary

IDXSTATS = ("chr1\t300000\t900\t0\n"
            "spike1\t1000\t10\t0\n"
            "chrc\t150000\t50\t0\n"
            "MTR2\t8000\t30\t0\n"
            "mitochondrion_A\t8000\t20\t0\n"
            "*\t*\t5\t0\n")


def run_organelle(tmp_path, monkeypatch, patterns, idx_text=IDXSTATS):
    idx = tmp_path / "s1_idxstats.tsv"
    idx.write_text(idx_text, encoding="utf-8")
    out = tmp_path / "organelle_summary.tsv"
    argv = ["organelle_summary.py", "--patterns", patterns, "--out", str(out),
            "--mqc", str(tmp_path / "mqc.tsv"), str(idx)]
    monkeypatch.setattr(sys, "argv", argv)
    organelle_summary.main()
    lines = out.read_text().splitlines()
    header = lines[0].split("\t")
    return dict(zip(header, lines[1].split("\t")))


def test_short_patterns_match_exactly(tmp_path, monkeypatch):
    row = run_organelle(tmp_path, monkeypatch, "chrc,chrm,pt,mt")
    # chrc (50) counts; MTR2 must NOT be swallowed by the 2-char 'mt' pattern;
    # mapped denominator = 900 + 10 + 50 + 30 + 20 = 1010
    assert row["organelle_reads"] == "50"
    assert row["organelle_fraction"] == "0.0495"
    assert "chrc:50" in row["matched_contigs"]


def test_long_patterns_match_substrings(tmp_path, monkeypatch):
    row = run_organelle(tmp_path, monkeypatch, "chrc,mitochondr")
    assert row["organelle_reads"] == "70"   # chrc 50 + mitochondrion_A 20


def test_unmapped_star_row_excluded(tmp_path, monkeypatch):
    # the '*' idxstats row (5 unmapped) must not enter the mapped denominator
    # (900 + 10 + 50 + 30 + 20 = 1010, not 1015)
    row = run_organelle(tmp_path, monkeypatch, "chrc")
    assert row["mapped_reads"] == "1010"


def test_fraction_zero_when_no_organelle(tmp_path, monkeypatch):
    row = run_organelle(tmp_path, monkeypatch, "pltd",
                        "chr1\t300000\t900\t0\n*\t*\t5\t0\n")
    assert row["organelle_fraction"] == "0.0000"


def run_spikein_summary(tmp_path, monkeypatch, idx_text=IDXSTATS):
    idx = tmp_path / "s1_idxstats.tsv"
    idx.write_text(idx_text, encoding="utf-8")
    out = tmp_path / "summary.tsv"
    argv = ["spikein_summary.py", "--patterns", "spike", "--name", "lambda",
            "--out", str(out), "--mqc", str(tmp_path / "mqc.tsv"), str(idx)]
    monkeypatch.setattr(sys, "argv", argv)
    spikein_summary.main()
    lines = out.read_text().splitlines()
    header = lines[0].split("\t")
    return dict(zip(header, lines[1].split("\t")))


def test_spike_pattern_matches_spike1(tmp_path, monkeypatch):
    # regression pin: the `spike` substring pattern must match `spike1`
    # (the chrS-vs-spike rename was the root cause of the all-NA summary)
    row = run_spikein_summary(tmp_path, monkeypatch)
    assert row["spike_mapped"] == "10"
    assert row["host_mapped"] == "1000"   # chrc/MTR2/mitochondrion are host
    assert row["spike_rate"] == "0.0099"
    assert row["scale_factor"] == "100000.0000"


def test_spikein_summary_na_for_missing_file(tmp_path, monkeypatch):
    out = tmp_path / "summary.tsv"
    argv = ["spikein_summary.py", "--patterns", "spike", "--name", "lambda",
            "--out", str(out), "--mqc", str(tmp_path / "mqc.tsv"),
            str(tmp_path / "missing_idxstats.tsv")]
    monkeypatch.setattr(sys, "argv", argv)
    spikein_summary.main()
    lines = out.read_text().splitlines()[1].split("\t")
    assert lines[1:] == ["NA", "NA", "NA", "NA"]


def test_spikein_factor_output(tmp_path, monkeypatch, capsys):
    idx = tmp_path / "s1_idxstats.tsv"
    idx.write_text(IDXSTATS, encoding="utf-8")
    argv = ["spikein_factor.py", "--patterns", "spike", str(idx)]
    monkeypatch.setattr(sys, "argv", argv)
    spikein_factor.main()
    assert capsys.readouterr().out == "s1\t100000.000000\n"


def test_spikein_factor_zero_spike_clamps_to_1e6(tmp_path, monkeypatch, capsys):
    idx = tmp_path / "s1_idxstats.tsv"
    idx.write_text("chr1\t300000\t900\t0\n*\t*\t5\t0\n", encoding="utf-8")
    argv = ["spikein_factor.py", "--patterns", "spike", str(idx)]
    monkeypatch.setattr(sys, "argv", argv)
    spikein_factor.main()
    assert capsys.readouterr().out == "s1\t1000000.000000\n"
