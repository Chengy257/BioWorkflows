"""Unit tests for seacr_to_narrowpeak.py: 6-column SEACR bed to the 10-column
narrowPeak contract, score clamping, and the malformed-input error paths."""
import sys

import pytest

import seacr_to_narrowpeak

SRC = ("chr1\t1000\t2000\t2500.75\t500\t1200\n"
       "chr1\t5000\t5400\t0.4\t1\t5200\n"
       "chr1\t7000\t7600\t999.4\t100\t7300\n")


def test_convert_clamps_and_names(tmp_path):
    src = tmp_path / "in.bed"
    dst = tmp_path / "out.narrowPeak"
    src.write_text(SRC, encoding="utf-8")
    n = seacr_to_narrowpeak.convert(str(src), str(dst), "H3K27ac_WT")
    assert n == 3
    rows = [line.split("\t") for line in dst.read_text().splitlines()]
    assert rows[0][0:4] == ["chr1", "1000", "2000", "H3K27ac_WT_1"]
    assert rows[0][4] == "1000"       # AUC 2500.75 saturates at 1000
    assert rows[0][6] == "2500.75"    # signalValue passes the AUC through
    assert rows[0][5] == "."
    assert rows[0][7:10] == ["0", "0", "0"]
    assert rows[1][4] == "0"          # AUC 0.4 rounds to 0
    assert rows[2][4] == "999"        # round(999.4) = 999


def test_convert_rejects_malformed(tmp_path):
    src = tmp_path / "bad.bed"
    dst = tmp_path / "out.narrowPeak"
    src.write_text("chr1\t1000\n", encoding="utf-8")
    with pytest.raises(ValueError):
        seacr_to_narrowpeak.convert(str(src), str(dst), "g")
    src.write_text("chr1\tx\t2000\t5.0\n", encoding="utf-8")
    with pytest.raises(ValueError):
        seacr_to_narrowpeak.convert(str(src), str(dst), "g")
    src.write_text("chr1\t2000\t1000\t5.0\n", encoding="utf-8")
    with pytest.raises(ValueError):
        seacr_to_narrowpeak.convert(str(src), str(dst), "g")


def test_main_arg_validation(tmp_path, monkeypatch):
    assert seacr_to_narrowpeak.main(["seacr_to_narrowpeak.py", "a", "b"]) == 2
    assert seacr_to_narrowpeak.main(
        ["seacr_to_narrowpeak.py", "a", "b", "name with space"]) == 2
    monkeypatch.setattr(sys, "argv", ["seacr_to_narrowpeak.py"])
    missing = tmp_path / "missing.bed"
    assert seacr_to_narrowpeak.main(
        ["seacr_to_narrowpeak.py", str(missing), str(tmp_path / "o"), "g"]) == 1
