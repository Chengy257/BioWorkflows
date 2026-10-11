"""Unit tests for qc_gates_summary.py: the longest-prefix --samples join
(real upstream file names carry infixes such as .mLb.mkD.sorted between the
sample id and the extension), the hard error on unmatched files, flagstat
parsing across samtools formats, and the gate logic."""
import sys

import qc_gates_summary

FLAGSTAT = ("1000 + 0 in total (QC-passed reads + QC-failed reads)\n"
            "900 + 0 mapped (90.00% : N/A)\n")
FLAGSTAT_ALT = ("1000 + 0 in total (QC-passed reads + QC-failed reads)\n"
                "900 + 0 mapped\n")  # samtools <= 1.17 style
DUP = ("## METRICS CLASS\tpicard.sam.DuplicationMetrics\n"
       "LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\t"
       "SECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\t"
       "UNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\t"
       "READ_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\n"
       "s1\t100\t400\t0\t100\t20\t80\t0\t0.2000\n")


def run(tmp_path, monkeypatch, flagstat_names, dup_names, samples, extra=()):
    flag_dir = tmp_path / "flagstat"
    dup_dir = tmp_path / "dup"
    flag_dir.mkdir()
    dup_dir.mkdir()
    for name in flagstat_names:
        (flag_dir / name).write_text(FLAGSTAT, encoding="utf-8")
    for name in dup_names:
        (dup_dir / name).write_text(DUP, encoding="utf-8")
    out = tmp_path / "qc_gates.tsv"
    mqc = tmp_path / "qc_gates_mqc.tsv"
    argv = ["qc_gates_summary.py", "--mapping-rate-min", "0.70",
            "--dup-rate-max", "0.50", "--organelle-max", "",
            "--dup-suffix", ".metrics.txt",
            "--flagstat-suffix", ".mLb.clN.sorted.bam.flagstat",
            "--samples", samples, "--dup-metrics", f"{dup_dir}/*",
            "--out", str(out), "--mqc", str(mqc)] + list(extra) + \
        [str(flag_dir / n) for n in flagstat_names]
    monkeypatch.setattr(sys, "argv", argv)
    rc = qc_gates_summary.main()
    if rc:
        return rc, None, None
    lines = out.read_text().splitlines()
    header = lines[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in lines[1:]]
    mqc_lines = mqc.read_text().splitlines()
    return 0, rows, mqc_lines


def test_longest_prefix_join_real_names(tmp_path, monkeypatch):
    # real upstream names: flagstat {id}.mLb.clN.sorted.bam.flagstat,
    # picard {id}.mLb.mkD.sorted.metrics.txt - both resolve to the bare id
    rc, rows, _ = run(tmp_path, monkeypatch,
                      ["s1.mLb.clN.sorted.bam.flagstat",
                       "s11.mLb.clN.sorted.bam.flagstat"],
                      ["s1.mLb.mkD.sorted.metrics.txt",
                       "s11.mLb.mkD.sorted.metrics.txt"],
                      samples="s1,s11")
    assert rc == 0
    by_sample = {r["sample"]: r for r in rows}
    assert set(by_sample) == {"s1", "s11"}
    # s11 must not be swallowed by the s1 prefix
    assert by_sample["s11"]["dup_rate"] == "0.2000"
    assert by_sample["s1"]["dup_rate"] == "0.2000"
    assert all(r["mapping_rate"] == "0.9000" for r in rows)
    assert all(r["gate"] == "PASS" for r in rows)


def test_unmatched_file_hard_error(tmp_path, monkeypatch, capsys):
    rc, _, _ = run(tmp_path, monkeypatch,
                   ["stray.mLb.clN.sorted.bam.flagstat"],
                   ["s1.mLb.mkD.sorted.metrics.txt"],
                   samples="s1")
    assert rc == 1
    assert "stray.mLb.clN.sorted.bam.flagstat" in capsys.readouterr().err


def test_no_samples_falls_back_to_suffix_strip(tmp_path, monkeypatch):
    rc, rows, _ = run(tmp_path, monkeypatch,
                      ["s1.mLb.clN.sorted.bam.flagstat"], [], samples="")
    assert rc == 0
    assert rows[0]["sample"] == "s1"
    assert rows[0]["dup_rate"] == "NA"


def test_gate_thresholds(tmp_path, monkeypatch):
    flag_dir = tmp_path / "flagstat"
    flag_dir.mkdir()
    (flag_dir / "s1.mLb.clN.sorted.bam.flagstat").write_text(
        FLAGSTAT_ALT, encoding="utf-8")
    out = tmp_path / "qc_gates.tsv"
    argv = ["qc_gates_summary.py", "--mapping-rate-min", "0.95",
            "--dup-rate-max", "0.50", "--organelle-max", "",
            "--flagstat-suffix", ".mLb.clN.sorted.bam.flagstat",
            "--samples", "s1", "--out", str(out), "--mqc", str(out) + ".mqc",
            str(flag_dir / "s1.mLb.clN.sorted.bam.flagstat")]
    monkeypatch.setattr(sys, "argv", argv)
    assert qc_gates_summary.main() == 0
    lines = [line.split("\t") for line in out.read_text().splitlines()]
    rows = {row[0]: row for row in lines[1:]}
    assert rows["s1"] == ["s1", "0.9000", "NA", "NA", "FAIL"]


def test_mqc_header(tmp_path, monkeypatch):
    rc, _, mqc_lines = run(tmp_path, monkeypatch,
                           ["s1.mLb.clN.sorted.bam.flagstat"], [],
                           samples="s1")
    assert rc == 0
    assert mqc_lines[0] == "# id: 'qc_gates_table'"
    assert mqc_lines[3] == "# plot_type: 'table'"
