"""Unit tests for diffbind_sheet.py: the authoritative control exclusion
(--exclude-ids carries the union of declared control ids - the workflow's
groupControl mapping), per-condition replicate numbering independent of
sheet row order, and the error paths."""
import glob
import sys

import diffbind_sheet  # noqa: E402  (conftest.py puts bin/ on sys.path)

SHEET = """sample_id,group,condition,control,batch
ctl1,INPUT,input,,
s2,H3K27ac_MUT,mut,ctl1,b2
s1,H3K27ac_WT,WT,ctl1,b1
s3,H3K27ac_MUT,mut,ctl1,
s4,H3K27ac_WT,WT,ctl1,
"""
# rows deliberately NOT grouped and with the contrast condition B (mut) listed
# first: replicate numbering and the Fold direction must not depend on it


def write_stage(tmp_path):
    bam_dir = tmp_path / "bams"
    peak_dir = tmp_path / "peaks"
    bam_dir.mkdir()
    peak_dir.mkdir()
    for sid in ("s1", "s2", "s3", "s4", "ctl1"):
        (bam_dir / f"{sid}.mLb.clN.sorted.bam").write_bytes(b"")
        (peak_dir / f"{sid}_peaks.narrowPeak").write_text("")
    return bam_dir, peak_dir


def run_sheet(monkeypatch, tmp_path, sheet_text, exclude_ids):
    sheet = tmp_path / "sheet.csv"
    sheet.write_text(sheet_text, encoding="utf-8")
    bam_dir, peak_dir = write_stage(tmp_path)
    out = tmp_path / "out.tsv"
    # the pipeline stages expanded file lists (the shell expands bams/* in a
    # real task); pytest passes no shell, so expand here
    staged = sorted(glob.glob(f"{bam_dir}/*")) + sorted(glob.glob(f"{peak_dir}/*"))
    argv = ["diffbind_sheet.py", "--contrast", "WT_vs_mut", "--sheet", str(sheet),
            "--bam-suffix", ".mLb.clN.sorted.bam", "--peak-suffix",
            "_peaks.narrowPeak", "--peak-type", "narrowpeak",
            "--use-controls", "F", "--exclude-ids", exclude_ids,
            "--out", str(out)] + staged
    monkeypatch.setattr(sys, "argv", argv)
    diffbind_sheet.main()
    lines = out.read_text().splitlines()
    header = lines[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in lines[1:]]
    return rows


def test_exclude_ids_removes_shared_control(tmp_path, monkeypatch):
    # ctl1 sits in its own group AND carries a condition that is not part of
    # the contrast, but a sheet could easily give a control row the contrast
    # condition: the authoritative exclusion set must win over row heuristics
    sheet = SHEET.replace("ctl1,INPUT,input,,", "ctl1,INPUT,WT,,")
    rows = run_sheet(monkeypatch, tmp_path, sheet, exclude_ids="ctl1")
    # output rows follow the sheet's own row order (s2 first by design)
    assert [r["SampleID"] for r in rows] == ["s2", "s1", "s3", "s4"]
    assert all(r["SampleID"] != "ctl1" for r in rows)


def test_replicate_numbering_per_condition(tmp_path, monkeypatch):
    # s4 (WT) precedes s1 in the sheet; mut rows come first entirely
    rows = run_sheet(monkeypatch, tmp_path, SHEET, exclude_ids="ctl1")
    numbering = {(r["Condition"], r["SampleID"]): r["Replicate"] for r in rows}
    assert numbering == {("WT", "s1"): "1", ("WT", "s4"): "2",
                         ("mut", "s2"): "1", ("mut", "s3"): "2"}


def test_use_controls_attaches_control_bam(tmp_path, monkeypatch):
    sheet = tmp_path / "sheet.csv"
    sheet.write_text(SHEET, encoding="utf-8")
    bam_dir, peak_dir = write_stage(tmp_path)
    out = tmp_path / "out.tsv"
    staged = sorted(glob.glob(f"{bam_dir}/*")) + sorted(glob.glob(f"{peak_dir}/*"))
    argv = ["diffbind_sheet.py", "--contrast", "WT_vs_mut", "--sheet", str(sheet),
            "--bam-suffix", ".mLb.clN.sorted.bam", "--peak-suffix",
            "_peaks.narrowPeak", "--peak-type", "narrowpeak",
            "--use-controls", "T", "--exclude-ids", "ctl1",
            "--out", str(out)] + staged
    monkeypatch.setattr(sys, "argv", argv)
    diffbind_sheet.main()
    rows = [line.split("\t") for line in out.read_text().splitlines()[1:]]
    # columns: SampleID, Condition, Replicate, bamReads, bamControl, ...
    control_cells = {row[0]: row[4] for row in rows}
    assert all("ctl1.mLb.clN.sorted.bam" in c for c in control_cells.values())


def test_missing_staged_file_fails(tmp_path, monkeypatch, capsys):
    sheet = tmp_path / "sheet.csv"
    sheet.write_text(SHEET, encoding="utf-8")
    bam_dir, peak_dir = write_stage(tmp_path)
    (peak_dir / "s4_peaks.narrowPeak").unlink()
    out = tmp_path / "out.tsv"
    staged = sorted(glob.glob(f"{bam_dir}/*")) + sorted(glob.glob(f"{peak_dir}/*"))
    argv = ["diffbind_sheet.py", "--contrast", "WT_vs_mut", "--sheet", str(sheet),
            "--bam-suffix", ".mLb.clN.sorted.bam", "--peak-suffix",
            "_peaks.narrowPeak", "--peak-type", "narrowpeak",
            "--use-controls", "F", "--exclude-ids", "ctl1",
            "--out", str(out)] + staged
    monkeypatch.setattr(sys, "argv", argv)
    assert diffbind_sheet.main() == 1
    assert "s4" in capsys.readouterr().err


def test_no_matching_conditions_fails(tmp_path, monkeypatch, capsys):
    sheet = tmp_path / "sheet.csv"
    sheet.write_text(SHEET, encoding="utf-8")
    bam_dir, peak_dir = write_stage(tmp_path)
    out = tmp_path / "out.tsv"
    argv = ["diffbind_sheet.py", "--contrast", "A_vs_B", "--sheet", str(sheet),
            "--bam-suffix", ".mLb.clN.sorted.bam", "--peak-suffix",
            "_peaks.narrowPeak", "--peak-type", "narrowpeak",
            "--use-controls", "F", "--exclude-ids", "ctl1",
            "--out", str(out), f"{bam_dir}/*", f"{peak_dir}/*"]
    monkeypatch.setattr(sys, "argv", argv)
    assert diffbind_sheet.main() == 1
    assert "no samples match" in capsys.readouterr().err
