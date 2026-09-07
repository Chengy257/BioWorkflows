"""Unit tests for the pure-stdlib scripts in workflow/scripts/ (docs/TODO.md
item 3).

Each script module is loaded straight from its file via importlib (they are
standalone scripts, not a package) and main() is invoked with a patched
sys.argv; outputs are asserted under pytest's tmp_path.

Run: python -m pytest tests -q  (paths derive from __file__, so the working
directory does not matter).
"""
import importlib.util
import os
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(REPO, "workflow", "scripts")

_MODULES = {}


def _load(name):
    """Import one workflow script by file name (without .py)."""
    if name not in _MODULES:
        spec = importlib.util.spec_from_file_location(
            f"srna_{name}", os.path.join(SCRIPTS, f"{name}.py"))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        _MODULES[name] = module
    return _MODULES[name]


def _run(monkeypatch, script, argv):
    """Run script.main() with sys.argv patched to the given arguments."""
    monkeypatch.setattr(sys, "argv", [f"{script}.py"] + argv)
    _load(script).main()


def _write_counts(path, rows):
    """Write one per-sample count file: header + (feature, count) rows."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("feature\treads\n")
        for feature, count in rows:
            fh.write(f"{feature}\t{count}\n")


# ---------------------------------------------------------------------
# count_features.py
# ---------------------------------------------------------------------
SAM_LINE = "r{n}\t{flag}\t{rname}\t100\t255\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\tIIIII"


def _sam(*rows):
    return "@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:1000\n" + \
        "\n".join(SAM_LINE.format(**row) for row in rows) + "\n"


def test_count_features_skips_header_secondary_unmapped_and_sorts(tmp_path, monkeypatch):
    sam = tmp_path / "sample.sam"
    sam.write_text(_sam(
        {"n": 1, "flag": 0, "rname": "chr1"},     # counted
        {"n": 2, "flag": 256, "rname": "chr1"},   # secondary (0x100): skipped
        {"n": 3, "flag": 0, "rname": "chr2"},     # counted
        {"n": 4, "flag": 4, "rname": "*"},        # unmapped (rname *): skipped
        {"n": 5, "flag": 0, "rname": "chr1"},     # counted
        {"n": 6, "flag": 0, "rname": "chr1"},     # counted
    ), encoding="utf-8")
    out = tmp_path / "counts.tsv"
    _run(monkeypatch, "count_features", [str(sam), "-o", str(out)])
    assert out.read_text(encoding="utf-8") == (
        "feature\treads\n"
        "chr1\t3\n"
        "chr2\t1\n"
        "__mapped_total\t4\n"
    )


def test_count_features_sorts_features_lexicographically(tmp_path, monkeypatch):
    sam = tmp_path / "sample.sam"
    sam.write_text(_sam(
        {"n": 1, "flag": 0, "rname": "chr2"},
        {"n": 2, "flag": 0, "rname": "chr10"},
        {"n": 3, "flag": 0, "rname": "chrB"},
        {"n": 4, "flag": 0, "rname": "chrA"},
    ), encoding="utf-8")
    out = tmp_path / "counts.tsv"
    _run(monkeypatch, "count_features", [str(sam), "-o", str(out)])
    assert out.read_text(encoding="utf-8") == (
        "feature\treads\n"
        "chr10\t1\n"
        "chr2\t1\n"
        "chrA\t1\n"
        "chrB\t1\n"
        "__mapped_total\t4\n"
    )


def test_count_features_writes_only_the_sentinel_when_nothing_maps(tmp_path, monkeypatch):
    sam = tmp_path / "sample.sam"
    sam.write_text(
        "@SQ\tSN:chr1\tLN:1000\n"
        "r1\t256\tchr1\t100\t255\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\tIIIII\n"
        "r2\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTACGTACGTACGT\tIIIII\n",
        encoding="utf-8")
    out = tmp_path / "counts.tsv"
    _run(monkeypatch, "count_features", [str(sam), "-o", str(out)])
    assert out.read_text(encoding="utf-8") == "feature\treads\n__mapped_total\t0\n"


# ---------------------------------------------------------------------
# merge_counts.py
# ---------------------------------------------------------------------
def test_read_counts_skips_header_and_parses_integers(tmp_path):
    path = tmp_path / "s1_counts.txt"
    _write_counts(path, [("b", 2), ("a", 1), ("__mapped_total", 3)])
    assert _load("merge_counts").read_counts(str(path)) == {
        "b": 2, "a": 1, "__mapped_total": 3}


def test_merge_counts_writes_matrix_rpm_and_all_classes_tables(tmp_path, monkeypatch):
    indir = tmp_path / "4.expression"
    _write_counts(indir / "rRNA" / "s1_counts.txt",
                  [("f2", 5), ("f1", 10), ("__mapped_total", 15)])
    _write_counts(indir / "rRNA" / "s2_counts.txt",
                  [("f1", 7), ("f3", 3), ("__mapped_total", 10)])
    _write_counts(indir / "miRNA" / "s1_counts.txt",
                  [("__mapped_total", 0)])          # class mapped nothing for s1
    _write_counts(indir / "miRNA" / "s2_counts.txt",
                  [("m1", 4), ("__mapped_total", 4)])
    _run(monkeypatch, "merge_counts",
         ["--indir", str(indir), "--classes", "rRNA", "miRNA",
          "--samples", "s1", "s2"])

    counts_tsv = (indir / "rRNA" / "rRNA_counts.tsv").read_text(encoding="utf-8")
    assert counts_tsv == (
        "feature\ts1\ts2\n"
        "f1\t10\t7\n"
        "f2\t5\t0\n"
        "f3\t0\t3\n"
    )
    # RPM = count / class mapped total * 1e6, %.3f; 0.000 when the class
    # mapped nothing for that sample (miRNA s1 has __mapped_total 0).
    rpm_tsv = (indir / "rRNA" / "rRNA_RPM.tsv").read_text(encoding="utf-8")
    assert rpm_tsv == (
        "feature\ts1\ts2\n"
        "f1\t666666.667\t700000.000\n"
        "f2\t333333.333\t0.000\n"
        "f3\t0.000\t300000.000\n"
    )
    miRNA_rpm = (indir / "miRNA" / "miRNA_RPM.tsv").read_text(encoding="utf-8")
    assert miRNA_rpm == (
        "feature\ts1\ts2\n"
        "m1\t0.000\t1000000.000\n"
    )
    # __mapped_total never appears as a feature row; classes keep their order.
    all_tsv = (indir / "all_classes_counts.tsv").read_text(encoding="utf-8")
    assert all_tsv == (
        "class\tfeature\ts1\ts2\n"
        "rRNA\tf1\t10\t7\n"
        "rRNA\tf2\t5\t0\n"
        "rRNA\tf3\t0\t3\n"
        "miRNA\tm1\t0\t4\n"
    )


def test_merge_counts_empty_class_list_writes_header_only_summary(tmp_path, monkeypatch):
    indir = tmp_path / "4.expression"
    indir.mkdir(parents=True)
    _run(monkeypatch, "merge_counts",
         ["--indir", str(indir), "--classes", "--samples", "s1"])
    assert (indir / "all_classes_counts.tsv").read_text(encoding="utf-8") == \
        "class\tfeature\ts1\n"


# ---------------------------------------------------------------------
# cascade_summary.py helpers
# ---------------------------------------------------------------------
def test_fq_read_count_counts_four_line_records(tmp_path):
    fq = tmp_path / "s1_unmapped.fq"
    with open(fq, "w", encoding="utf-8", newline="\n") as fh:
        for i in range(3):
            fh.write(f"@r{i}\nACGT\n+\nIIII\n")
    assert _load("cascade_summary").fq_read_count(str(fq)) == 3


def test_trimmed_reads_parses_trim_galore_report(tmp_path):
    report = tmp_path / "s1_trimming_report.txt"
    report.write_text(
        "=== Summary ===\n"
        "Total reads processed: 2,000\n"
        "Reads with adapters: 500 (25.0%)\n"
        "Reads written (passing filters): 1,850 (92.5%)\n",
        encoding="utf-8")
    assert _load("cascade_summary").trimmed_reads(str(report)) == 1850


def test_trimmed_reads_returns_zero_without_matching_line(tmp_path):
    report = tmp_path / "s1_trimming_report.txt"
    report.write_text("=== Summary ===\nTotal reads processed: 10\n", encoding="utf-8")
    assert _load("cascade_summary").trimmed_reads(str(report)) == 0


def test_mapped_total_reads_the_sentinel_row(tmp_path):
    counts = tmp_path / "s1_counts.txt"
    _write_counts(counts, [("f1", 7), ("__mapped_total", 9)])
    assert _load("cascade_summary").mapped_total(str(counts)) == 9


def test_mapped_total_defaults_to_zero_without_sentinel(tmp_path):
    counts = tmp_path / "s1_counts.txt"
    _write_counts(counts, [("f1", 7)])
    assert _load("cascade_summary").mapped_total(str(counts)) == 0


# ---------------------------------------------------------------------
# cascade_summary.py main()
# ---------------------------------------------------------------------
def _write_fastq(path, reads):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for i in range(reads):
            fh.write(f"@r{i}\nACGT\n+\nIIII\n")


def _write_trim_report(path, written):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "=== Summary ===\n"
        f"Total reads processed: {written + 10:,}\n"
        f"Reads written (passing filters): {written:,} (90.0%)\n",
        encoding="utf-8")


def _build_summary_inputs(tmp_path):
    """Assemble 4.expression counts, 3.align unmapped fastqs and trim reports.

    rRNA unmapped fastqs exist (s1: 2 reads, s2: 1 read); miRNA unmapped
    fastqs are missing (-1); the genome unmapped fastq exists only for s1.
    """
    indir = tmp_path / "4.expression"
    mapped = {"s1": {"rRNA": 15, "miRNA": 0}, "s2": {"rRNA": 10, "miRNA": 4}}
    for sample in ("s1", "s2"):
        for klass in ("rRNA", "miRNA"):
            rows = ([("f1", mapped[sample][klass])]
                    if mapped[sample][klass] else [])
            rows.append(("__mapped_total", mapped[sample][klass]))
            _write_counts(indir / klass / f"{sample}_counts.txt", rows)
        _write_fastq(tmp_path / "3.align" / "filter" / "rRNA" /
                     f"{sample}_unmapped.fq", 2 if sample == "s1" else 1)
        _write_trim_report(tmp_path / "logs" / "trim" /
                           f"{sample}_trimming_report.txt",
                           90 if sample == "s1" else 180)
    _write_fastq(tmp_path / "3.align" / "genome" / "s1_unmapped.fq", 3)
    return indir


def test_cascade_summary_wide_and_mqc_outputs_genome_configured(tmp_path, monkeypatch):
    indir = _build_summary_inputs(tmp_path)
    outdir = tmp_path / "5.QC"
    _run(monkeypatch, "cascade_summary",
         ["--indir", str(indir), "--classes", "rRNA", "miRNA",
          "--samples", "s1", "s2",
          "--trim-reports", str(tmp_path / "logs" / "trim" / "s1_trimming_report.txt"),
          str(tmp_path / "logs" / "trim" / "s2_trimming_report.txt"),
          "--raw-counts", "100", "200",
          "--outdir", str(outdir), "--genome-configured", "true"])
    wide = (outdir / "cascade_summary.tsv").read_text(encoding="utf-8")
    assert wide == (
        "sample\traw\ttrimmed\trRNA_mapped\trRNA_unmapped"
        "\tmiRNA_mapped\tmiRNA_unmapped\tgenome_unmapped\n"
        "s1\t100\t90\t15\t2\t0\t-1\t3\n"
        "s2\t200\t180\t10\t1\t4\t-1\t-1\n"
    )
    mqc = (outdir / "cascade_summary_mqc.tsv").read_text(encoding="utf-8")
    assert mqc == (
        "# id: 'srna_cascade'\n"
        "# section_name: 'sRNA cascade read fate'\n"
        "# plot_type: 'bargraph'\n"
        "# pconfig:\n"
        "#     id: 'srna_cascade_bg'\n"
        "#     title: 'sRNA cascade read fate'\n"
        "sample\tstage\treads\n"
        "s1\traw\t100\n"
        "s1\ttrimmed\t90\n"
        "s1\trRNA_mapped\t15\n"
        "s1\trRNA_unmapped\t2\n"
        "s1\tmiRNA_mapped\t0\n"
        "s1\tmiRNA_unmapped\t-1\n"
        "s1\tgenome_unmapped\t3\n"
        "s2\traw\t200\n"
        "s2\ttrimmed\t180\n"
        "s2\trRNA_mapped\t10\n"
        "s2\trRNA_unmapped\t1\n"
        "s2\tmiRNA_mapped\t4\n"
        "s2\tmiRNA_unmapped\t-1\n"
        "s2\tgenome_unmapped\t-1\n"
    )


def test_cascade_summary_genome_not_configured_writes_na(tmp_path, monkeypatch):
    indir = _build_summary_inputs(tmp_path)
    outdir = tmp_path / "5.QC"
    _run(monkeypatch, "cascade_summary",
         ["--indir", str(indir), "--classes", "rRNA",
          "--samples", "s1",
          "--trim-reports", str(tmp_path / "logs" / "trim" / "s1_trimming_report.txt"),
          "--raw-counts", "100",
          "--outdir", str(outdir), "--genome-configured", "false"])
    wide = (outdir / "cascade_summary.tsv").read_text(encoding="utf-8")
    assert wide == (
        "sample\traw\ttrimmed\trRNA_mapped\trRNA_unmapped\tgenome_unmapped\n"
        "s1\t100\t90\t15\t2\tNA\n"
    )
    mqc = (outdir / "cascade_summary_mqc.tsv").read_text(encoding="utf-8")
    assert mqc.endswith(
        "s1\trRNA_mapped\t15\n"
        "s1\trRNA_unmapped\t2\n"
        "s1\tgenome_unmapped\tNA\n"
    )


def test_cascade_summary_empty_cascade_still_writes_zero_class_summary(tmp_path, monkeypatch):
    outdir = tmp_path / "5.QC"
    report = tmp_path / "logs" / "trim" / "s1_trimming_report.txt"
    _write_trim_report(report, 90)
    _run(monkeypatch, "cascade_summary",
         ["--indir", str(tmp_path / "4.expression"), "--classes",
          "--samples", "s1",
          "--trim-reports", str(report),
          "--raw-counts", "100",
          "--outdir", str(outdir), "--genome-configured", "false"])
    assert (outdir / "cascade_summary.tsv").read_text(encoding="utf-8") == (
        "sample\traw\ttrimmed\tgenome_unmapped\n"
        "s1\t100\t90\tNA\n"
    )
