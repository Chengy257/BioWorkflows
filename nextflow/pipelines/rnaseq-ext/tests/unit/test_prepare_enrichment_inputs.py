"""Pin the prepare_enrichment_inputs.py semantics: the legacy DEG split
(padj<=max AND |log2FC|>=min), the ranked list ordered by the stat column,
and the loud missing-column failure."""
import pytest

import prepare_enrichment_inputs as pei


def run(tmp_path, monkeypatch, table, **kw):
    argv = [
        "prepare_enrichment_inputs.py",
        str(tmp_path / "de.tsv"), str(tmp_path / "out"), "drugA_vs_control",
        "--id-col", kw.get("id_col", "gene_id"),
        "--log2fc-col", kw.get("log2fc_col", "log2FoldChange"),
        "--padj-col", kw.get("padj_col", "padj"),
        "--stat-col", kw.get("stat_col", "stat"),
        "--padj-max", str(kw.get("padj_max", 0.05)),
        "--log2fc-min", str(kw.get("log2fc_min", 1.0)),
    ]
    (tmp_path / "de.tsv").write_text(table)
    monkeypatch.setattr("sys.argv", argv)
    pei.main()
    out = tmp_path / "out"
    read = lambda name: (out / name).read_text()
    return read("drugA_vs_control_UP.DEGs.txt"), read("drugA_vs_control_DOWN.DEGs.txt"), \
        read("drugA_vs_control_ALL.DEGs.txt"), read("drugA_vs_control_ranked.tsv")


DE_TABLE = """gene_id\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj
g_up1\t2.5\t0.5\t5.0\t0.0001\t0.001
g_up2\t1.0\t0.5\t2.0\t0.01\t0.04
g_edge\t0.9\t0.5\t1.8\t0.01\t0.04
g_down\t-2.0\t0.5\t-4.0\t0.0001\t0.003
g_ns\t0.1\t0.5\t0.2\t0.9\t0.9
g_nonsig_padj\t5.0\t0.5\t10.0\t0.01\t0.20
"""


def test_legacy_split_thresholds(tmp_path, monkeypatch):
    up, down, all_, ranked = run(tmp_path, monkeypatch, DE_TABLE)
    assert up.split() == ["g_up1", "g_up2"]      # |log2FC| >= 1.0 and padj <= 0.05
    assert "g_edge" not in up.split()            # 0.9 < 1.0 cutoff
    assert "g_nonsig_padj" not in up.split()     # padj 0.20 > 0.05
    assert down.split() == ["g_down"]
    assert all_.split() == ["g_up1", "g_up2", "g_down"]


def test_ranked_ordered_by_stat_desc(tmp_path, monkeypatch):
    _, _, _, ranked = run(tmp_path, monkeypatch, DE_TABLE)
    rows = [l.split("\t") for l in ranked.strip().split("\n")[1:]]
    stats = [float(r[1]) for r in rows]
    assert stats == sorted(stats, reverse=True)
    assert rows[0][0] == "g_nonsig_padj"         # stat 10.0 ranks first


def test_missing_column_fails(tmp_path, monkeypatch):
    with pytest.raises(SystemExit):
        run(tmp_path, monkeypatch, "gene_id\tlog2FC\tpadj\ng1\t1\t0.01\n",
            log2fc_col="log2FoldChange")
