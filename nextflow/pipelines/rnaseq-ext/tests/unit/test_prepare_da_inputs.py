"""Pin the prepare_da_inputs.py semantics: metadata subsetting/ordering,
contrast derivation with the CONTROL group as reference (legacy direction),
blocking factors, and the loud failure modes."""
import pytest

import prepare_da_inputs as pda


def write(path, text):
    path.write_text(text)
    return str(path)


def run(tmp_path, monkeypatch, metadata, matrix, route="salmon",
        control="control", blocking="", extra_contrasts=()):
    argv = [
        "prepare_da_inputs.py",
        "--metadata", write(tmp_path / "meta.csv", metadata),
        "--matrix", write(tmp_path / "matrix.tsv", matrix),
        "--route", route,
        "--control-group", control,
        "--outdir", str(tmp_path / "out"),
    ]
    if blocking:
        argv += ["--blocking", blocking]
    for c in extra_contrasts:
        argv += ["--contrast", c]
    monkeypatch.setattr("sys.argv", argv)
    pda.main()
    out = tmp_path / "out"
    return (out / "matrix.tsv").read_text(), (out / "samplesheet.csv").read_text(), \
        (out / "contrasts.yaml").read_text()


def test_salmon_route_subsets_and_orders_samples(tmp_path, monkeypatch):
    matrix, sheet, contrasts = run(
        tmp_path, monkeypatch,
        metadata="sample,group,batch\nc1,control,b1\nt1,drugA,b1\nt2,drugA,b2\n",
        matrix="gene_id\tgene_name\tjunk\ts1\tc1\tt1\tt2\ng1\tg1\tx\t9\t10\t20\t30\n",
    )
    # junk column dropped; columns ordered to the metadata sheet
    assert matrix.split("\n")[0] == "gene_id\tc1\tt1\tt2"
    assert matrix.split("\n")[1] == "g1\t10\t20\t30"
    # samplesheet keeps sample + variable; batch joins only with --blocking
    assert sheet.split("\n")[0] == "sample,group"
    assert "c1,control" in sheet


def test_contrasts_derive_treat_vs_control(tmp_path, monkeypatch):
    _, _, contrasts = run(
        tmp_path, monkeypatch,
        metadata="sample,group\nc1,control\nt1,drugA\nt2,drugB\n",
        matrix="gene_id\ts1\tc1\tt1\tt2\ng1\t1\t2\t3\t4\n",
    )
    # legacy direction: reference = control, target = the treatment group
    assert "- id: drugA_vs_control" in contrasts
    assert "comparison: ['group', 'control', 'drugA']" in contrasts
    assert "- id: drugB_vs_control" in contrasts
    assert "comparison: ['group', 'control', 'drugB']" in contrasts


def test_blocking_factors_emitted_per_contrast(tmp_path, monkeypatch):
    _, _, contrasts = run(
        tmp_path, monkeypatch,
        metadata="sample,group,batch\nc1,control,b1\nt1,drugA,b1\n",
        matrix="gene_id\tc1\tt1\ng1\t1\t2\n",
        blocking="batch",
    )
    assert "blocking_factors: ['batch']" in contrasts


def test_explicit_contrast_appended(tmp_path, monkeypatch):
    _, _, contrasts = run(
        tmp_path, monkeypatch,
        metadata="sample,group\nc1,control\nt1,drugA\n",
        matrix="gene_id\tc1\tt1\ng1\t1\t2\n",
        extra_contrasts=("drugA_vs_control_custom:control:drugA",),
    )
    assert "- id: drugA_vs_control_custom" in contrasts
    assert "comparison: ['group', 'control', 'drugA']" in contrasts


def test_matrix_missing_sample_fails(tmp_path, monkeypatch):
    with pytest.raises(SystemExit):
        run(
            tmp_path, monkeypatch,
            metadata="sample,group\nc1,control\nzz,drugA\n",
            matrix="gene_id\tc1\tt1\ng1\t1\t2\n",
        )


def test_duplicate_metadata_ids_fail(tmp_path, monkeypatch):
    with pytest.raises(SystemExit):
        run(
            tmp_path, monkeypatch,
            metadata="sample,group\nc1,control\nc1,drugA\n",
            matrix="gene_id\tc1\ng1\t1\n",
        )


def test_no_control_samples_fails(tmp_path, monkeypatch):
    with pytest.raises(SystemExit):
        run(
            tmp_path, monkeypatch,
            metadata="sample,group\nt1,drugA\n",
            matrix="gene_id\tt1\ng1\t1\n",
        )
