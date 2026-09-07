"""Unit tests for the seclip-seq peak-annotation helper scripts.

Covers workflow/scripts/gtf_to_gene_regions.py (GTF attribute parsing with
quoted/unquoted values, gene/exon interval extraction, output contracts,
error paths) and workflow/scripts/annotate_peaks.py (the intersect/closest
merge logic: exon vs gene vs intergenic classification, score-column
selection for BED6 and consensus BEDs, nearest-gene table lookup).

Standard library + pytest only; all paths derive from __file__, so the suite
is independent of the working directory.
"""
import importlib.util
import os

import pytest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.join(os.path.dirname(TESTS_DIR), "workflow", "scripts")


def _load_module(name):
    path = os.path.join(SCRIPTS_DIR, name + ".py")
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


gtf = _load_module("gtf_to_gene_regions")
ann = _load_module("annotate_peaks")


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


# ---------------------------------------------------------------------------
# gtf_to_gene_regions.parse_attributes
# ---------------------------------------------------------------------------

def test_parse_attributes_quoted():
    attrs = gtf.parse_attributes(
        'gene_id "ENSG1"; gene_name "Alpha one"; gene_biotype "protein_coding";')
    assert attrs == {"gene_id": "ENSG1", "gene_name": "Alpha one",
                     "gene_biotype": "protein_coding"}


def test_parse_attributes_unquoted():
    attrs = gtf.parse_attributes("gene_id G1; transcript_id T1.2")
    assert attrs == {"gene_id": "G1", "transcript_id": "T1.2"}


def test_parse_attributes_empty_and_valueless():
    assert gtf.parse_attributes("") == {}
    assert gtf.parse_attributes('gene_id "X";; extra') == {"gene_id": "X", "extra": ""}


# ---------------------------------------------------------------------------
# gtf_to_gene_regions.parse_gtf
# ---------------------------------------------------------------------------

def test_parse_gtf_gene_and_exon_intervals(tmp_path):
    gtf_text = (
        'chr1\tsim\tgene\t101\t500\t.\t+\t.\tgene_id "g1"; gene_name "Alpha"; gene_biotype "protein_coding";\n'
        'chr1\tsim\texon\t101\t200\t.\t+\t.\tgene_id "g1";\n'
        'chr1\tsim\texon\t301\t500\t.\t+\t.\tgene_id "g1";\n'
    )
    genes, exons = gtf.parse_gtf(_write(tmp_path, "a.gtf", gtf_text))
    gene = genes["g1"]
    # 1-based closed GTF -> 0-based half-open BED
    assert (gene["chrom"], gene["start"], gene["end"], gene["strand"]) == ("chr1", 100, 500, "+")
    assert gene["name"] == "Alpha" and gene["biotype"] == "protein_coding"
    assert sorted((e["start"], e["end"]) for e in exons) == [(100, 200), (300, 500)]


def test_parse_gtf_gene_span_from_exons_without_gene_row(tmp_path):
    gtf_text = (
        'chr2\tsim\texon\t50\t100\t.\t-\t.\tgene_id "g2"; gene_biotype "miRNA";\n'
        'chr2\tsim\texon\t300\t400\t.\t-\t.\tgene_id "g2";\n'
    )
    genes, exons = gtf.parse_gtf(_write(tmp_path, "b.gtf", gtf_text))
    gene = genes["g2"]
    assert (gene["start"], gene["end"], gene["strand"]) == (49, 400, "-")
    # attributes fall back to the exon records, then to gene_id / '.'
    assert gene["name"] == "g2" and gene["biotype"] == "miRNA"
    assert len(exons) == 2


def test_parse_gtf_ignores_transcripts_comments_and_blanks(tmp_path):
    gtf_text = (
        "# a comment line\n"
        "\n"
        'chr1\tsim\ttranscript\t1\t900\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        'chr1\tsim\tCDS\t10\t50\t.\t+\t.\tgene_id "g1";\n'
        'chr1\tsim\texon\t1\t100\t.\t+\t.\tgene_id "g1";\n'
    )
    genes, exons = gtf.parse_gtf(_write(tmp_path, "c.gtf", gtf_text))
    assert list(genes) == ["g1"]
    assert (genes["g1"]["start"], genes["g1"]["end"]) == (0, 100)
    assert len(exons) == 1


def test_parse_gtf_gene_row_attributes_win_and_spans_merge(tmp_path):
    # exon-only span first, then a narrower gene row: the span must stay the
    # union and the gene-row name/biotype must override the exon fallbacks.
    gtf_text = (
        'chr1\tsim\texon\t50\t100\t.\t+\t.\tgene_id "g3"; gene_name "ExonName"; gene_biotype "x";\n'
        'chr1\tsim\tgene\t60\t90\t.\t+\t.\tgene_id "g3"; gene_name "GeneName"; gene_biotype "y";\n'
    )
    genes, _exons = gtf.parse_gtf(_write(tmp_path, "d.gtf", gtf_text))
    gene = genes["g3"]
    assert (gene["start"], gene["end"]) == (49, 100)
    assert gene["name"] == "GeneName" and gene["biotype"] == "y"


def test_parse_gtf_duplicate_exons_deduplicated(tmp_path):
    exon = 'chr1\tsim\texon\t10\t20\t.\t+\t.\tgene_id "g4";\n'
    genes, exons = gtf.parse_gtf(_write(tmp_path, "e.gtf", exon + exon))
    assert list(genes) == ["g4"]
    assert len(exons) == 1


def test_parse_gtf_missing_gene_id_on_gene_row_raises(tmp_path):
    with pytest.raises(gtf.GtfParseError):
        gtf.parse_gtf(_write(tmp_path, "f.gtf", 'chr1\tsim\tgene\t1\t9\t.\t+\t.\tgene_name "x";\n'))


def test_parse_gtf_missing_gene_id_on_exon_row_raises(tmp_path):
    with pytest.raises(gtf.GtfParseError):
        gtf.parse_gtf(_write(tmp_path, "g.gtf", 'chr1\tsim\texon\t1\t9\t.\t+\t.\ttranscript_id "t";\n'))


def test_parse_gtf_short_line_raises(tmp_path):
    with pytest.raises(gtf.GtfParseError):
        gtf.parse_gtf(_write(tmp_path, "h.gtf", "chr1\tsim\tgene\t1\t9\n"))


def test_parse_gtf_non_integer_coordinates_raise(tmp_path):
    line = 'chr1\tsim\tgene\tabc\t9\t.\t+\t.\tgene_id "g5";\n'
    with pytest.raises(gtf.GtfParseError):
        gtf.parse_gtf(_write(tmp_path, "i.gtf", line))


def test_parse_gtf_end_before_start_raises(tmp_path):
    line = 'chr1\tsim\tgene\t50\t9\t.\t+\t.\tgene_id "g6";\n'
    with pytest.raises(gtf.GtfParseError):
        gtf.parse_gtf(_write(tmp_path, "j.gtf", line))


def test_parse_gtf_chrom_conflict_for_same_gene_id_raises(tmp_path):
    gtf_text = (
        'chr1\tsim\tgene\t1\t9\t.\t+\t.\tgene_id "g7";\n'
        'chr2\tsim\texon\t1\t9\t.\t+\t.\tgene_id "g7";\n'
    )
    with pytest.raises(gtf.GtfParseError):
        gtf.parse_gtf(_write(tmp_path, "k.gtf", gtf_text))


def test_write_outputs_contract(tmp_path):
    gtf_text = (
        'chr1\tsim\tgene\t101\t500\t.\t+\t.\tgene_id "gA"; gene_name "Alpha";\n'
        'chr1\tsim\texon\t101\t200\t.\t+\t.\tgene_id "gA";\n'
        'chr2\tsim\tgene\t1\t100\t.\t-\t.\tgene_id "gB"; gene_biotype "miRNA";\n'
    )
    genes, exons = gtf.parse_gtf(_write(tmp_path, "l.gtf", gtf_text))
    genes_out = str(tmp_path / "genes.bed")
    exons_out = str(tmp_path / "exons.bed")
    table_out = str(tmp_path / "genes.tsv")
    gtf.write_outputs(genes, exons, genes_out, exons_out, table_out)

    gene_lines = open(genes_out).read().splitlines()
    assert gene_lines == [
        "chr1\t100\t500\tgA\t0\t+",
        "chr2\t0\t100\tgB\t0\t-",
    ]
    exon_lines = open(exons_out).read().splitlines()
    assert exon_lines == ["chr1\t100\t200\tgA\t0\t+"]
    table_lines = open(table_out).read().splitlines()
    assert table_lines == [
        "gene_id\tgene_name\tgene_biotype",
        "gA\tAlpha\t.",
        "gB\tgB\tmiRNA",
    ]


def test_gtf_main_end_to_end(tmp_path):
    gtf_text = 'chr1\tsim\tgene\t101\t500\t.\t+\t.\tgene_id "gZ"; gene_name "Zeta"; gene_biotype "lncRNA";\n'
    gtf_file = _write(tmp_path, "m.gtf", gtf_text)
    outs = [str(tmp_path / n) for n in ("g.bed", "e.bed", "t.tsv")]
    assert gtf.main(["--gtf", gtf_file, "--genes-out", outs[0],
                     "--exons-out", outs[1], "--table-out", outs[2]]) == 0
    assert all(os.path.exists(p) for p in outs)


# ---------------------------------------------------------------------------
# annotate_peaks
# ---------------------------------------------------------------------------

PEAKS_BED6 = (
    "chr1\t100\t200\tsite1\t42\t+\n"     # exon-overlapping
    "chr1\t500\t600\tsite2\t7\t-\n"      # intronic (gene body, no exon)
    "chr1\t5000\t5100\tsite3\t3\t+\n"    # intergenic, nearest gene left
    "chrX\t1\t50\tsite4\t9\t.\n"         # contig without any gene
)
EXON_HITS = "chr1\t100\t200\tsite1\t42\t+\n"
GENE_HITS = (
    "chr1\t100\t200\tsite1\t42\t+\n"
    "chr1\t500\t600\tsite2\t7\t-\n"
)
GENES_BED = "chr1\t90\t700\tg1\t0\t+\n"
EXONS_BED = "chr1\t90\t200\tg1\t0\t+\n"
GENE_TABLE = "gene_id\tgene_name\tgene_biotype\ng1\tAlpha\tprotein_coding\n"
CLOSEST_BED6 = (
    "chr1\t100\t200\tsite1\t42\t+\tchr1\t90\t700\tg1\t0\t+\t0\n"
    "chr1\t500\t600\tsite2\t7\t-\tchr1\t90\t700\tg1\t0\t+\t0\n"
    "chr1\t5000\t5100\tsite3\t3\t+\tchr1\t90\t700\tg1\t0\t+\t-4300\n"
    "chrX\t1\t50\tsite4\t9\t.\t-1\t-1\t-1\t-1\t-1\t-1\t-1\n"
)


def _annotation_inputs(tmp_path):
    return {
        "peaks": _write(tmp_path, "peaks.bed", PEAKS_BED6),
        "exon_hits": _write(tmp_path, "exons.u.bed", EXON_HITS),
        "gene_hits": _write(tmp_path, "genes.u.bed", GENE_HITS),
        "closest": _write(tmp_path, "closest.tsv", CLOSEST_BED6),
        "table": _write(tmp_path, "genes.tsv", GENE_TABLE),
    }


def _build_rows(tmp_path):
    paths = _annotation_inputs(tmp_path)
    return ann.annotate(
        ann.read_peaks(paths["peaks"], 5),
        ann.read_hits(paths["exon_hits"]),
        ann.read_hits(paths["gene_hits"]),
        ann.read_closest(paths["closest"], 6),
        ann.read_gene_table(paths["table"]))


def test_annotate_feature_classes(tmp_path):
    rows = _build_rows(tmp_path)
    assert [r["feature_class"] for r in rows] == \
        ["exon", "gene", "intergenic", "intergenic"]


def test_annotate_nearest_gene_fields(tmp_path):
    rows = _build_rows(tmp_path)
    assert (rows[0]["nearest_gene"], rows[0]["nearest_gene_id"],
            rows[0]["gene_biotype"], rows[0]["distance"]) == \
        ("Alpha", "g1", "protein_coding", "0")
    assert rows[0]["score"] == "42"


def test_annotate_signed_distance_and_ungene_contig(tmp_path):
    rows = _build_rows(tmp_path)
    assert rows[2]["distance"] == "-4300"      # carried through unchanged
    assert rows[3]["nearest_gene"] == "."
    assert rows[3]["nearest_gene_id"] == "."
    assert rows[3]["distance"] == "NA"
    assert rows[3]["gene_biotype"] == "."


def test_annotate_gene_without_table_entry_falls_back_to_id(tmp_path):
    paths = _annotation_inputs(tmp_path)
    paths["table"] = _write(tmp_path, "bare_table.tsv",
                            "gene_id\tgene_name\tgene_biotype\ng9\tOther\t.\n")
    rows = ann.annotate(
        ann.read_peaks(paths["peaks"], 5),
        ann.read_hits(paths["exon_hits"]),
        ann.read_hits(paths["gene_hits"]),
        ann.read_closest(paths["closest"], 6),
        ann.read_gene_table(paths["table"]))
    assert rows[0]["nearest_gene"] == "g1"      # unknown id -> used verbatim
    assert rows[0]["gene_biotype"] == "."


def test_annotate_consensus_bed4_score_is_support(tmp_path):
    paths = _annotation_inputs(tmp_path)
    peaks4 = _write(tmp_path, "consensus.bed",
                    "chr1\t100\t200\t2\nchr1\t5000\t5100\t3\n")
    closest4 = _write(tmp_path, "closest4.tsv",
                      "chr1\t100\t200\t2\tchr1\t90\t700\tg1\t0\t+\t0\n"
                      "chr1\t5000\t5100\t3\tchr1\t90\t700\tg1\t0\t+\t-4300\n")
    rows = ann.annotate(
        ann.read_peaks(peaks4, 4),
        ann.read_hits(paths["exon_hits"]),
        ann.read_hits(paths["gene_hits"]),
        ann.read_closest(closest4, 4),
        ann.read_gene_table(paths["table"]))
    assert [(r["score"], r["feature_class"]) for r in rows] == \
        [("2", "exon"), ("3", "intergenic")]


def test_annotate_missing_closest_row_raises(tmp_path):
    paths = _annotation_inputs(tmp_path)
    trimmed = _write(tmp_path, "trimmed.tsv",
                     "\n".join(CLOSEST_BED6.strip().splitlines()[:3]) + "\n")
    with pytest.raises(ann.AnnotationError):
        ann.annotate(ann.read_peaks(paths["peaks"], 5),
                     ann.read_hits(paths["exon_hits"]),
                     ann.read_hits(paths["gene_hits"]),
                     ann.read_closest(trimmed, 6),
                     ann.read_gene_table(paths["table"]))


def test_annotate_bad_closest_width_raises(tmp_path):
    with pytest.raises(ann.AnnotationError):
        ann.read_closest(_write(tmp_path, "short.tsv",
                                "chr1\t100\t200\tsite1\t42\t+\tchr1\t90\t700\n"), 6)


def test_read_peaks_rejects_short_lines(tmp_path):
    with pytest.raises(ann.AnnotationError):
        ann.read_peaks(_write(tmp_path, "bad.bed", "chr1\t100\t200\n"), 5)


def test_read_gene_table_skips_header_and_fills_blanks(tmp_path):
    table = ann.read_gene_table(_write(
        tmp_path, "t.tsv",
        "gene_id\tgene_name\tgene_biotype\ng1\t\tprotein_coding\n"))
    assert table == {"g1": {"name": "g1", "biotype": "protein_coding"}}


def test_annotate_main_end_to_end(tmp_path):
    paths = _annotation_inputs(tmp_path)
    out = str(tmp_path / "out.annotation.tsv")
    code = ann.main(["--peaks", paths["peaks"],
                     "--exon-hits", paths["exon_hits"],
                     "--gene-hits", paths["gene_hits"],
                     "--closest", paths["closest"],
                     "--gene-table", paths["table"],
                     "--peak-cols", "6", "--score-col", "5",
                     "--out", out])
    assert code == 0
    lines = open(out).read().splitlines()
    assert lines[0] == "\t".join(ann.TSV_HEADER)
    assert len(lines) == 5                     # header + 4 peaks
    assert lines[1].split("\t")[7] == "exon"   # feature_class column
    assert lines[4].split("\t")[6] == "NA"     # distance column
