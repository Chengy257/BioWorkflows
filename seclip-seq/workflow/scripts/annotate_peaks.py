#!/usr/bin/env python3
"""Merge bedtools intersect/closest outputs into the final peak annotation TSV.

Inputs (produced by the annotate_*_peaks rules in rules/consensus.smk):

  --peaks       the peak BED: PureCLIP output is 7 columns (BED6 chrom,
                start, end, site name, crosslink-site score, strand — plus a
                trailing score-attributes field, verified against a 40k-read
                real run 2026-09-08); consensus BEDs are 4-column (chrom,
                start, end, support)
  --exon-hits   bedtools intersect -a peaks -b exons.bed -u
  --gene-hits   bedtools intersect -a peaks -b genes.bed -u
  --closest     bedtools closest -a peaks -b genes.bed -d -t first:
                <peak columns>, then the 6 gene BED columns, then distance
  --gene-table  gene_id/gene_name/gene_biotype table written by
                gtf_to_gene_regions.py

Output TSV (header row included), one row per input peak in file order:

  chrom  start  end  score  nearest_gene  nearest_gene_id  distance
  feature_class  gene_biotype

feature_class: "exon" (the peak overlaps an exon), "gene" (inside a gene
body but not an exon), "intergenic" (neither). distance is the signed
bedtools closest -d distance to the nearest gene (0 = overlapping; negative
= the gene lies upstream of the peak). Peaks whose contig carries no gene at
all (bedtools closest reports -1 fields) get "." identifiers and an "NA"
distance.

--peak-cols is the number of leading BED columns the peak occupies in the
bedtools outputs (7 for PureCLIP's BED6 + attributes, 4 for consensus
BEDs); --score-col is the 1-based peak column reported as the annotation
score (5 = PureCLIP crosslink-site score, 4 = consensus support).

Pure standard library; exits non-zero on IO/parse errors.
"""
import argparse
import sys

TSV_HEADER = ("chrom", "start", "end", "score", "nearest_gene",
              "nearest_gene_id", "distance", "feature_class", "gene_biotype")


class AnnotationError(Exception):
    """Malformed or mismatched annotation inputs (message is user-facing)."""


def read_peaks(path, score_col):
    """Peak BED -> list of {chrom, start, end, score} in file order.

    score_col is the 1-based column reported as the annotation score."""
    if score_col < 1:
        raise AnnotationError("--score-col must be >= 1, got %d" % score_col)
    peaks = []
    with open(path, "r") as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.rstrip("\r\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < max(3, score_col):
                raise AnnotationError(
                    "%s line %d: expected >= %d columns, got %d"
                    % (path, lineno, max(3, score_col), len(fields)))
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                raise AnnotationError(
                    "%s line %d: non-integer coordinates %r/%r"
                    % (path, lineno, fields[1], fields[2]))
            peaks.append({"chrom": fields[0], "start": start, "end": end,
                          "score": fields[score_col - 1]})
    return peaks


def read_hits(path):
    """bedtools intersect -u output -> set of (chrom, start, end) peak keys."""
    hits = set()
    with open(path, "r") as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.rstrip("\r\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise AnnotationError(
                    "%s line %d: expected >= 3 columns, got %d"
                    % (path, lineno, len(fields)))
            try:
                hits.add((fields[0], int(fields[1]), int(fields[2])))
            except ValueError:
                raise AnnotationError(
                    "%s line %d: non-integer coordinates %r/%r"
                    % (path, lineno, fields[1], fields[2]))
    return hits


def read_closest(path, peak_cols):
    """bedtools closest -d output -> {peak key: (gene_id | None, distance)}.

    Layout: <peak_cols> peak columns, the 6 gene BED columns, then the signed
    distance. Rows without any gene on the peak's contig carry -1 gene
    fields; they map to (None, "NA")."""
    if peak_cols < 1:
        raise AnnotationError("--peak-cols must be >= 1, got %d" % peak_cols)
    expected = peak_cols + 7
    closest = {}
    with open(path, "r") as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != expected:
                raise AnnotationError(
                    "%s line %d: expected %d columns (--peak-cols %d + 6 gene "
                    "columns + distance), got %d"
                    % (path, lineno, expected, peak_cols, len(fields)))
            try:
                key = (fields[0], int(fields[1]), int(fields[2]))
            except ValueError:
                raise AnnotationError(
                    "%s line %d: non-integer coordinates %r/%r"
                    % (path, lineno, fields[1], fields[2]))
            gene_chrom = fields[peak_cols]
            gene_id = fields[peak_cols + 3]
            distance_s = fields[peak_cols + 6]
            if gene_chrom == "-1":
                closest[key] = (None, "NA")
                continue
            try:
                distance = int(distance_s)
            except ValueError:
                raise AnnotationError(
                    "%s line %d: non-integer distance %r" % (path, lineno, distance_s))
            closest[key] = (gene_id, str(distance))
    return closest


def read_gene_table(path):
    """Gene table -> {gene_id: {name, biotype}} (skips the header row)."""
    table = {}
    with open(path, "r") as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if lineno == 1 and fields and fields[0] == "gene_id":
                continue
            if len(fields) < 3:
                raise AnnotationError(
                    "%s line %d: expected 3 columns (gene_id, gene_name, "
                    "gene_biotype), got %d" % (path, lineno, len(fields)))
            table[fields[0]] = {"name": fields[1] or fields[0],
                                "biotype": fields[2] or "."}
    return table


def annotate(peaks, exon_hits, gene_hits, closest, table):
    """Merge everything into output rows (dicts) in peak file order."""
    rows = []
    for p in peaks:
        key = (p["chrom"], p["start"], p["end"])
        if key not in closest:
            raise AnnotationError(
                "peak %s:%d-%d has no bedtools closest row; the --closest file "
                "does not match the --peaks file" % key)
        gene_id, distance = closest[key]
        if gene_id is None:
            nearest_gene = nearest_id = "."
            biotype = "."
        else:
            entry = table.get(gene_id)
            nearest_id = gene_id
            nearest_gene = entry["name"] if entry else gene_id
            biotype = entry["biotype"] if entry else "."
        if key in exon_hits:
            feature_class = "exon"
        elif key in gene_hits:
            feature_class = "gene"
        else:
            feature_class = "intergenic"
        rows.append({"chrom": p["chrom"], "start": p["start"], "end": p["end"],
                     "score": p["score"], "nearest_gene": nearest_gene,
                     "nearest_gene_id": nearest_id, "distance": distance,
                     "feature_class": feature_class, "gene_biotype": biotype})
    return rows


def write_tsv(rows, out_path):
    with open(out_path, "w") as fh:
        fh.write("\t".join(TSV_HEADER) + "\n")
        for r in rows:
            fh.write("\t".join([r["chrom"], str(r["start"]), str(r["end"]),
                                r["score"], r["nearest_gene"],
                                r["nearest_gene_id"], r["distance"],
                                r["feature_class"], r["gene_biotype"]]) + "\n")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Merge bedtools intersect/closest outputs into the "
                    "seclip-seq peak annotation TSV.")
    parser.add_argument("--peaks", required=True,
                        help="peak BED (7-column PureCLIP BED6+attributes, or 4-column consensus)")
    parser.add_argument("--exon-hits", required=True,
                        help="bedtools intersect -a peaks -b exons.bed -u output")
    parser.add_argument("--gene-hits", required=True,
                        help="bedtools intersect -a peaks -b genes.bed -u output")
    parser.add_argument("--closest", required=True,
                        help="bedtools closest -a peaks -b genes.bed -d -t first output")
    parser.add_argument("--gene-table", required=True,
                        help="gene_id/gene_name/gene_biotype TSV from gtf_to_gene_regions.py")
    parser.add_argument("--peak-cols", type=int, default=7,
                        help="peak columns in the bedtools outputs (7 = PureCLIP BED6+attributes, 4 = consensus)")
    parser.add_argument("--score-col", type=int, default=5,
                        help="1-based peak column used as the score (5 = PureCLIP score, 4 = support)")
    parser.add_argument("--out", required=True, help="output annotation TSV")
    args = parser.parse_args(argv)
    try:
        peaks = read_peaks(args.peaks, args.score_col)
        exon_hits = read_hits(args.exon_hits)
        gene_hits = read_hits(args.gene_hits)
        closest = read_closest(args.closest, args.peak_cols)
        table = read_gene_table(args.gene_table)
        rows = annotate(peaks, exon_hits, gene_hits, closest, table)
        write_tsv(rows, args.out)
    except (AnnotationError, OSError) as exc:
        print("annotate_peaks: error: %s" % exc, file=sys.stderr)
        return 1
    counts = {}
    for r in rows:
        counts[r["feature_class"]] = counts.get(r["feature_class"], 0) + 1
    summary = ", ".join("%s %d" % (k, counts[k]) for k in sorted(counts))
    print("[annotate_peaks] peaks %d (%s) -> %s" % (len(rows), summary or "none", args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
