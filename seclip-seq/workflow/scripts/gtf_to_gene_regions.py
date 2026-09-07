#!/usr/bin/env python3
"""Convert a gene-annotation GTF into gene/exon BED files for peak annotation.

Pure standard library. The GTF coordinate system is 1-based closed; the BED
output is 0-based half-open (BED start = GTF start - 1). Outputs:

  genes.bed  BED6: chrom, start, end, name (gene_id), score (0), strand
  exons.bed  BED6: chrom, start, end, name (gene_id), score (0), strand
  genes.tsv  one row per gene: gene_id <TAB> gene_name <TAB> gene_biotype,
             with '.' placeholders for missing attributes (header row
             included)

The gene_id is used as the BED name so the annotation merge step can join on
a stable key; the human-readable gene_name travels in genes.tsv. Gene spans
combine the gene records and their exon records (min start / max end), so
GTFs without explicit gene rows still yield complete gene bodies. gene_name
and gene_biotype are taken from the gene record when present, otherwise from
the exon records; missing gene_id attributes are a hard error.

Usage:
  python3 gtf_to_gene_regions.py --gtf annotation.gtf \
      --genes-out genes.bed --exons-out exons.bed --table-out genes.tsv
"""
import argparse
import os
import sys

TABLE_HEADER = "gene_id\tgene_name\tgene_biotype"


class GtfParseError(Exception):
    """Malformed GTF input (the message carries the file and line number)."""


def parse_attributes(text):
    """Parse a GTF attribute column into a dict.

    Attributes are semicolon-separated "key value" pairs; values may be
    double-quoted (spaces preserved) or bare. Empty pairs are skipped; a key
    without a value maps to the empty string."""
    attrs = {}
    for pair in text.split(";"):
        pair = pair.strip()
        if not pair:
            continue
        parts = pair.split(None, 1)
        key = parts[0]
        value = parts[1].strip() if len(parts) > 1 else ""
        if len(value) >= 2 and value.startswith('"') and value.endswith('"'):
            value = value[1:-1]
        attrs[key] = value
    return attrs


def _merge_span(gene, chrom, strand, start, end, path, lineno, gene_id):
    """Grow a gene span, rejecting chrom/strand conflicts for one gene_id."""
    if gene["chrom"] != chrom or gene["strand"] != strand:
        raise GtfParseError(
            "%s line %d: gene_id %r re-declared on %s (%s) but first seen on "
            "%s (%s)" % (path, lineno, gene_id, chrom, strand,
                         gene["chrom"], gene["strand"]))
    gene["start"] = min(gene["start"], start)
    gene["end"] = max(gene["end"], end)


def parse_gtf(path):
    """Parse the GTF; returns (genes, exons).

    genes: dict gene_id -> {chrom, strand, start, end, name, biotype} with a
           0-based half-open span, name = gene_name or gene_id and
           biotype = gene_biotype or '.'.
    exons: deduplicated list of {chrom, strand, start, end, gene_id}."""
    genes = {}
    exons = {}
    with open(path, "r") as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                raise GtfParseError(
                    "%s line %d: expected >= 9 tab-separated GTF columns, got %d"
                    % (path, lineno, len(fields)))
            chrom, _source, feature, start_s, end_s, _score, strand, _frame, attr_text = fields[:9]
            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                raise GtfParseError(
                    "%s line %d: non-integer coordinates %r/%r"
                    % (path, lineno, start_s, end_s))
            if end < start:
                raise GtfParseError(
                    "%s line %d: end (%d) < start (%d)" % (path, lineno, end, start))
            if strand not in ("+", "-", "."):
                raise GtfParseError(
                    "%s line %d: illegal strand %r" % (path, lineno, strand))
            if feature not in ("gene", "exon"):
                continue
            attrs = parse_attributes(attr_text)
            gene_id = attrs.get("gene_id", "")
            if not gene_id:
                raise GtfParseError(
                    "%s line %d: %s record without a gene_id attribute"
                    % (path, lineno, feature))
            bed_start = start - 1
            gene = genes.get(gene_id)
            if feature == "exon":
                exons[(chrom, bed_start, end, gene_id, strand)] = None
                if gene is None:
                    genes[gene_id] = {"chrom": chrom, "strand": strand,
                                      "start": bed_start, "end": end,
                                      "name": gene_id,
                                      "biotype": attrs.get("gene_biotype", ".")}
                else:
                    _merge_span(gene, chrom, strand, bed_start, end,
                                path, lineno, gene_id)
                    if gene["biotype"] == "." and attrs.get("gene_biotype"):
                        gene["biotype"] = attrs["gene_biotype"]
            else:
                name = attrs.get("gene_name") or gene_id
                biotype = attrs.get("gene_biotype") or "."
                if gene is None:
                    genes[gene_id] = {"chrom": chrom, "strand": strand,
                                      "start": bed_start, "end": end,
                                      "name": name, "biotype": biotype}
                else:
                    _merge_span(gene, chrom, strand, bed_start, end,
                                path, lineno, gene_id)
                    gene["name"] = name          # gene-row attributes win
                    if biotype != ".":
                        gene["biotype"] = biotype
    exon_list = [{"chrom": c, "start": s, "end": e, "gene_id": g, "strand": st}
                 for (c, s, e, g, st) in exons]
    return genes, exon_list


def write_outputs(genes, exons, genes_out, exons_out, table_out):
    """Write the gene/exon BEDs (sorted by chrom, start) and the gene table."""
    with open(genes_out, "w") as fh:
        for gid in sorted(genes):
            g = genes[gid]
            fh.write("%s\t%d\t%d\t%s\t0\t%s\n"
                     % (g["chrom"], g["start"], g["end"], gid, g["strand"]))
    with open(exons_out, "w") as fh:
        for e in sorted(exons, key=lambda x: (x["chrom"], x["start"], x["end"],
                                              x["gene_id"])):
            fh.write("%s\t%d\t%d\t%s\t0\t%s\n"
                     % (e["chrom"], e["start"], e["end"], e["gene_id"], e["strand"]))
    with open(table_out, "w") as fh:
        fh.write(TABLE_HEADER + "\n")
        for gid in sorted(genes):
            g = genes[gid]
            fh.write("%s\t%s\t%s\n" % (gid, g["name"], g["biotype"]))


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Convert a GTF into gene/exon BEDs plus a gene attribute "
                    "table for the seclip-seq peak annotation stage.")
    parser.add_argument("--gtf", required=True, help="input GTF (1-based)")
    parser.add_argument("--genes-out", required=True,
                        help="output gene-body BED6 (0-based half-open)")
    parser.add_argument("--exons-out", required=True, help="output exon BED6")
    parser.add_argument("--table-out", required=True,
                        help="output gene_id/gene_name/gene_biotype TSV")
    args = parser.parse_args(argv)
    try:
        genes, exons = parse_gtf(args.gtf)
        write_outputs(genes, exons, args.genes_out, args.exons_out, args.table_out)
    except (GtfParseError, OSError) as exc:
        print("gtf_to_gene_regions: error: %s" % exc, file=sys.stderr)
        return 1
    print("[gtf_to_gene_regions] genes %d, exons %d; wrote %s, %s, %s"
          % (len(genes), len(exons), args.genes_out, args.exons_out, args.table_out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
