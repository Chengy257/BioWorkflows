#!/usr/bin/env python3
"""Prepare nf-core/differentialabundance inputs from the rnaseq route.

Emits the three DA input files for one project:
    <outdir>/matrix.tsv         gene x sample count matrix (switch-selected)
    <outdir>/samplesheet.csv    DA samplesheet (sample + analysis columns)
    <outdir>/contrasts.yaml     one contrast per non-control group (legacy
                                semantics: each treatment vs the control
                                group; reference = control, target = group)

Matrix sources (decision record D1 quantification-route switch):
    --route salmon        the adopted run's salmon.merged.gene_counts.tsv
                          (gene_id, gene_name, <samples>...) - passed through
    --route featurecounts the rnaseq-ext count.matrix.tsv (gene_id, <samples>...)

Usage:
    python prepare_da_inputs.py --metadata meta.csv --matrix matrix.tsv \
        --route salmon --control-group control --outdir da_inputs \
        [--variable group] [--blocking batch] [--contrast groupA_vs_groupB]
"""
import argparse
import csv
import os
import sys


def read_table(path):
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t" if path.endswith((".tsv", ".txt")) else ",")
        header = next(reader)
        rows = [r for r in reader if any(c.strip() for c in r)]
    return header, rows


def fail(msg):
    sys.exit(f"[ERROR] {msg}")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("--metadata", required=True, help="analysis metadata CSV: sample,group[,batch]")
    p.add_argument("--matrix", required=True, help="gene x sample count matrix (route-selected)")
    p.add_argument("--route", required=True, choices=["salmon", "featurecounts"],
                   help="authoritative-matrix route (decision record D1)")
    p.add_argument("--control-group", required=True, help="the control group name (legacy control_group)")
    p.add_argument("--outdir", required=True)
    p.add_argument("--variable", default="group", help="samplesheet column holding the contrast variable")
    p.add_argument("--blocking", default="", help="optional samplesheet column used as the blocking factor (batch)")
    p.add_argument("--contrast", action="append", default=[],
                   help="explicit extra contrast '<id>:<reference>:<target>' (repeatable); "
                        "ids must not repeat the auto-derived ones")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ---- Metadata (DA samplesheet source) -----------------------------------
    meta_header, meta_rows = read_table(args.metadata)
    for col in ("sample", args.variable):
        if col not in meta_header:
            fail(f"metadata {args.metadata}: missing column '{col}' (header: {meta_header})")
    if args.blocking and args.blocking not in meta_header:
        fail(f"metadata {args.metadata}: --blocking column '{args.blocking}' not in header {meta_header}")

    ids = [r[meta_header.index("sample")].strip() for r in meta_rows]
    if len(ids) != len(set(ids)):
        fail("duplicate sample ids in the metadata sheet")
    groups = {}
    for r in meta_rows:
        s = r[meta_header.index("sample")].strip()
        groups[s] = r[meta_header.index(args.variable)].strip()

    control_samples = sorted(s for s, g in groups.items() if g == args.control_group)
    if not control_samples:
        fail(f"no samples in the control group '{args.control_group}'")

    # ---- Matrix: validate coverage, subset + order to the metadata ---------
    m_header, m_rows = read_table(args.matrix)
    matrix_samples = m_header[1:]
    missing = [s for s in ids if s not in matrix_samples]
    if missing:
        fail(f"matrix {args.matrix}: no column for metadata sample(s) {missing} "
             f"(matrix columns: {matrix_samples})")
    keep_idx = [m_header.index(s) for s in ids]

    out_matrix = os.path.join(args.outdir, "matrix.tsv")
    with open(out_matrix, "w") as fh:
        fh.write("\t".join(m_header[:1] + ids) + "\n")
        for row in m_rows:
            fh.write("\t".join(row[0:1] + [row[i] for i in keep_idx]) + "\n")

    # ---- DA samplesheet -----------------------------------------------------
    out_sheet = os.path.join(args.outdir, "samplesheet.csv")
    with open(out_sheet, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample", args.variable] + ([args.blocking] if args.blocking else []))
        for r in meta_rows:
            writer.writerow([r[meta_header.index("sample")], r[meta_header.index(args.variable)]]
                            + ([r[meta_header.index(args.blocking)]] if args.blocking else []))

    # ---- Contrasts: each non-control group vs the control (legacy semantics)
    contrasts = []
    for g in sorted({groups[s] for s in ids} - {args.control_group}):
        contrasts.append((f"{g}_vs_{args.control_group}", args.control_group, g))
    for spec in args.contrast:
        parts = spec.split(":")
        if len(parts) != 3:
            fail(f"--contrast '{spec}': expected '<id>:<reference>:<target>'")
        contrasts.append(tuple(parts))
    seen = set()
    for cid, _, _ in contrasts:
        if cid in seen:
            fail(f"duplicate contrast id '{cid}'")
        seen.add(cid)

    out_contrasts = os.path.join(args.outdir, "contrasts.yaml")
    with open(out_contrasts, "w") as fh:
        fh.write("contrasts:\n")
        for cid, ref, tgt in contrasts:
            fh.write(f"  - id: {cid}\n")
            fh.write(f"    comparison: ['{args.variable}', '{ref}', '{tgt}']\n")
            if args.blocking:
                fh.write(f"    blocking_factors: ['{args.blocking}']\n")

    print(f"[prepare_da_inputs] route={args.route} samples={len(ids)} contrasts={len(contrasts)} -> {args.outdir}")
    print(f"  matrix:       {out_matrix}")
    print(f"  samplesheet:  {out_sheet}")
    print(f"  contrasts:    {out_contrasts}")


if __name__ == "__main__":
    main()
