#!/usr/bin/env python3
"""Validate the RNA-seq sample metadata table before workflow execution.

The table must contain ``id`` and ``group`` columns. Optional columns include
``layout`` (PE, SE, or auto) and ``batch``. The legacy invocation
``validate_samples.py samples.csv control`` remains supported.
"""

import argparse
import csv
import sys
from collections import Counter
from itertools import combinations
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate RNA-seq sample metadata and report actionable errors."
    )
    parser.add_argument("sample_table", help="CSV sample metadata file")
    parser.add_argument(
        "control_group",
        nargs="?",
        default="control",
        help="Expected control-group name (default: control)",
    )
    parser.add_argument(
        "--require-batch",
        action="store_true",
        help="Require a non-empty batch column for every sample",
    )
    parser.add_argument(
        "--strict-warnings",
        action="store_true",
        help="Return a non-zero exit status when warnings are present",
    )
    return parser.parse_args()


def emit(errors, warnings, summary):
    for message in warnings:
        print(f"[WARN ] {message}", file=sys.stderr)
    for message in errors:
        print(f"[ERROR] {message}", file=sys.stderr)
    print(f"[INFO ] {summary}", file=sys.stderr)
    print(
        f"[validate_samples] {len(errors)} error(s), {len(warnings)} warning(s)",
        file=sys.stderr,
    )


def main():
    args = parse_args()
    path = Path(args.sample_table)
    if not path.is_file():
        print(f"[ERROR] Sample table not found: {path}", file=sys.stderr)
        return 1

    errors, warnings = [], []
    try:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames or []
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as exc:
        print(f"[ERROR] Unable to read sample table: {exc}", file=sys.stderr)
        return 1

    required = ["id", "group"]
    if args.require_batch:
        required.append("batch")
    for column in required:
        if column not in fields:
            errors.append(f"Missing required column '{column}'. Current header: {fields}")

    if errors:
        emit(errors, warnings, f"Validation stopped before row checks for {path}")
        return 1

    if not rows:
        errors.append("Sample table contains no data rows")

    ids = []
    groups = []
    layouts = Counter()
    batches = []

    for line_number, row in enumerate(rows, start=2):
        sample_id = (row.get("id") or "").strip()
        group = (row.get("group") or "").strip()
        layout = (row.get("layout") or "auto").strip().lower() or "auto"
        batch = (row.get("batch") or "").strip()

        if not sample_id:
            errors.append(f"Line {line_number}: sample id is empty")
            continue
        if sample_id in ids:
            errors.append(f"Line {line_number}: duplicate sample id '{sample_id}'")
        ids.append(sample_id)

        if "-" in sample_id:
            errors.append(
                f"Line {line_number}: sample id contains '-': '{sample_id}'. "
                "Hyphens can be renamed by downstream scripts and cause mismatches."
            )
        if any(char.isspace() or char == "/" for char in sample_id):
            errors.append(
                f"Line {line_number}: sample id contains whitespace or '/': '{sample_id}'"
            )

        if not group:
            errors.append(f"Line {line_number}: group is empty for sample '{sample_id}'")
        else:
            groups.append(group)

        if layout not in {"pe", "se", "auto"}:
            warnings.append(
                f"Line {line_number}: layout should be PE, SE, or auto; got '{row.get('layout')}'"
            )
        layouts[layout] += 1

        if args.require_batch and not batch:
            errors.append(
                f"Line {line_number}: batch is required but empty for sample '{sample_id}'"
            )
        if batch:
            batches.append(batch)

    group_set = set(groups)
    if args.control_group not in group_set:
        errors.append(
            f"Control group '{args.control_group}' is absent from the group column. "
            f"Available groups: {sorted(group_set)}"
        )

    for left, right in combinations(ids, 2):
        if left.startswith(right) or right.startswith(left):
            warnings.append(
                f"Sample ids are prefix-related and may be confused by loose filename matching: "
                f"'{left}' vs '{right}'"
            )

    for group, count in sorted(Counter(groups).items()):
        if count < 2:
            warnings.append(
                f"Group '{group}' has only {count} sample(s); differential-expression "
                "analysis without biological replication is unreliable."
            )

    if args.require_batch and batches and len(set(batches)) < 2:
        warnings.append(
            "Batch correction is enabled, but all non-empty batch values are identical."
        )

    layout_summary = ", ".join(f"{key.upper()}={value}" for key, value in sorted(layouts.items()))
    summary = (
        f"Validated {len(rows)} row(s), {len(set(ids))} unique sample(s), "
        f"{len(group_set)} group(s); layouts: {layout_summary or 'none'}"
    )
    emit(errors, warnings, summary)

    if errors:
        return 1
    if warnings and args.strict_warnings:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
