#!/usr/bin/env python3
"""Combine per-sample ABRicate tables without mixing databases or summaries."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


OUTPUT_FIELDS = [
    "sample",
    "database",
    "sequence",
    "start",
    "end",
    "strand",
    "gene",
    "coverage",
    "percent_coverage",
    "percent_identity",
    "accession",
    "product",
    "category",
]


def get(row: dict[str, str], name: str) -> str:
    return (row.get(name) or "").strip()


def category_for(database: str, row: dict[str, str]) -> str:
    if database.lower() == "vfdb":
        codes = sorted(set(re.findall(r"\bVFC\d+\b", get(row, "PRODUCT"))))
        return ";".join(codes)
    if database.lower() == "card":
        return get(row, "RESISTANCE") or get(row, "PRODUCT")
    return get(row, "RESISTANCE")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    files = sorted(args.input.glob("*.tab"))
    if not files:
        raise FileNotFoundError(f"no .tab files found in {args.input}")
    args.output.parent.mkdir(parents=True, exist_ok=True)

    total = 0
    with args.output.open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(
            output,
            fieldnames=OUTPUT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for path in files:
            sample = path.stem
            with path.open(newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames is None:
                    raise ValueError(f"missing ABRicate header: {path}")
                required = {"SEQUENCE", "START", "END", "GENE", "%COVERAGE", "%IDENTITY"}
                missing = required.difference(reader.fieldnames)
                if missing:
                    raise ValueError(
                        f"{path} is missing ABRicate columns: {', '.join(sorted(missing))}"
                    )
                for row in reader:
                    gene = get(row, "GENE")
                    if not gene:
                        continue
                    writer.writerow(
                        {
                            "sample": sample,
                            "database": args.database,
                            "sequence": get(row, "SEQUENCE"),
                            "start": get(row, "START"),
                            "end": get(row, "END"),
                            "strand": get(row, "STRAND"),
                            "gene": gene,
                            "coverage": get(row, "COVERAGE"),
                            "percent_coverage": get(row, "%COVERAGE"),
                            "percent_identity": get(row, "%IDENTITY"),
                            "accession": get(row, "ACCESSION"),
                            "product": get(row, "PRODUCT"),
                            "category": category_for(args.database, row),
                        }
                    )
                    total += 1
    print(f"Collated {total} hits from {len(files)} {args.database} tables.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
