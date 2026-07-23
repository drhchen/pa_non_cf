#!/usr/bin/env python3
"""Subset a TSV metadata table to samples in a validated manifest."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with args.manifest.open(newline="", encoding="utf-8") as handle:
        samples = [row["sample"] for row in csv.DictReader(handle, delimiter="\t")]
    with args.metadata.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "sample" not in reader.fieldnames:
            raise ValueError("metadata must contain a sample column")
        rows = {row["sample"]: row for row in reader}

    missing = sorted(set(samples).difference(rows))
    if missing:
        raise ValueError("manifest samples absent from metadata: " + ", ".join(missing))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=reader.fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for sample in samples:
            writer.writerow(rows[sample])
    print(f"Wrote metadata for {len(samples)} retained samples.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

