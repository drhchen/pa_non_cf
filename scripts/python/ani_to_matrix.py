#!/usr/bin/env python3
"""Convert standard five-column FastANI output to a symmetric ANI matrix."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def sample_name(value: str) -> str:
    name = Path(value).name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    for suffix in (".fasta", ".fna", ".fa"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return name


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    values: dict[tuple[str, str], list[float]] = {}
    samples: set[str] = set()
    with args.input.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"line {line_number} has fewer than three columns")
            query = sample_name(fields[0])
            reference = sample_name(fields[1])
            ani = float(fields[2])
            key = tuple(sorted((query, reference)))
            values.setdefault(key, []).append(ani)
            samples.update((query, reference))

    ordered = sorted(samples)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", *ordered])
        for row_sample in ordered:
            row = [row_sample]
            for column_sample in ordered:
                if row_sample == column_sample:
                    row.append("100")
                    continue
                pair_values = values.get(tuple(sorted((row_sample, column_sample))))
                row.append("NA" if not pair_values else f"{sum(pair_values) / len(pair_values):.6f}")
            writer.writerow(row)
    print(f"Wrote a {len(ordered)} x {len(ordered)} ANI matrix.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

