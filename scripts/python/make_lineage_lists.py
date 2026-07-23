#!/usr/bin/env python3
"""Create lineage-specific assembly lists from validated metadata."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--lineages", nargs="+", default=["BAPS1", "BAPS2"])
    args = parser.parse_args()

    with args.manifest.open(newline="", encoding="utf-8") as handle:
        manifest = {
            row["sample"]: row["assembly"]
            for row in csv.DictReader(handle, delimiter="\t")
        }
    grouped: dict[str, list[tuple[str, str]]] = defaultdict(list)
    with args.metadata.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            group = row["phylogroup"]
            sample = row["sample"]
            if group in args.lineages:
                if sample not in manifest:
                    raise ValueError(f"sample missing from manifest: {sample}")
                grouped[group].append((sample, manifest[sample]))

    args.output.mkdir(parents=True, exist_ok=True)
    for lineage in args.lineages:
        records = sorted(grouped[lineage])
        if not records:
            raise ValueError(f"no samples found for {lineage}")
        with (args.output / f"{lineage}.input.tsv").open(
            "w", newline="", encoding="utf-8"
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerows(records)
        with (args.output / f"{lineage}.samples.txt").open(
            "w", encoding="utf-8"
        ) as handle:
            for sample, _ in records:
                handle.write(f"{sample}\n")
        print(f"{lineage}: {len(records)} samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

