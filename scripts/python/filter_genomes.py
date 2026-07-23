#!/usr/bin/env python3
"""Apply CheckM and reference-ANI thresholds to an assembly manifest."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def normalized_sample(value: str) -> str:
    name = Path(value).name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    for suffix in (".fasta", ".fna", ".fa"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return name


def read_checkm(path: Path) -> dict[str, tuple[float, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("CheckM table has no header")
        lowered = {name.lower().replace(" ", "_"): name for name in reader.fieldnames}
        bin_column = lowered.get("bin_id")
        completeness_column = lowered.get("completeness")
        contamination_column = lowered.get("contamination")
        if not all((bin_column, completeness_column, contamination_column)):
            raise ValueError(
                "CheckM table must contain Bin Id, Completeness, and Contamination"
            )
        result = {}
        for row in reader:
            sample = normalized_sample(row[bin_column])
            result[sample] = (
                float(row[completeness_column]),
                float(row[contamination_column]),
            )
        return result


def read_ani(path: Path) -> dict[tuple[str, str], list[float]]:
    result: dict[tuple[str, str], list[float]] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"FastANI line {line_number} has fewer than three fields")
            query = normalized_sample(fields[0])
            reference = normalized_sample(fields[1])
            key = tuple(sorted((query, reference)))
            result.setdefault(key, []).append(float(fields[2]))
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--checkm", type=Path, required=True)
    parser.add_argument("--ani", type=Path, required=True)
    parser.add_argument("--retained", type=Path, required=True)
    parser.add_argument("--excluded", type=Path, required=True)
    parser.add_argument("--min-completeness", type=float, default=90.0)
    parser.add_argument("--max-contamination", type=float, default=10.0)
    parser.add_argument("--min-ani", type=float, default=95.0)
    args = parser.parse_args()

    checkm = read_checkm(args.checkm)
    ani_pairs = read_ani(args.ani)
    with args.manifest.open(newline="", encoding="utf-8") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))

    args.retained.parent.mkdir(parents=True, exist_ok=True)
    args.excluded.parent.mkdir(parents=True, exist_ok=True)
    retained_count = 0
    excluded_count = 0
    with args.retained.open("w", newline="", encoding="utf-8") as retained_handle, args.excluded.open(
        "w", newline="", encoding="utf-8"
    ) as excluded_handle:
        retained_writer = csv.writer(retained_handle, delimiter="\t", lineterminator="\n")
        excluded_writer = csv.writer(excluded_handle, delimiter="\t", lineterminator="\n")
        retained_writer.writerow(["sample", "assembly"])
        excluded_writer.writerow(
            ["sample", "assembly", "completeness", "contamination", "ani", "reason"]
        )
        manifest_samples = [row["sample"] for row in manifest]
        for row in manifest:
            sample = row["sample"]
            completeness, contamination = checkm.get(sample, (float("nan"), float("nan")))
            comparison_values = []
            missing_comparisons = []
            for other in manifest_samples:
                if other == sample:
                    continue
                values = ani_pairs.get(tuple(sorted((sample, other))))
                if values:
                    comparison_values.append(sum(values) / len(values))
                else:
                    missing_comparisons.append(other)
            ani_value = min(comparison_values) if comparison_values else float("nan")
            reasons = []
            if sample not in checkm:
                reasons.append("missing_checkm")
            else:
                if completeness < args.min_completeness:
                    reasons.append("completeness_below_threshold")
                if contamination > args.max_contamination:
                    reasons.append("contamination_above_threshold")
            if missing_comparisons or not comparison_values:
                reasons.append("missing_ani")
            if comparison_values and ani_value < args.min_ani:
                reasons.append("ani_below_threshold")

            if reasons:
                excluded_writer.writerow(
                    [
                        sample,
                        row["assembly"],
                        completeness,
                        contamination,
                        ani_value,
                        ";".join(reasons),
                    ]
                )
                excluded_count += 1
            else:
                retained_writer.writerow([sample, row["assembly"]])
                retained_count += 1

    if retained_count == 0:
        raise ValueError("all genomes were excluded; inspect threshold tables")
    print(f"Retained {retained_count} genomes; excluded {excluded_count}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
