#!/usr/bin/env python3
"""Validate assembly basenames and public sample metadata."""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path


FASTA_SUFFIXES = (".fa", ".fasta", ".fna")


def sample_from_path(path: Path) -> str:
    name = path.name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    for suffix in FASTA_SUFFIXES:
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    raise ValueError(f"unsupported assembly suffix: {path}")


def find_assemblies(directory: Path) -> list[Path]:
    files: list[Path] = []
    for path in sorted(directory.iterdir()):
        lower = path.name.lower()
        if path.is_file() and any(
            lower.endswith(suffix) or lower.endswith(f"{suffix}.gz")
            for suffix in FASTA_SUFFIXES
        ):
            files.append(path.resolve())
    return files


def validate_fasta(path: Path) -> None:
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    with opener(path, "rt", encoding="ascii", errors="strict") as handle:
        first_nonempty = next((line for line in handle if line.strip()), "")
    if not first_nonempty.startswith(">"):
        raise ValueError(f"not a FASTA file (first non-empty line is not a header): {path}")


def read_metadata(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    raw = path.read_bytes()
    if b"\r\n" in raw:
        raise ValueError(f"Windows CRLF line endings detected: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"metadata has no header: {path}")
        rows = list(reader)
        return reader.fieldnames, rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assemblies", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument(
        "--manifest",
        type=Path,
        help="write a validated sample-to-assembly manifest",
    )
    args = parser.parse_args()

    errors: list[str] = []
    if not args.assemblies.is_dir():
        errors.append(f"assembly directory does not exist: {args.assemblies}")
    if not args.metadata.is_file():
        errors.append(f"metadata file does not exist: {args.metadata}")
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1

    assemblies = find_assemblies(args.assemblies)
    if not assemblies:
        errors.append(f"no supported assemblies found in {args.assemblies}")

    assembly_by_sample: dict[str, Path] = {}
    for assembly in assemblies:
        sample = sample_from_path(assembly)
        if not sample.strip():
            errors.append(f"empty sample ID derived from {assembly}")
        if sample in assembly_by_sample:
            errors.append(
                f"duplicate normalized assembly basename '{sample}': "
                f"{assembly_by_sample[sample]} and {assembly}"
            )
        assembly_by_sample[sample] = assembly
        try:
            validate_fasta(assembly)
        except (OSError, UnicodeError, ValueError) as exc:
            errors.append(str(exc))

    try:
        fields, rows = read_metadata(args.metadata)
    except (OSError, UnicodeError, ValueError) as exc:
        errors.append(str(exc))
        fields, rows = [], []

    required = {"sample", "accession", "phylogroup"}
    missing_fields = sorted(required.difference(fields))
    if missing_fields:
        errors.append(f"metadata is missing columns: {', '.join(missing_fields)}")

    metadata_samples: list[str] = []
    for line_number, row in enumerate(rows, start=2):
        sample = (row.get("sample") or "").strip()
        metadata_samples.append(sample)
        if not sample:
            errors.append(f"metadata line {line_number}: empty sample")
        if any(character.isspace() for character in sample):
            errors.append(f"metadata line {line_number}: sample contains whitespace: {sample!r}")
        group = (row.get("phylogroup") or "").strip()
        if group and group not in {"BAPS1", "BAPS2", "BAPS3", "BAPS4", "BAPS5"}:
            errors.append(
                f"metadata line {line_number}: invalid phylogroup {group!r}; "
                "expected BAPS1-BAPS5 or blank"
            )

    duplicates = sorted({sample for sample in metadata_samples if metadata_samples.count(sample) > 1})
    if duplicates:
        errors.append(f"duplicate metadata sample IDs: {', '.join(duplicates)}")

    missing_assemblies = sorted(set(metadata_samples).difference(assembly_by_sample))
    extra_assemblies = sorted(set(assembly_by_sample).difference(metadata_samples))
    if missing_assemblies:
        errors.append(
            "metadata samples without assemblies: " + ", ".join(missing_assemblies[:20])
        )
    if extra_assemblies:
        errors.append(
            "assemblies absent from metadata: " + ", ".join(extra_assemblies[:20])
        )

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1

    if args.manifest:
        args.manifest.parent.mkdir(parents=True, exist_ok=True)
        with args.manifest.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["sample", "assembly"])
            for sample in metadata_samples:
                writer.writerow([sample, assembly_by_sample[sample]])

    print(
        f"Validated {len(assemblies)} assemblies and {len(rows)} metadata records.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

