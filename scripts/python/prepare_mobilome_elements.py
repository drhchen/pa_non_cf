#!/usr/bin/env python3
"""Extract mobilome records into stable one-sequence files and a manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
from pathlib import Path
from typing import Iterator


def fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    header: str | None = None
    sequence: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header = line[1:]
                sequence = []
            elif header is None:
                raise ValueError(f"sequence before FASTA header in {path}")
            else:
                sequence.append(line)
    if header is not None:
        yield header, "".join(sequence)


def safe_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._-")
    return token or "unnamed"


def classify_element(header: str) -> str:
    lowered = header.lower()
    if "conjugative_integron" in lowered or re.search(r"\bice\b", lowered):
        return "ICE"
    if "mobilizable_integron" in lowered or re.search(r"\bime\b", lowered):
        return "IME"
    return "unknown"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--file-pattern", default="*mobilome.fasta")
    parser.add_argument(
        "--record-pattern",
        default="integron",
        help="case-insensitive header substring; use an empty string for all records",
    )
    args = parser.parse_args()

    source_files = sorted(args.input.rglob(args.file_pattern))
    if not source_files:
        raise FileNotFoundError(
            f"no files matching {args.file_pattern!r} under {args.input}"
        )

    args.output.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output / "elements.tsv"
    matched = 0
    with manifest_path.open("w", newline="", encoding="utf-8") as manifest:
        writer = csv.writer(manifest, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "element_id",
                "sample",
                "element_type",
                "source_file",
                "original_record",
                "length_bp",
                "sha256",
            ]
        )
        for source in source_files:
            sample = safe_token(source.parent.name)
            for header, sequence in fasta_records(source):
                if args.record_pattern.lower() not in header.lower():
                    continue
                if not sequence:
                    raise ValueError(f"empty FASTA record {header!r} in {source}")
                record_token = safe_token(header.split()[0])
                identity = f"{source.resolve()}\0{header}".encode()
                short_hash = hashlib.sha256(identity).hexdigest()[:10]
                element_id = f"{sample}__{record_token}__{short_hash}"
                sequence_hash = hashlib.sha256(sequence.upper().encode()).hexdigest()
                output_fasta = args.output / f"{element_id}.fasta"
                with output_fasta.open("w", encoding="ascii") as handle:
                    handle.write(f">{element_id}\n")
                    for start in range(0, len(sequence), 80):
                        handle.write(sequence[start : start + 80].upper() + "\n")
                writer.writerow(
                    [
                        element_id,
                        sample,
                        classify_element(header),
                        source.relative_to(args.input),
                        header,
                        len(sequence),
                        sequence_hash,
                    ]
                )
                matched += 1

    if matched == 0:
        manifest_path.unlink(missing_ok=True)
        raise ValueError(
            f"no FASTA records matched header substring {args.record_pattern!r}"
        )
    print(f"Extracted {matched} elements from {len(source_files)} files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

