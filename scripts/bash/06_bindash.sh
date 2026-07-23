#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: 06_bindash.sh --elements DIR --output DIR [--kmer 21] [--force]

Sketch one-sequence FASTA files from prepare_mobilome_elements.py and calculate
all pairwise BinDash distances. The installed BinDash build must support
--kmerlen; verify its help output before the publication run.
EOF
}

ELEMENTS=""
OUTPUT=""
KMER=21
FORCE=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --elements) ELEMENTS=${2:?}; shift 2 ;;
        --output) OUTPUT=${2:?}; shift 2 ;;
        --kmer) KMER=${2:?}; shift 2 ;;
        --force) FORCE=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done
[[ -d "$ELEMENTS" && -n "$OUTPUT" ]] || { usage >&2; exit 2; }
command -v bindash >/dev/null 2>&1 || { echo "bindash not found" >&2; exit 127; }

if [[ -d "$OUTPUT" && -n "$(find "$OUTPUT" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
    if [[ "$FORCE" -eq 1 ]]; then
        rm -rf "$OUTPUT"
    else
        echo "Output directory is not empty; use --force: $OUTPUT" >&2
        exit 1
    fi
fi
mkdir -p "$OUTPUT"
find "$ELEMENTS" -maxdepth 1 -type f -name '*.fasta' -print | LC_ALL=C sort > "$OUTPUT/elements.list"
[[ -s "$OUTPUT/elements.list" ]] || { echo "No element FASTA files found in $ELEMENTS" >&2; exit 1; }

bindash sketch \
    --listfname="$OUTPUT/elements.list" \
    --outfname="$OUTPUT/elements.sketch" \
    --kmerlen="$KMER"
bindash dist \
    --outfname="$OUTPUT/dist.tsv" \
    "$OUTPUT/elements.sketch"
bindash --version > "$OUTPUT/bindash_version.txt" 2>&1 || bindash --help > "$OUTPUT/bindash_version.txt" 2>&1

