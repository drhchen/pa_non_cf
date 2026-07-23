#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: 00_validate_inputs.sh --assemblies DIR --metadata FILE [--manifest FILE]

Validate assembly FASTA files and their one-to-one mapping to a tab-separated
public metadata table. The default manifest is results/provenance/assemblies.tsv.
EOF
}

ASSEMBLIES=""
METADATA=""
MANIFEST="results/provenance/assemblies.tsv"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --assemblies) ASSEMBLIES=${2:?}; shift 2 ;;
        --metadata) METADATA=${2:?}; shift 2 ;;
        --manifest) MANIFEST=${2:?}; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n "$ASSEMBLIES" && -n "$METADATA" ]] || { usage >&2; exit 2; }

python3 scripts/python/validate_inputs.py \
    --assemblies "$ASSEMBLIES" \
    --metadata "$METADATA" \
    --manifest "$MANIFEST"

