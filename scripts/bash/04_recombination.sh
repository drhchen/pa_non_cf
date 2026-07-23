#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: 04_recombination.sh --config FILE [--lineage BAPS1|BAPS2] [--force]

Create lineage sample lists and run Gubbins on reference-coordinate alignments.
Set BAPS1_ALIGNMENT and BAPS2_ALIGNMENT in the config, or install
generate_ska_alignment.py and provide REFERENCE_FASTA to generate alignments.
EOF
}

CONFIG=""
ONLY_LINEAGE=""
FORCE=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config) CONFIG=${2:?}; shift 2 ;;
        --lineage) ONLY_LINEAGE=${2:?}; shift 2 ;;
        --force) FORCE=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done
[[ -f "$CONFIG" ]] || { echo "Config file not found: $CONFIG" >&2; exit 2; }
# shellcheck source=/dev/null
source "$CONFIG"

: "${RESULTS_DIR:?Set RESULTS_DIR in the config}"
: "${METADATA:?Set METADATA in the config}"
: "${REFERENCE_FASTA:?Set REFERENCE_FASTA in the config}"
: "${THREADS:?Set THREADS in the config}"
: "${GUBBINS_FILTER_PERCENTAGE:?Set GUBBINS_FILTER_PERCENTAGE in the config}"

command -v run_gubbins.py >/dev/null 2>&1 || { echo "run_gubbins.py not found" >&2; exit 127; }
MANIFEST="$RESULTS_DIR/provenance/retained_assemblies.tsv"
[[ -s "$MANIFEST" ]] || { echo "Validated assembly manifest not found: $MANIFEST" >&2; exit 1; }
LIST_DIR="$RESULTS_DIR/recombination/lineage_lists"
python3 scripts/python/make_lineage_lists.py \
    --manifest "$MANIFEST" \
    --metadata "$METADATA" \
    --output "$LIST_DIR"

LINEAGES=(BAPS1 BAPS2)
[[ -z "$ONLY_LINEAGE" ]] || LINEAGES=("$ONLY_LINEAGE")
for lineage in "${LINEAGES[@]}"; do
    [[ "$lineage" == "BAPS1" || "$lineage" == "BAPS2" ]] || {
        echo "Unsupported lineage: $lineage" >&2
        exit 2
    }
    variable_name="${lineage}_ALIGNMENT"
    alignment=${!variable_name:-}
    lineage_out="$RESULTS_DIR/recombination/$lineage"
    if [[ -z "$alignment" ]]; then
        command -v generate_ska_alignment.py >/dev/null 2>&1 || {
            echo "$variable_name is unset and generate_ska_alignment.py is unavailable" >&2
            exit 127
        }
        mkdir -p "$lineage_out"
        alignment="$lineage_out/${lineage}.reference.aln"
        generate_ska_alignment.py \
            --reference "$REFERENCE_FASTA" \
            --input "$LIST_DIR/${lineage}.input.tsv" \
            --out "$alignment" \
            --threads "$THREADS"
    fi
    [[ -s "$alignment" ]] || { echo "Alignment not found or empty: $alignment" >&2; exit 1; }

    if [[ -d "$lineage_out/gubbins" && -n "$(find "$lineage_out/gubbins" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
        if [[ "$FORCE" -eq 1 ]]; then
            rm -rf "$lineage_out/gubbins"
        else
            echo "Gubbins output exists; use --force: $lineage_out/gubbins" >&2
            exit 1
        fi
    fi
    mkdir -p "$lineage_out/gubbins"
    alignment=$(cd "$(dirname "$alignment")" && pwd)/$(basename "$alignment")
    (
        cd "$lineage_out/gubbins"
        run_gubbins.py \
            --prefix "$lineage" \
            --filter-percentage "$GUBBINS_FILTER_PERCENTAGE" \
            --threads "$THREADS" \
            "$alignment"
    )
done
