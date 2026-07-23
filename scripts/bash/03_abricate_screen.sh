#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: 03_abricate_screen.sh --config FILE [--force]

Screen every normalized assembly separately with each configured ABRicate
database, summarize hits, and perform lineage-aware statistical analysis.
EOF
}

CONFIG=""
FORCE=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config) CONFIG=${2:?}; shift 2 ;;
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
: "${ABRICATE_DATABASES:?Set ABRICATE_DATABASES in the config}"
: "${ABRICATE_MIN_ID:?Set ABRICATE_MIN_ID in the config}"
: "${ABRICATE_MIN_COVERAGE:?Set ABRICATE_MIN_COVERAGE in the config}"

for executable in abricate python3 Rscript; do
    command -v "$executable" >/dev/null 2>&1 || {
        echo "Required executable not found: $executable" >&2
        exit 127
    }
done

ASSEMBLIES="$RESULTS_DIR/qc_annotation/assemblies"
RETAINED="$RESULTS_DIR/provenance/retained_assemblies.tsv"
RETAINED_METADATA="$RESULTS_DIR/provenance/retained_metadata.tsv"
OUT="$RESULTS_DIR/abricate"
[[ -d "$ASSEMBLIES" ]] || { echo "Run module 01 first: $ASSEMBLIES not found" >&2; exit 1; }
[[ -s "$RETAINED" ]] || { echo "Retained assembly manifest not found: $RETAINED" >&2; exit 1; }
[[ -s "$RETAINED_METADATA" ]] || { echo "Retained metadata not found: $RETAINED_METADATA" >&2; exit 1; }
mkdir -p "$OUT" "$RESULTS_DIR/provenance"
abricate --list > "$RESULTS_DIR/provenance/abricate_databases.tsv"
abricate --version > "$RESULTS_DIR/provenance/abricate_version.txt" 2>&1

for database in $ABRICATE_DATABASES; do
    DB_OUT="$OUT/$database"
    if [[ -d "$DB_OUT" && -n "$(find "$DB_OUT" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]; then
        if [[ "$FORCE" -eq 1 ]]; then
            rm -rf "$DB_OUT"
        else
            echo "Output directory is not empty; use --force: $DB_OUT" >&2
            exit 1
        fi
    fi
    mkdir -p "$DB_OUT/hits"
    while IFS=$'\t' read -r sample _assembly; do
        [[ "$sample" == "sample" ]] && continue
        assembly="$ASSEMBLIES/$sample.fna"
        [[ -s "$assembly" ]] || { echo "Normalized retained assembly not found: $assembly" >&2; exit 1; }
        abricate \
            --db "$database" \
            --minid "$ABRICATE_MIN_ID" \
            --mincov "$ABRICATE_MIN_COVERAGE" \
            --quiet \
            "$assembly" > "$DB_OUT/hits/$sample.tab"
    done < "$RETAINED"

    HIT_FILES=()
    while IFS= read -r hit_file; do
        HIT_FILES+=("$hit_file")
    done < <(find "$DB_OUT/hits" -type f -name '*.tab' -print | LC_ALL=C sort)
    [[ ${#HIT_FILES[@]} -gt 0 ]] || { echo "No hit tables generated for $database" >&2; exit 1; }
    abricate --summary "${HIT_FILES[@]}" > "$DB_OUT/summary.tsv"
    python3 scripts/python/collate_abricate_hits.py \
        --database "$database" \
        --input "$DB_OUT/hits" \
        --output "$DB_OUT/hits_long.tsv"
    ANALYSIS_ARGS=(
        "$DB_OUT/hits_long.tsv"
        "$RETAINED_METADATA"
        "$DB_OUT/analysis"
    )
    if [[ -f config/category_overrides.tsv ]]; then
        ANALYSIS_ARGS+=("config/category_overrides.tsv")
    fi
    Rscript scripts/R/analyze_abricate.R "${ANALYSIS_ARGS[@]}"
done
