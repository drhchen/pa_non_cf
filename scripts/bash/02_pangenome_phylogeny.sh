#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: 02_pangenome_phylogeny.sh --config FILE [--force]

Run Panaroo, SNP-sites, and IQ-TREE on Prokka GFF files produced by module 01.
The configured Panaroo mode must match the final manuscript analysis.
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
: "${THREADS:?Set THREADS in the config}"
: "${PANAROO_CLEAN_MODE:?Set PANAROO_CLEAN_MODE in the config}"
: "${CORE_THRESHOLD:?Set CORE_THRESHOLD in the config}"
: "${IQTREE_MODEL:?Set IQTREE_MODEL in the config}"
: "${SEED:?Set SEED in the config}"

for executable in panaroo snp-sites; do
    command -v "$executable" >/dev/null 2>&1 || {
        echo "Required executable not found: $executable" >&2
        exit 127
    }
done
if command -v iqtree2 >/dev/null 2>&1; then
    IQTREE=iqtree2
elif command -v iqtree >/dev/null 2>&1; then
    IQTREE=iqtree
else
    echo "Required executable not found: iqtree2 or iqtree" >&2
    exit 127
fi

GFF_ROOT="$RESULTS_DIR/qc_annotation/prokka"
OUT="$RESULTS_DIR/pangenome"
GFFS=()
while IFS= read -r gff; do
    GFFS+=("$gff")
done < <(find "$GFF_ROOT" -type f -name '*.gff' -print | LC_ALL=C sort)
[[ ${#GFFS[@]} -gt 1 ]] || { echo "Fewer than two Prokka GFF files found under $GFF_ROOT" >&2; exit 1; }

if [[ -d "$OUT" && -n "$(find "$OUT" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]; then
    if [[ "$FORCE" -eq 1 ]]; then
        rm -rf "$OUT"
    else
        echo "Output directory is not empty; use --force: $OUT" >&2
        exit 1
    fi
fi
mkdir -p "$OUT" "$RESULTS_DIR/provenance"

panaroo \
    --input "${GFFS[@]}" \
    --out_dir "$OUT" \
    --clean-mode "$PANAROO_CLEAN_MODE" \
    --remove-invalid-genes \
    --alignment core \
    --aligner mafft \
    --core_threshold "$CORE_THRESHOLD" \
    --threads "$THREADS"

CORE_ALIGNMENT="$OUT/core_gene_alignment_filtered.aln"
[[ -s "$CORE_ALIGNMENT" ]] || {
    echo "Expected Panaroo core alignment not found: $CORE_ALIGNMENT" >&2
    exit 1
}
snp-sites -o "$OUT/core_snps.aln" "$CORE_ALIGNMENT"
"$IQTREE" \
    -s "$OUT/core_snps.aln" \
    -m "$IQTREE_MODEL" \
    -T "$THREADS" \
    -seed "$SEED" \
    -pre "$OUT/core_snps"

{
    date -u '+run_completed_utc\t%Y-%m-%dT%H:%M:%SZ'
    panaroo --version 2>&1 | head -n 1
    snp-sites -V 2>&1 | head -n 1
    "$IQTREE" -version 2>&1 | head -n 1
} > "$RESULTS_DIR/provenance/02_pangenome_phylogeny_versions.txt"
