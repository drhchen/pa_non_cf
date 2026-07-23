#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: 01_qc_annotation.sh --config FILE [--force]

Validate assemblies, create normalized uncompressed FASTA inputs, run QUAST,
CheckM, and FastANI, apply publication thresholds, then annotate retained
assemblies with Prokka.
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

: "${ASSEMBLY_DIR:?Set ASSEMBLY_DIR in the config}"
: "${METADATA:?Set METADATA in the config}"
: "${RESULTS_DIR:?Set RESULTS_DIR in the config}"
: "${THREADS:?Set THREADS in the config}"

for executable in python3 prokka fastANI checkm; do
    command -v "$executable" >/dev/null 2>&1 || {
        echo "Required executable not found: $executable" >&2
        exit 127
    }
done
if command -v quast.py >/dev/null 2>&1; then
    QUAST=quast.py
elif command -v quast >/dev/null 2>&1; then
    QUAST=quast
else
    echo "Required executable not found: quast.py or quast" >&2
    exit 127
fi

OUT="$RESULTS_DIR/qc_annotation"
MANIFEST="$RESULTS_DIR/provenance/assemblies.tsv"
mkdir -p "$RESULTS_DIR/provenance" "$OUT/assemblies" "$OUT/quast" "$OUT/prokka" "$OUT/fastani"

bash scripts/bash/00_validate_inputs.sh \
    --assemblies "$ASSEMBLY_DIR" \
    --metadata "$METADATA" \
    --manifest "$MANIFEST"

while IFS=$'\t' read -r sample assembly; do
    [[ "$sample" == "sample" ]] && continue
    normalized="$OUT/assemblies/$sample.fna"
    if [[ -e "$normalized" && "$FORCE" -ne 1 ]]; then
        echo "Output exists; use --force to replace generated files: $normalized" >&2
        exit 1
    fi
    if [[ "$assembly" == *.gz ]]; then
        gzip -dc "$assembly" > "$normalized"
    else
        cp "$assembly" "$normalized"
    fi

    quast_out="$OUT/quast/$sample"
    if [[ "$FORCE" -eq 1 ]]; then
        rm -rf "$quast_out"
    fi
    "$QUAST" --threads "$THREADS" --output-dir "$quast_out" "$normalized"
done < "$MANIFEST"

find "$OUT/assemblies" -type f -name '*.fna' -print | LC_ALL=C sort > "$OUT/fastani/assemblies.list"
cp "$OUT/fastani/assemblies.list" "$OUT/fastani/query.list"
cp "$OUT/fastani/assemblies.list" "$OUT/fastani/reference.list"

fastANI \
    --ql "$OUT/fastani/query.list" \
    --rl "$OUT/fastani/reference.list" \
    --threads "$THREADS" \
    -o "$OUT/fastani/pairwise_ani.tsv"
python3 scripts/python/ani_to_matrix.py \
    --input "$OUT/fastani/pairwise_ani.tsv" \
    --output "$OUT/fastani/pairwise_ani_matrix.tsv"

[[ -s "$REFERENCE_FASTA" ]] || { echo "Reference FASTA not found: $REFERENCE_FASTA" >&2; exit 1; }
fastANI \
    --ql "$OUT/fastani/query.list" \
    --ref "$REFERENCE_FASTA" \
    --threads "$THREADS" \
    -o "$OUT/fastani/to_reference.tsv"

if [[ "$FORCE" -eq 1 ]]; then
    rm -rf "$OUT/checkm/workflow"
fi
mkdir -p "$OUT/checkm"
checkm lineage_wf \
    --extension fna \
    --threads "$THREADS" \
    --tab_table \
    --file "$OUT/checkm/lineage_wf.tsv" \
    "$OUT/assemblies" \
    "$OUT/checkm/workflow"

RETAINED="$RESULTS_DIR/provenance/retained_assemblies.tsv"
python3 scripts/python/filter_genomes.py \
    --manifest "$MANIFEST" \
    --checkm "$OUT/checkm/lineage_wf.tsv" \
    --ani "$OUT/fastani/pairwise_ani.tsv" \
    --retained "$RETAINED" \
    --excluded "$RESULTS_DIR/provenance/excluded_assemblies.tsv" \
    --min-completeness "$MIN_COMPLETENESS" \
    --max-contamination "$MAX_CONTAMINATION" \
    --min-ani "$MIN_ANI"
python3 scripts/python/subset_metadata.py \
    --metadata "$METADATA" \
    --manifest "$RETAINED" \
    --output "$RESULTS_DIR/provenance/retained_metadata.tsv"

index=1
while IFS=$'\t' read -r sample _assembly; do
    [[ "$sample" == "sample" ]] && continue
    normalized="$OUT/assemblies/$sample.fna"
    prokka_out="$OUT/prokka/$sample"
    if [[ "$FORCE" -eq 1 ]]; then
        rm -rf "$prokka_out"
    fi
    locus_tag=$(printf 'PA%06d' "$index")
    prokka \
        --outdir "$prokka_out" \
        --prefix "$sample" \
        --locustag "$locus_tag" \
        --genus Pseudomonas \
        --species aeruginosa \
        --cpus "$THREADS" \
        --compliant \
        "$normalized"
    index=$((index + 1))
done < "$RETAINED"

{
    date -u '+run_completed_utc\t%Y-%m-%dT%H:%M:%SZ'
    prokka --version 2>&1 | head -n 1
    fastANI --version 2>&1 | head -n 1
    checkm tree --help 2>&1 | head -n 1 || true
    "$QUAST" --version 2>&1 | head -n 1
} > "$RESULTS_DIR/provenance/01_qc_annotation_versions.txt"
