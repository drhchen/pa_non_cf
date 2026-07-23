#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: 05_mcorr.sh --config FILE --baps1-xmfa FILE --baps2-xmfa FILE [--force]

Run mcorr-xmfa and mcorr-fit for the two lineages. Additional build-specific
arguments can be supplied with MCORR_XMFA_ARGS and MCORR_FIT_ARGS in the config.
Archive the installed executables and help/version output with the final run.
EOF
}

CONFIG=""
BAPS1_XMFA=""
BAPS2_XMFA=""
FORCE=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config) CONFIG=${2:?}; shift 2 ;;
        --baps1-xmfa) BAPS1_XMFA=${2:?}; shift 2 ;;
        --baps2-xmfa) BAPS2_XMFA=${2:?}; shift 2 ;;
        --force) FORCE=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done
[[ -f "$CONFIG" ]] || { echo "Config file not found: $CONFIG" >&2; exit 2; }
[[ -s "$BAPS1_XMFA" && -s "$BAPS2_XMFA" ]] || { usage >&2; exit 2; }
# shellcheck source=/dev/null
source "$CONFIG"
: "${RESULTS_DIR:?Set RESULTS_DIR in the config}"

for executable in mcorr-xmfa mcorr-fit; do
    command -v "$executable" >/dev/null 2>&1 || { echo "$executable not found" >&2; exit 127; }
done

OUT="$RESULTS_DIR/mcorr"
mkdir -p "$OUT" "$RESULTS_DIR/provenance"
read -r -a XMFA_ARGS <<< "${MCORR_XMFA_ARGS:-}"
read -r -a FIT_ARGS <<< "${MCORR_FIT_ARGS:-}"
for lineage in BAPS1 BAPS2; do
    input_variable="${lineage}_XMFA"
    input=${!input_variable}
    raw_prefix="$OUT/$lineage"
    fit_prefix="$OUT/$lineage"
    if [[ -e "${raw_prefix}.csv" && "$FORCE" -ne 1 ]]; then
        echo "mcorr output exists; use --force: ${raw_prefix}.csv" >&2
        exit 1
    fi
    mcorr-xmfa "${XMFA_ARGS[@]}" "$input" "$raw_prefix"
    mcorr-fit "${FIT_ARGS[@]}" "${raw_prefix}.csv" "$fit_prefix"
done

{
    date -u '+run_completed_utc\t%Y-%m-%dT%H:%M:%SZ'
    mcorr-xmfa --help 2>&1 | head -n 20 || true
    mcorr-fit --help 2>&1 | head -n 20 || true
} > "$RESULTS_DIR/provenance/05_mcorr_build.txt"
