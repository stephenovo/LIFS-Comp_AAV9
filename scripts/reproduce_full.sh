#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: bash scripts/reproduce_full.sh \
  --fit4function-dir PATH --manifest PATH --data-root PATH \
  --reconstructed PATH --screen-csv PATH [options]

Runs the data-dependent audit, benchmarks, and virtual screen without
downloading or modifying external research data.

Options:
  --output-dir PATH       Artifact directory (default: artifacts/full-data)
  --pool-size N           Virtual-screen pool size (default: 1000000)
  --ensemble-size N       MLP ensemble size (default: 5)
  --max-iter N             MLP maximum iterations (default: 80)
EOF
}

FIT4FUNCTION_DIR=""
MANIFEST=""
DATA_ROOT=""
RECONSTRUCTED=""
SCREEN_CSV=""
OUTPUT_DIR="artifacts/full-data"
POOL_SIZE=1000000
ENSEMBLE_SIZE=5
MAX_ITER=80

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fit4function-dir) FIT4FUNCTION_DIR="$2"; shift 2 ;;
    --manifest) MANIFEST="$2"; shift 2 ;;
    --data-root) DATA_ROOT="$2"; shift 2 ;;
    --reconstructed) RECONSTRUCTED="$2"; shift 2 ;;
    --screen-csv) SCREEN_CSV="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --pool-size) POOL_SIZE="$2"; shift 2 ;;
    --ensemble-size) ENSEMBLE_SIZE="$2"; shift 2 ;;
    --max-iter) MAX_ITER="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

for required in FIT4FUNCTION_DIR MANIFEST DATA_ROOT RECONSTRUCTED SCREEN_CSV; do
  if [[ -z "${!required}" ]]; then
    echo "Missing required option for ${required}." >&2
    usage >&2
    exit 2
  fi
done

for path in "$FIT4FUNCTION_DIR" "$MANIFEST" "$RECONSTRUCTED" "$SCREEN_CSV"; do
  [[ -e "$path" ]] || { echo "Missing input: $path" >&2; exit 1; }
done

expected_commit="$(python -c 'import json,sys; print(json.load(open(sys.argv[1], encoding="utf-8"))["fit4function"]["commit"])' docs/source_manifest.json)"
actual_commit="$(git -C "$FIT4FUNCTION_DIR" rev-parse HEAD 2>/dev/null)" || {
  echo "Fit4Function directory is not a git checkout: $FIT4FUNCTION_DIR" >&2
  exit 1
}
if [[ "$actual_commit" != "$expected_commit" ]]; then
  echo "Fit4Function commit mismatch: expected $expected_commit, found $actual_commit" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"
echo "[1/5] Verifying external-data checksums"
aav9-sma verify-manifest "$MANIFEST" --root "$DATA_ROOT" > "$OUTPUT_DIR/manifest_verification.json"

echo "[2/5] Auditing pinned Fit4Function release"
aav9-sma audit-fit4function "$FIT4FUNCTION_DIR" \
  --output "$OUTPUT_DIR/fit4function_release_audit.json"

echo "[3/5] Running frozen animal-4 development benchmarks"
aav9-sma benchmark-multiorgan "$RECONSTRUCTED" \
  --models ridge random_forest \
  --output "$OUTPUT_DIR/multiorgan_baseline_metrics.csv"
aav9-sma benchmark-multitask-ensemble "$RECONSTRUCTED" \
  --ensemble-size "$ENSEMBLE_SIZE" --max-iter "$MAX_ITER" \
  --output "$OUTPUT_DIR/multitask_ensemble_metrics.csv"

echo "[4/5] Running virtual screen"
aav9-sma screen-virtual "$SCREEN_CSV" "$RECONSTRUCTED" \
  --pool-size "$POOL_SIZE" --ensemble-size "$ENSEMBLE_SIZE" --max-iter "$MAX_ITER" \
  --output-ranked "$OUTPUT_DIR/virtual_screen_ranked.csv.gz" \
  --output-pareto "$OUTPUT_DIR/virtual_screen_pareto.csv" \
  --output-shortlist "$OUTPUT_DIR/virtual_screen_shortlist.csv" \
  --output-summary "$OUTPUT_DIR/virtual_screen_summary.json"

echo "[5/5] Auditing gate sensitivity and residue composition"
aav9-sma audit-screen-funnel \
  "$OUTPUT_DIR/virtual_screen_ranked.csv.gz" \
  "$OUTPUT_DIR/virtual_screen_shortlist.csv" \
  "$SCREEN_CSV" \
  --output-gates "$OUTPUT_DIR/packaging_gate_sensitivity.csv" \
  --output-composition "$OUTPUT_DIR/funnel_composition_audit.csv" \
  --output-summary "$OUTPUT_DIR/funnel_sensitivity_summary.json"

echo "Full-data run completed in $OUTPUT_DIR"
