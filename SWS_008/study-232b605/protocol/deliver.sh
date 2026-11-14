#!/usr/bin/env bash
# Run the protocol on a samplesheet and write the delivery into OUT.
#
# usage: deliver.sh --samplesheet CSV --annotations DIR --out OUT
#                   [--settings JSON] [--assay insitu|capture] [--baitmap FILE] [--keep-rejected]
#
# Environment: NEXTFLOW (workflow engine launcher; nextflow if unset), PYTHON (python3 if unset).
# The analysis package `hicparams` must be importable by PYTHON.
set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NEXTFLOW="${NEXTFLOW:-nextflow}"
PYTHON="${PYTHON:-python3}"
SETTINGS="$HERE/assay_settings.json"
ASSAY="insitu"
SHEET=""
ANN=""
OUT=""
BAITMAP=""
KEEP="false"

while [ $# -gt 0 ]; do
  case "$1" in
    --samplesheet) SHEET="$2"; shift 2 ;;
    --annotations) ANN="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --settings) SETTINGS="$2"; shift 2 ;;
    --assay) ASSAY="$2"; shift 2 ;;
    --baitmap) BAITMAP="$2"; shift 2 ;;
    --keep-rejected) KEEP="true"; shift ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done
if [ -z "$SHEET" ] || [ -z "$ANN" ] || [ -z "$OUT" ]; then
  echo "usage: deliver.sh --samplesheet CSV --annotations DIR --out OUT [--settings JSON] [--assay insitu|capture] [--baitmap FILE] [--keep-rejected]" >&2
  exit 2
fi

mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"
SHEET="$(cd "$(dirname "$SHEET")" && pwd)/$(basename "$SHEET")"
ANN="$(cd "$ANN" && pwd)"
if [ -n "$BAITMAP" ]; then
  BAITMAP="$(cd "$(dirname "$BAITMAP")" && pwd)/$(basename "$BAITMAP")"
fi

cd "$OUT" || exit 1
ARGS=(run "$HERE/main.nf" -c "$HERE/nextflow.config" -ansi-log false
      -output-dir "$OUT/results" -work-dir "$OUT/work"
      --samplesheet "$SHEET" --settings "$SETTINGS" --annotations "$ANN"
      --assay "$ASSAY" --keep_rejected "$KEEP")
if [ -n "$BAITMAP" ]; then
  ARGS+=(--baitmap "$BAITMAP")
fi

"$NEXTFLOW" "${ARGS[@]}"
RC_RUN=$?
rm -rf "$OUT/work"

"$PYTHON" -m hicparams "$OUT/results/pairs.csv" filtered "$OUT/parameters.tsv" --annotations "$ANN"
RC_PARAMS=$?
"$PYTHON" "$HERE/bin/qc_table.py" "$OUT/results/qc.csv" "$OUT/qc.tsv"
RC_QC=$?
cp "$OUT/.nextflow.log" "$OUT/run.log" 2>/dev/null

for rc in "$RC_RUN" "$RC_PARAMS" "$RC_QC"; do
  if [ "$rc" -ne 0 ]; then
    exit "$rc"
  fi
done
exit 0
