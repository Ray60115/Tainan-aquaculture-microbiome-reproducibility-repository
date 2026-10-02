#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python}"
RESULTS="$ROOT/reproduced_results"

mkdir -p "$RESULTS/functional"

TAXONOMIC_TURNOVER="$RESULTS/taxonomy/tables/turnover_values.csv"
if [ ! -f "$TAXONOMIC_TURNOVER" ]; then
  TAXONOMIC_TURNOVER="$ROOT/reference_results/taxonomy/turnover_values.csv"
fi

"$PYTHON_BIN" "$ROOT/code/01_prepare_picrust2_inputs.py" \
  --asv "$ROOT/data/primary/asv_table_168_profiles.csv.gz" \
  --controls "$ROOT/data/primary/asv_table_9_nonformal_columns.csv.gz" \
  --metadata "$ROOT/data/primary/sample_metadata_168.csv" \
  --outdir "$ROOT"

# This command reruns the downstream dependence-aware functional statistics
# from the archived primary PICRUSt2 outputs.  A full phylogenetic-placement
# rerun is documented separately in code/03_run_picrust2_primary.sh.
"$PYTHON_BIN" "$ROOT/code/05_analyze_predicted_functions.py" \
  --project "$RESULTS/functional" \
  --metadata "$ROOT/data/primary/sample_metadata_168.csv" \
  --picrust "$ROOT/picrust2_primary" \
  --abundance-input "$ROOT/input/study_abundance_168_prokaryotic.tsv" \
  --taxonomy-turnover "$TAXONOMIC_TURNOVER" \
  --module-definitions "$ROOT/input/pathway_module_definitions.csv" \
  --pathway-name-map "$ROOT/reference_results/functional/Table_S11_pathway_module_definitions.csv"

echo "Functional analysis complete: $RESULTS/functional"
