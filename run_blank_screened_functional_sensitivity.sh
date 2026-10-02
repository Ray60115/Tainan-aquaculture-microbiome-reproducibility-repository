#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python}"
RESULTS="${ROOT}/reproduced_results"

"${PYTHON_BIN}" "${ROOT}/code/01_prepare_picrust2_inputs.py" \
  --asv "${ROOT}/data/primary/asv_table_168_profiles.csv.gz" \
  --controls "${ROOT}/data/primary/asv_table_9_nonformal_columns.csv.gz" \
  --metadata "${ROOT}/data/primary/sample_metadata_168.csv" \
  --outdir "${ROOT}"

PICRUST2_PROCESSES="${PICRUST2_PROCESSES:-1}" \
  bash "${ROOT}/code/06_run_blank_screened_functional_sensitivity.sh"

TAXONOMIC_TURNOVER="${RESULTS}/taxonomy/tables/turnover_values.csv"
if [ ! -f "${TAXONOMIC_TURNOVER}" ]; then
  TAXONOMIC_TURNOVER="${ROOT}/reference_results/taxonomy/turnover_values.csv"
fi

for cohort in primary blank_screened; do
  if [ "${cohort}" = "primary" ]; then
    PICRUST_DIR="${ROOT}/picrust2_primary"
    ABUNDANCE="${ROOT}/input/study_abundance_168_prokaryotic.tsv"
    PROJECT="${RESULTS}/functional"
  else
    PICRUST_DIR="${RESULTS}/picrust2_blank_screened"
    ABUNDANCE="${ROOT}/input/study_abundance_168_prokaryotic_blank_screened.tsv"
    PROJECT="${RESULTS}/functional_blank_screened"
  fi

  "${PYTHON_BIN}" "${ROOT}/code/05_analyze_predicted_functions.py" \
    --project "${PROJECT}" \
    --metadata "${ROOT}/data/primary/sample_metadata_168.csv" \
    --picrust "${PICRUST_DIR}" \
    --abundance-input "${ABUNDANCE}" \
    --taxonomy-turnover "${TAXONOMIC_TURNOVER}" \
    --module-definitions "${ROOT}/input/pathway_module_definitions.csv" \
    --pathway-name-map \
      "${ROOT}/reference_results/functional/Table_S11_pathway_module_definitions.csv"
done

"${PYTHON_BIN}" "${ROOT}/code/07_compare_blank_screened_functional_results.py" \
  --project "${RESULTS}/functional_sensitivity_comparison" \
  --primary-picrust "${ROOT}/picrust2_primary" \
  --screened-picrust "${RESULTS}/picrust2_blank_screened" \
  --primary-results "${RESULTS}/functional" \
  --screened-results "${RESULTS}/functional_blank_screened"

echo "Blank-screened functional sensitivity analysis complete: ${RESULTS}"
