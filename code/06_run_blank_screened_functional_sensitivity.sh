#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PRIMARY="${ANALYSIS_ROOT}/picrust2_primary"
SCREENED="${BLANK_SCREENED_OUTPUT_DIR:-${ANALYSIS_ROOT}/reproduced_results/picrust2_blank_screened}"
INPUT="${ANALYSIS_ROOT}/input/study_abundance_168_prokaryotic_blank_screened.tsv"
LOG_DIR="${ANALYSIS_ROOT}/reproduced_results/logs"
PROCESSES="${PICRUST2_PROCESSES:-1}"

mkdir -p "${LOG_DIR}"

for required_command in metagenome_pipeline.py pathway_pipeline.py; do
  if ! command -v "${required_command}" >/dev/null 2>&1; then
    echo "Required PICRUSt2 command not found on PATH: ${required_command}" >&2
    exit 2
  fi
done
for required_file in \
  "${INPUT}" \
  "${PRIMARY}/combined_marker_predicted_and_nsti.tsv.gz" \
  "${PRIMARY}/combined_EC_predicted.tsv.gz"; do
  if [ ! -f "${required_file}" ]; then
    echo "Required file not found: ${required_file}" >&2
    echo "Verify that the v1.0.2 archive is complete, or rerun the primary PICRUSt2 prediction step to recreate the missing file." >&2
    exit 2
  fi
done

if [ -e "${SCREENED}" ]; then
  echo "Sensitivity output already exists: ${SCREENED}" >&2
  exit 2
fi

mkdir -p "${SCREENED}"
ln -s "${PRIMARY}/combined_marker_predicted_and_nsti.tsv.gz" \
  "${SCREENED}/combined_marker_predicted_and_nsti.tsv.gz"
ln -s "${PRIMARY}/combined_EC_predicted.tsv.gz" \
  "${SCREENED}/combined_EC_predicted.tsv.gz"

PYTHONUNBUFFERED=1 metagenome_pipeline.py \
  --input "${INPUT}" \
  --function "${PRIMARY}/combined_EC_predicted.tsv.gz" \
  --marker "${PRIMARY}/combined_marker_predicted_and_nsti.tsv.gz" \
  --max_nsti 2 \
  --out_dir "${SCREENED}/EC_metagenome_out" \
  > "${LOG_DIR}/picrust2_blank_screened_metagenome.stdout.log" \
  2> "${LOG_DIR}/picrust2_blank_screened_metagenome.stderr.log"

PYTHONUNBUFFERED=1 pathway_pipeline.py \
  --input "${SCREENED}/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz" \
  --out_dir "${SCREENED}/pathways_out" \
  --processes "${PROCESSES}" \
  --intermediate "${SCREENED}/intermediate_pathways" \
  --verbose \
  > "${LOG_DIR}/picrust2_blank_screened_pathways.stdout.log" \
  2> "${LOG_DIR}/picrust2_blank_screened_pathways.stderr.log"

echo "Blank-screened functional sensitivity run completed."
