#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
INPUT_FASTA="${ANALYSIS_ROOT}/input/study_seqs_168_prokaryotic.fna"
INPUT_TABLE="${ANALYSIS_ROOT}/input/study_abundance_168_prokaryotic.tsv"
OUTPUT_DIR="${PICRUST_OUTPUT_DIR:-${ANALYSIS_ROOT}/reproduced_results/picrust2_primary_rerun}"
ENVIRONMENT_DIR="${ANALYSIS_ROOT}/environment"
LOG_DIR="${ANALYSIS_ROOT}/reproduced_results/logs"

if [ -z "${BAC_EC_TABLE:-}" ] || [ -z "${ARC_EC_TABLE:-}" ]; then
  echo "Set BAC_EC_TABLE and ARC_EC_TABLE to the installed official PICRUSt2-SC EC reference tables." >&2
  exit 2
fi
for required_file in "${INPUT_FASTA}" "${INPUT_TABLE}" "${BAC_EC_TABLE}" "${ARC_EC_TABLE}"; do
  if [ ! -f "${required_file}" ]; then
    echo "Required file not found: ${required_file}" >&2
    exit 2
  fi
done
for required_command in picrust2_pipeline.py run_sepp.py epa-ng gappa hmmsearch; do
  if ! command -v "${required_command}" >/dev/null 2>&1; then
    echo "Required command not found on PATH: ${required_command}" >&2
    exit 2
  fi
done

mkdir -p "${ENVIRONMENT_DIR}" "${LOG_DIR}" "$(dirname "${OUTPUT_DIR}")"

if [ -n "${SEPP_CONFIG:-}" ]; then
  python "${SCRIPT_DIR}/configure_sepp_compat.py" \
    --config "${SEPP_CONFIG}" \
    --adapter "${SCRIPT_DIR}/pplacer_epa_ng_compat.py" \
    --record "${ENVIRONMENT_DIR}/sepp_main_config_rerun.ini"
fi

picrust2_pipeline.py --version > "${ENVIRONMENT_DIR}/picrust2_version_rerun.txt"

if [ -d "${OUTPUT_DIR}" ]; then
  if find "${OUTPUT_DIR}" -mindepth 1 -maxdepth 1 | read -r; then
    echo "Output directory is not empty: ${OUTPUT_DIR}" >&2
    exit 2
  fi
  rmdir "${OUTPUT_DIR}"
fi

{
  echo "PICRUSt2 version: $(picrust2_pipeline.py --version)"
  echo "SEPP version: $(run_sepp.py --version 2>&1 | tail -n 1)"
  echo "EPA-ng version: $(epa-ng --version 2>&1 | head -n 1)"
  echo "gappa version: $(gappa --version 2>&1 | head -n 1)"
  echo "HMMER version: $(hmmsearch -h 2>&1 | sed -n '2p')"
  echo "Compatibility adapter SHA-256: $(sha256sum "${SCRIPT_DIR}/pplacer_epa_ng_compat.py" | cut -d' ' -f1)"
  echo "Input FASTA: ${INPUT_FASTA}"
  echo "Input table: ${INPUT_TABLE}"
  echo "Output directory: ${OUTPUT_DIR}"
  echo "Processes: 6"
  echo "Predicted gene-family set: EC (the gene-family input used for MetaCyc pathway inference; optional KO output was not requested)."
  echo "EC references: the installed official PICRUSt2-SC bacterial and archaeal EC tables were passed explicitly as custom trait tables."
  echo "Implementation note: explicit EC-table arguments avoid a PICRUSt2 2.6.3 preflight routine that otherwise loads every installed default trait category, including unused GO/PFAM/KO tables, before respecting --in_traits."
  echo "Placement: SEPP hierarchical decomposition and HMM alignment with six processes; likelihood placement within each SEPP reference subtree was performed by EPA-ng 0.3.8 through scripts/pplacer_epa_ng_compat.py, with EPA-ng jobs serialized by a memory lock, followed by SEPP's standard jplace merger."
  echo "Compatibility rationale: the pplacer executable bundled with SEPP 4.5.5 segfaulted before reading data on the current Linux 6.12 runtime. Direct whole-reference EPA-ng placement exceeded the 14-GiB cgroup limit, whereas EPA-ng placement on SEPP's decomposed reference subtrees completed within the available resources."
  echo "Defaults retained: bacterial and archaeal PICRUSt2-SC references; marker-copy normalization; max NSTI 2; MinPath and pathway gap filling."
  echo "Resource note: initial default-table attempts were terminated during preflight because PICRUSt2 loaded unused trait tables. The primary run therefore passes only the official EC tables required for MetaCyc inference."
  echo "Command:"
  echo "picrust2_pipeline.py -s ${INPUT_FASTA} -i ${INPUT_TABLE} -o ${OUTPUT_DIR} -p 6 -t sepp --custom_trait_tables_ref1 ${BAC_EC_TABLE} --custom_trait_tables_ref2 ${ARC_EC_TABLE} --verbose"
} > "${ENVIRONMENT_DIR}/picrust2_run_manifest_rerun.txt"

SECONDS=0
PYTHONUNBUFFERED=1 picrust2_pipeline.py \
  -s "${INPUT_FASTA}" \
  -i "${INPUT_TABLE}" \
  -o "${OUTPUT_DIR}" \
  -p 6 \
  -t sepp \
  --custom_trait_tables_ref1 "${BAC_EC_TABLE}" \
  --custom_trait_tables_ref2 "${ARC_EC_TABLE}" \
  --verbose \
  > "${LOG_DIR}/picrust2_primary.stdout.log" \
  2> "${LOG_DIR}/picrust2_primary.stderr.log"

echo "Elapsed_seconds=${SECONDS}" \
  >> "${ENVIRONMENT_DIR}/picrust2_run_manifest_rerun.txt"
echo "PICRUSt2 primary run completed."
