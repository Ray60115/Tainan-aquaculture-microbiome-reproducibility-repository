#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python}"
RESULTS="$ROOT/reproduced_results"

mkdir -p "$RESULTS/manifests" "$RESULTS/rarefaction" "$RESULTS/taxonomy" "$RESULTS/model_sensitivity"

"$PYTHON_BIN" - <<'PY'
import importlib
required = ["numpy", "pandas", "scipy", "sklearn", "matplotlib", "seaborn", "xgboost", "pygam"]
missing = []
for package in required:
    try:
        importlib.import_module(package)
    except ImportError:
        missing.append(package)
if missing:
    raise SystemExit(
        "Missing required Python packages: " + ", ".join(missing)
        + ". Install requirements.txt before running the integrated workflow."
    )
PY

"$PYTHON_BIN" "$ROOT/code/01_build_sample_manifest.py" \
  --workflow-manifest "$ROOT/workflow/sample_input_204.tsv" \
  --sample-metadata "$ROOT/data/primary/library_to_profile_metadata_190.csv" \
  --dada2-qc "$ROOT/workflow/qc/dada2_qc.tsv" \
  --read-qc "$ROOT/workflow/qc/read_qc_summary.tsv" \
  --trim-qc "$ROOT/workflow/qc/trim_qc_summary.tsv" \
  --output-dir "$RESULTS/manifests"

if [ -f "$ROOT/data/accessions/NCBI_SRA_metadata_190.tsv" ]; then
  "$PYTHON_BIN" "$ROOT/code/01_attach_sra_accessions.py" \
    --sra-metadata "$ROOT/data/accessions/NCBI_SRA_metadata_190.tsv" \
    --library-metadata "$ROOT/data/primary/library_to_profile_metadata_190.csv" \
    --manifest-dir "$RESULTS/manifests" \
    --output-metadata "$RESULTS/manifests/NCBI_SRA_metadata_190.tsv"
fi

"$PYTHON_BIN" "$ROOT/code/02_reconstruct_asv_filter.py" \
  --raw-asv "$ROOT/data/primary/asv_table_workflow_feature_table_199_columns.csv.gz" \
  --averaged-asv "$ROOT/data/primary/asv_table_168_profiles.csv.gz" \
  --sample-metadata "$ROOT/data/primary/library_to_profile_metadata_190.csv" \
  --output-dir "$RESULTS/manifests"

"$PYTHON_BIN" "$ROOT/code/03_run_rarefaction_sensitivity.py" \
  --library-asv "$ROOT/data/primary/asv_table_workflow_feature_table_199_columns.csv.gz" \
  --metadata "$ROOT/data/primary/library_to_profile_metadata_190.csv" \
  --output-dir "$RESULTS/rarefaction"

"$PYTHON_BIN" "$ROOT/code/02_reanalyze_prokaryotic_taxonomy.py" \
  --project-root "$ROOT/code" \
  --outdir "$RESULTS/taxonomy"

"$PYTHON_BIN" "$ROOT/code/04_build_complete_lineage_table.py" \
  --asv "$ROOT/data/primary/asv_table_168_profiles.csv.gz" \
  --metadata "$ROOT/data/primary/sample_metadata_168.csv" \
  --output "$RESULTS/taxonomy/tables/Additional_file_1_complete_genus_lineages.csv" \
  --summary "$RESULTS/taxonomy/tables/complete_genus_lineage_summary.json"

"$PYTHON_BIN" "$ROOT/code/09_recompute_prokaryotic_taxon_tables.py" \
  --asv "$ROOT/data/primary/asv_table_168_profiles.csv.gz" \
  --metadata "$ROOT/data/primary/sample_metadata_168.csv" \
  --display-definitions "$ROOT/input/taxonomic_display_definitions.csv" \
  --outdir "$RESULTS/taxonomy/tables"

"$PYTHON_BIN" "$ROOT/code/10_run_prokaryotic_random_forest.py" \
  --helper-dir "$ROOT/code" \
  --results "$RESULTS/taxonomy" \
  --asv "$ROOT/data/primary/asv_table_168_profiles.csv.gz"

"$PYTHON_BIN" "$ROOT/code/11_run_multialgorithm_reconstruction.py" \
  --asv "$ROOT/data/primary/asv_table_168_profiles.csv.gz" \
  --metadata "$RESULTS/taxonomy/tables/validated_sample_metadata_with_shannon.csv" \
  --output-dir "$RESULTS/model_sensitivity"

"$PYTHON_BIN" "$ROOT/code/12_rebuild_prokaryotic_figures.py" \
  --helper-dir "$ROOT/code" \
  --results "$RESULTS/taxonomy" \
  --asv "$ROOT/data/primary/asv_table_168_profiles.csv.gz"

"$PYTHON_BIN" "$ROOT/code/11_make_merged_rf_figure.py" \
  --results "$RESULTS/taxonomy" \
  --output "$RESULTS/taxonomy/figures/Fig11_grouped_reconstruction_and_importance.png"

"$PYTHON_BIN" "$ROOT/code/13_plot_water_quality.py" \
  --metadata "$ROOT/data/primary/sample_metadata_168.csv" \
  --output "$RESULTS/taxonomy/figures/Fig6_water_quality_current.png"

"$PYTHON_BIN" "$ROOT/code/14_plot_current_seasonal_prevalence.py" \
  --asv "$ROOT/data/primary/asv_table_168_profiles.csv.gz" \
  --metadata "$ROOT/data/primary/sample_metadata_168.csv" \
  --permanova "$RESULTS/taxonomy/tables/permanova_dependence_aware_results.csv" \
  --output "$RESULTS/taxonomy/figures/Fig4_seasonal_prevalence_current.png" \
  --selected-table "$RESULTS/taxonomy/tables/seasonal_prevalence_top30.csv"

echo "Core analysis complete: $RESULTS"
