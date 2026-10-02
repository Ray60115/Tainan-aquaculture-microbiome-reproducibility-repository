# Code guide

Use `../run_core_analysis.sh` and `../run_functional_analysis.sh` for the
integrated execution order. The principal active scripts are:

- `01_build_sample_manifest.py`
- `01_attach_sra_accessions.py`
- `02_reconstruct_asv_filter.py`
- `03_run_rarefaction_sensitivity.py`
- `02_reanalyze_prokaryotic_taxonomy.py` with
  `reanalyze_dependence_aware.py`
- `10_run_prokaryotic_random_forest.py` with
  `run_random_forest_reconstruction.py`
- `11_run_multialgorithm_reconstruction.py`
- `12_rebuild_prokaryotic_figures.py`
- `13_plot_water_quality.py`
- `14_plot_current_seasonal_prevalence.py`
- `00_validate_manuscript_invariants.py`
- `04_build_complete_lineage_table.py`
- `01_prepare_picrust2_inputs.py`, `03_run_picrust2_primary.sh`, and
  `05_analyze_predicted_functions.py`

Compatibility scripts for the original PICRUSt2/SEPP/EPA-ng run are retained
because they document the actual constrained-runtime solution.

`06_run_blank_screened_functional_sensitivity.sh` and
`07_compare_blank_screened_functional_results.py` implement the negative-control
sensitivity path. The primary combined EC prediction and compact rerun outputs
are archived, and `../run_blank_screened_functional_sensitivity.sh` executes the
full preparation, PICRUSt2, downstream-statistics, and comparison sequence.

`prepare_complete_climate_2024.py` is a provenance helper for appending the
November–December 2024 station files. It now takes explicit command-line paths;
the original raw station downloads are external and are not redistributed,
whereas the finalized daily climate table used by every analysis is included
under `../data/primary/`.
