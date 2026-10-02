# Manuscript-to-code map

| Manuscript component | Primary code | Principal input/output |
|---|---|---|
| HiFi read QC, primer trimming, DADA2, chimera removal, taxonomy | `workflow/HiFi-16S-workflow_v0.9_commit_d29598f/` | 204 FASTQ inputs → workflow ASV/QC outputs |
| 204 → 190 → 168 sample mapping | `code/01_build_sample_manifest.py` | `workflow/qc/` and `data/primary/library_to_profile_metadata_190.csv` → `data/manifests/` |
| 29,461 → 23,896 ASV filtering and exact profile averaging | `code/02_reconstruct_asv_filter.py` | workflow and averaged ASV tables → ASV audit |
| Library depth and rarefaction sensitivity | `code/03_run_rarefaction_sensitivity.py` | 190-library feature table → rarefaction sample/aggregate tables |
| Shannon, PCoA, seasonal/layer/event tests, PERMANOVA, PERMDISP, turnover, vectors, blockwise R² | `code/02_reanalyze_prokaryotic_taxonomy.py`; `code/reanalyze_dependence_aware.py` | `reproduced_results/taxonomy/` |
| Random-forest nested reconstruction and grouped importance | `code/10_run_prokaryotic_random_forest.py`; `code/run_random_forest_reconstruction.py` | Table 3, Fig. 11, Table S13 |
| Six-algorithm and grow-out-only sensitivity | `code/11_run_multialgorithm_reconstruction.py` | Table S12a,b and grow-out rows of Table 3 |
| PICRUSt2 input preparation and primary inference | `code/01_prepare_picrust2_inputs.py`; `code/03_run_picrust2_primary.sh`; compatibility helpers | 23,896 ASVs → archived primary PICRUSt2 outputs |
| Functional PERMANOVA/PERMDISP, modules, turnover and plots | `code/05_analyze_predicted_functions.py` | Table 4, Tables S10–S11, Fig. 10 |
| Negative-control functional sensitivity | `run_blank_screened_functional_sensitivity.sh`; `code/06_run_blank_screened_functional_sensitivity.sh`; `code/07_compare_blank_screened_functional_results.py` | `picrust2_blank_screened/`; `reference_results/functional_blank_screened/`; `reference_results/functional_sensitivity/` |
| Complete resolved-genus lineage table | `code/04_build_complete_lineage_table.py` | Additional file 1 and its summary JSON |

## Current main figures

| Figure | Rebuilding code | Note |
|---|---|---|
| Fig. 1 | `code/12_rebuild_prokaryotic_figures.py` → `Fig1_meteorological_context_panels_b_to_f.png`; upstream site-map asset | Climate panels are code-derived; exact map panel awaits its source asset. |
| Fig. 2 | `code/12_rebuild_prokaryotic_figures.py` → `Fig2_community_and_alpha_diversity.png` | Community composition and Shannon diversity. |
| Fig. 3 | `code/12_rebuild_prokaryotic_figures.py` → `Fig3_taxonomic_descriptive_contrasts.png` | Descriptive nursery/grow-out and pond contrasts. |
| Fig. 4 | `code/14_plot_current_seasonal_prevalence.py` | Updated prevalence, not the older abundance heatmap. |
| Fig. 5 | `code/12_rebuild_prokaryotic_figures.py` → `Fig5_turnover_dependence_aware.png` | Final annotation placement may be adjusted at typesetting. |
| Fig. 6 | `code/13_plot_water_quality.py` | Pond-month water-quality series. |
| Fig. 7 | `code/12_rebuild_prokaryotic_figures.py` → `Fig7_environmental_vectors_dependence_aware.png` | Vector statistics are unchanged; final long-label placement may be adjusted manually. |
| Fig. 8 | `code/12_rebuild_prokaryotic_figures.py` → `Fig8_disturbance_dependence_aware.png` | Individual disturbance tests. |
| Fig. 9 | `code/12_rebuild_prokaryotic_figures.py` → `Fig9_dependence_aware_variance_partition.png` | Blockwise variation partitioning. |
| Fig. 10 | `code/05_analyze_predicted_functions.py` | Predicted functional potential. |
| Fig. 11 | `code/10_run_prokaryotic_random_forest.py`; `code/11_make_merged_rf_figure.py` | Nested reconstruction and grouped importance. |

The archived plotting functions preserve the analysis-generated values. Minor
label repositioning, transparent backgrounds, and final panel numbering are
presentation-only edits and do not alter statistical output.
