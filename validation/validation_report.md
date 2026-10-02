# Validation report (2026-10-01)

## Static checks

- All Python files compile successfully under Python 3.12.13.
- All shell scripts pass `bash -n`.
- The public downstream requirements resolve in a clean environment. SciPy is
  pinned to 1.16.3 because pyGAM 0.10.1 declares `scipy>=1.11.1,<1.17`;
  the earlier 1.17.0 pin was not installable by a standard dependency resolver.
- The integrated core workflow was run successfully in that clean environment.
  All dependence-aware taxonomy tables matched the archived reference tables
  exactly; random-forest performance/predictions matched to floating-point
  precision. Publication-facing RF labels now use stable legacy seed keys, so
  changing `Temperature` to `Water temperature` does not alter permutations.
- The downstream primary functional workflow was rerun from the archived
  PICRUSt2 outputs. Every functional CSV matched its archived reference result
  exactly.
- The archived upstream PacBio `main.nf` is byte-for-byte identical to the
  study's used `main.nf` snapshot (SHA-256
  `d90d5390904bcc26fe8e5223006abe7f66a8ff30b9de200bcebcbc4e1c680696`).

## Sample and ASV audit

The integrated scripts were executed from the repository inputs and confirmed:

- 204 workflow inputs;
- 199 sample columns in the exported workflow feature table; the five workflow
  inputs absent from that table are exactly Probiotic1–3 and SeaWater1–2, all
  outside the formal cohort;
- 190 formal biological libraries;
- 168 canonical pond-month-layer profiles;
- 146 profiles represented by one source library and 22 by two source
  libraries;
- exact reconstruction of the supplied 168-profile ASV matrix from the 190
  biological libraries (maximum absolute cell difference = 0; zero
  nonmatching cells);
- 29,461 starting ASVs;
- 4,242 ASVs absent from all 168 biological profiles;
- 728 mitochondrial, 590 chloroplast, and 5 explicitly eukaryotic positive
  ASVs removed;
- 23,896 retained prokaryotic ASVs; and
- 15 archaeal ASVs representing 0.007323% of retained reads.

The authoritative machine-readable audit outputs are retained under
`data/manifests/`.

## Repeated rarefaction

The deterministic recovered implementation reproduced every cohort size in
submitted Table S2c: 190/168/75 at 1,500 reads through 110/91/32 at 25,000
reads (libraries/profiles/pairs). All 100 iterations at every depth returned
`p < 0.05`, and the bottom-minus-surface effect remained positive at every
depth. Because the original code archive did not contain the rarefaction
script or exact RNG stream, medians and Wilcoxon statistics differ slightly at
some depths. Submitted values and deterministic rerun values are both retained.

## Six-model reconstruction sensitivity

The newly consolidated script was run with the fixed hyperparameters reported
in Supplemental Table S12. For the primary R² and MAE columns:

- 46 of 50 model/predictor/cohort rows match the submitted table at four
  decimal places;
- ridge regression, random forest, gradient boosting, XGBoost, and SVR match
  in both cohorts;
- GAM structure-only and season models also match; and
- the four GAM rows that add water-quality or disturbance smooths have small
  numerical differences (maximum absolute difference: R² 0.00472; MAE
  0.00125). The direction, predictor-block ordering, and manuscript
  interpretation are unchanged.

Submitted and rerun files are stored side by side under
`reference_results/sensitivity/`.

## Scope of figure reproducibility

Analysis-derived figures are covered by the scripts in `code/`. The current
seasonal prevalence panel and water-quality time series have dedicated updated
scripts. Figure 1a's site map cannot be regenerated exactly until its map
source/coordinates are added. Final manual label positioning and transparent
background adjustments are presentation-only changes and do not alter data or
statistics.

## Manuscript-facing numerical audit

The archived tables reproduce the core manuscript statistics for Shannon
diversity, Bray–Curtis PCoA, restricted PERMANOVA/PERMDISP, environmental
vectors, turnover, blockwise R², random-forest reconstruction and PICRUSt2
functional profiles to the reported precision. The automated audit is
`code/00_validate_manuscript_invariants.py`; its machine-readable output is
`validation/manuscript_invariant_audit.json`.

Four stale values/references remain in the manuscript text: Fig. 3b should be
Fig. 4b; 2.23% should be 2.29% in two places; 21 of 30 should be 23 of 30; and
30.5 should be 30.3. Two structured-cross-validation citations point to [11]
instead of [27], and the unnumbered Shade et al. (2012) citation needs a matching
reference. These are text-level corrections; the underlying repository tables
and calculations are internally consistent.

## Remaining external dependencies

Raw FASTQ files are external to this archive and will be retrieved from SRA.
The complete sample–SAMN–SRR crosswalk is included for all 190 formal
libraries. The identifiers are unique and match the formal library manifest
exactly. Public accessibility without login was confirmed on 2026-10-02;
BioProject `PRJNA1513996` lists 190 SRA experiments.

## Blank-screened functional sensitivity

The primary combined EC prediction is archived and the full conservative
blank-screened path was rerun with PICRUSt2 2.6.3 and GLPK 5.0. The rerun used
23,787 ASVs after removing 109 blank-enriched ASVs and produced 555 pathways.
Across 14,028 sample pairs, pathway Bray–Curtis distances correlated at
Spearman rho = 0.9998773042 (`0.9999` at manuscript precision). All 22 compared
tests retained the same significance classification at 0.05. The compact
rerun outputs, downstream tables, direct comparison files, and one-command
workflow are included in this repository.
