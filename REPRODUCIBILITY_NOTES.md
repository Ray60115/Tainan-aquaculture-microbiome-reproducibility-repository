# Reproducibility and release notes

## Completed accession check

BioProject `PRJNA1513996`, SRA study `SRP727697`, and the complete set of 190
SRA experiments were confirmed accessible without login on 2026-10-02. The
verification record is `validation/sra_public_access_verified.txt`.

The complete 190-library sample–SAMN–SRR crosswalk has already been added. It
contains one unique BioSample and SRA run accession for every formal library
and matches the formal library manifest exactly.

The manuscript Data Availability statement should cite the public SRA
accessions and the repository DOI after Zenodo issues it.

## Completed blank-screened functional sensitivity verification

The primary `combined_EC_predicted.tsv.gz` has been added and the conservative
blank-screened analysis was rerun with PICRUSt2 2.6.3 and GLPK 5.0. Removing
109 blank-enriched ASVs left 23,787 input ASVs and produced 555 MetaCyc
pathways, compared with 556 pathways in the primary analysis. Across 14,028
sample pairs, primary and screened pathway Bray–Curtis distances had Spearman
rho = 0.9998773042 (`0.9999` at manuscript precision). All 22 compared tests
retained the same significance classification at 0.05.

The compact PICRUSt2 rerun outputs are under `picrust2_blank_screened/`, the
downstream screened results are under
`reference_results/functional_blank_screened/`, and the direct comparison is
under `reference_results/functional_sensitivity/`. The complete rerun command
is `run_blank_screened_functional_sensitivity.sh`.

## Other scope note

Figure 1a (study-site map) uses a separately prepared geographic asset; the
farm coordinates/source map are not present in the analysis inputs. All
analysis-derived panels remain covered by the archived code.

The compact PICRUSt2 primary and blank-screened outputs are intentionally
retained so downstream functional analyses can be checked without repeating
phylogenetic placement.

## Zenodo metadata

Zenodo requires the depositor to select a license. The bundled PacBio workflow
retains its own BSD-3-Clause-Clear license; the depositor should select the
applicable license metadata for the study-authored code, data, and
documentation. Zenodo will issue a DOI upon publication unless one is reserved
in advance.

## Manuscript text corrections identified by the repository audit

1. Methods: `Fig. 3b` for the first 20 seasonal lineages should be `Fig. 4b`.
2. Results and Discussion: the selected 30 lineages total 2.2854%, which rounds
   to `2.29%`, not `2.23%`.
3. Results: at the stated median depth, 23 of 30 selected lineages—not 21—have
   expected counts below 20. The reported median expected count of 9.7 is
   otherwise reproduced using sample-level grow-out means.
4. Results: the 0.01–<0.03% abundance bin has a median seasonal prevalence
   range of 30.314 points, which rounds to `30.3`, not `30.5`.
5. The structured-cross-validation statements cite reference `[11]`, whereas
   reference `[27]` is the Roberts et al. cross-validation paper.
6. `(Shade et al., 2012)` does not have a matching 2012 entry in the numbered
   reference list; verify or replace that citation.
7. Replace “upon reasonable request” in Data Availability with the public
   GitHub/Zenodo DOI and final SRA accessions before resubmission.

## Numerical provenance notes

- The submitted rarefaction summaries are preserved in
  `data/supporting/Additional_file_3_rarefaction_summary.csv` and the manuscript
  table. The recovered archive did not contain the original rarefaction script
  or its exact random-number stream. `code/03_run_rarefaction_sensitivity.py`
  now provides a deterministic implementation of the documented procedure;
  cohort sizes, direction, and significance conclusions reproduce, while
  Monte Carlo summary values may differ slightly in the last reported digits.
- The submitted six-model values are preserved under
  `reference_results/sensitivity/*_submitted.csv`. The reconstructed model
  script reproduces ridge, random forest, gradient boosting, XGBoost, SVR, and
  the structure/season GAM rows to four decimals. Four GAM rows involving
  water-quality/disturbance smooths show small numerical differences under the
  reconstructed environment; the interpretation and model ordering are
  unchanged. Retaining both submitted and rerun tables makes this visible.
