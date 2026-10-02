# Manuscript–repository audit (2026-10-01)

## Overall status

**PASS.** The formal sample cohorts, ASV filtering, principal
taxonomic statistics, reconstruction results, and primary PICRUSt2 results are
internally consistent and match the submitted manuscript to the reported
precision. Public accessibility of all 190 assigned SAMN/SRR records was
confirmed without login on 2026-10-02. The primary combined EC predictions have
now been added and the blank-screened functional sensitivity analysis has been
rerun independently. Its pathway-distance Spearman correlation is 0.999877
(`0.9999` at manuscript precision), and none of the 22 compared significance
classifications changed.

## Verified data flow

| Stage | Verified count | Audit result |
|---|---:|---|
| Workflow manifest and QC | 204 inputs | Pass |
| Workflow feature-table columns | 199 | Pass; the five absent workflow inputs are Probiotic1–3 and SeaWater1–2 |
| Formal biological libraries | 190 | Pass; all are present in the feature table |
| Formal libraries with unique public SAMN / SRR | 190 / 190 | Pass; all records publicly accessible on 2026-10-02 |
| Surface / bottom libraries | 102 / 88 | Pass |
| Canonical pond-month-layer profiles | 168 | Pass |
| Profiles with one / two source occasions | 146 / 22 | Pass |
| Pond-months / complete layer pairs | 93 / 75 | Pass |
| Season-valid all-pond samples / groups | 160 / 89 | Pass |
| Season-valid T2–T4 samples / groups | 122 / 70 | Pass |
| Starting ASVs | 29,461 | Pass |
| ASVs absent from all formal profiles | 4,242 | Pass |
| Positive mitochondrial / chloroplast / explicit eukaryotic ASVs removed | 728 / 590 / 5 | Pass |
| Retained prokaryotic ASVs | 23,896 | Pass |
| Resolved genus-level lineages | 2,673 | Pass |
| Predicted MetaCyc pathways | 556 | Pass |
| Blank-screened ASVs / MetaCyc pathways | 23,787 / 555 | Pass |
| Blank-screened pathway-distance rho / classification changes | 0.999877 / 0 | Pass |

Arithmetic averaging of the 190 source libraries reconstructs all 168 formal
profiles exactly: the maximum absolute ASV-cell difference is zero. The strict
analysis table now contains only the 168 formal profiles. The eight `-G`
columns and water blank are stored separately so that they cannot be mistaken
for formal analysis samples.

## Verified manuscript-facing results

The archived outputs reproduce the reported values for:

- Shannon diversity layer and season tests;
- Bray–Curtis PCoA axis percentages;
- restricted PERMANOVA and PERMDISP analyses;
- environmental-vector fits for pH and water temperature;
- month-to-month turnover;
- blockwise unique and full-model R²;
- random-forest cross-validation and grouped importance;
- six-algorithm submitted reference tables; and
- primary PICRUSt2 placement, pathway, PCoA, functional turnover, and module
  analyses; and
- the conservative blank-screened functional sensitivity claim (`rho = 0.9999`
  at manuscript precision, with no significance-classification changes).

Run `python code/00_validate_manuscript_invariants.py` for the fast automated
audit. Its machine-readable result is
`validation/manuscript_invariant_audit.json`.

## Repository corrections made during this audit

- Split the former mixed 168-profile table into a strict 168-profile table and
  a separate nine-column non-formal/control table.
- Preserved the 199-column workflow feature table under an explicit filename.
- Replaced “technical replicate” language with “separate within-month sampling
  occasions,” matching the manuscript.
- Updated current figure and supplementary-table numbering and output names.
- Replaced Word-dependent taxon-table construction with a public,
  machine-readable lineage-definition file.
- Added generation of Additional file 1 and verified all three selected-taxon
  tables exactly against the archived reference results.
- Expanded dependency and input preflight checks so missing software or
  PICRUSt2 reference assets cause an explicit stop rather than a partial run.
- Updated plotted environmental labels to distinguish water temperature from
  maximum air temperature.
- Documented the workflow's legacy GTDB/SILVA display labels while preserving
  the used upstream snapshot unchanged. The configured databases are GTDB R220
  and SILVA 138.2, as stated in the manuscript.
- Added the primary combined EC predictions, reran the complete conservative
  blank-screened PICRUSt2 pathway analysis, archived compact rerun/reference
  outputs, and added a one-command sensitivity workflow.

## Manuscript text requiring correction

These are manuscript-text mismatches, not disagreements among repository data
files:

1. Methods: `Fig. 3b` should be `Fig. 4b`.
2. Results and Discussion: the combined grow-out mean abundance of the selected
   30 lineages is 2.2854%, so `2.23%` should be `2.29%` in both places.
3. Results: 23 of the 30 selected lineages, not 21, have expected counts below
   20 at the stated median depth. The stated median expected count of 9.7 is
   reproduced from sample-level grow-out means.
4. Results: the median seasonal prevalence range in the 0.01–<0.03% abundance
   bin is 30.314 points, so `30.5` should be `30.3`.
5. The two statements about structured cross-validation cite `[11]`; the
   matching Roberts et al. paper is reference `[27]`.
6. `(Shade et al., 2012)` has no matching 2012 reference-list entry and needs
   verification or replacement.
7. Data Availability must replace “upon reasonable request” with the public
   repository/Zenodo DOI and final SRA accessions.

## Transparent reproducibility qualifications

- The recovered archive did not contain the original rarefaction script or
  random-number stream. The supplied deterministic implementation reproduces
  cohort sizes, direction, and significance conclusions, but some Monte Carlo
  values differ in the last digits. Submitted and rerun tables are retained
  side by side.
- The reconstructed six-model script matches 46 of 50 submitted rows at four
  decimals. Four GAM rows containing water-quality or disturbance smooths have
  small environment-dependent differences; model ordering and interpretation
  are unchanged. Submitted and rerun tables are retained side by side.
- Figure 1a's map source asset is not archived, although the climate panels and
  all analysis-derived figures have code paths.

The complete 190-library accession crosswalk is included and public access is
verified. The blank-screened functional sensitivity is independently
reproducible. The repository is suitable for public GitHub release and Zenodo
archiving. The depositor must select the applicable license metadata in Zenodo,
and applicable manuscript text corrections above should be handled in the
manuscript.
