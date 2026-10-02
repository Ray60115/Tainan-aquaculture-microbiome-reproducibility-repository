# Tainan 16S and PICRUSt2 reanalysis decisions

## Formal biological cohort

- The analytical cohort contains 168 unique sample identifiers from ponds T1–T4.
- Columns ending in `-1` and `-2` identify separate within-month sampling
  occasions, not technical replicates. The supplied profile matrix contains
  their canonical pond-by-month-by-layer arithmetic means, so no source
  occasion is treated as an independent profile-level observation.
- Columns ending in `-G` are excluded.
- `Water` is a blank water sample and is excluded from all biological analyses.

## Feature filtering

- ASVs with zero total abundance across the 168 formal samples are removed.
- Mitochondrial, chloroplast, and explicitly eukaryotic assignments are removed
  before Shannon diversity, Bray–Curtis dissimilarity, ordination, taxon
  summaries, random-forest reconstruction, and PICRUSt2.
- The primary input contains 23,896 prokaryotic ASVs.
- The single blank cannot support a formal prevalence-based contaminant model.
  ASVs whose blank relative abundance exceeds their maximum relative abundance
  in every formal sample are therefore removed only in a conservative
  sensitivity analysis.

## PICRUSt2

- The legacy pathway TSV is not patched or used.
- Functional profiles are regenerated from the filtered 168-sample ASV matrix
  with PICRUSt2 2.6.3.
- Placement is performed independently against the installed official
  bacterial and archaeal PICRUSt2-SC references. For the larger bacterial
  reference, SEPP 4.5.5 performs hierarchical decomposition, HMM search, and
  fragment alignment; EPA-ng 0.3.8 performs likelihood placement within each
  reference subtree; and SEPP's standard merger combines the jplace files
  before grafting with gappa 0.8.5. This compatibility route is used because
  the bundled pplacer executable segfaults before reading inputs on the
  current Linux 6.12 runtime and whole-bacterial-reference EPA-ng exceeds the
  14-GiB memory limit.
- For the smaller archaeal reference, the PICRUSt2 whole-reference EPA-ng
  route is used after HMM alignment, with an accumulated likelihood-weight
  threshold of 0.99, at most 100 placements per query, and query chunks of
  250 sequences, followed by grafting with gappa 0.8.5.
- The minimum aligned fraction is 0.80 in each independent domain-specific
  step. Of 23,896 input ASVs, one failed this criterion against the bacterial
  reference and six failed it against the archaeal reference. These
  domain-specific failures do not by themselves remove an ASV from the other
  reference.
- When an ASV has predictions from both references, PICRUSt2 retains the
  domain-specific prediction with the lower NSTI. Concordance with the
  supplied taxonomic domain is audited both before and after the prespecified
  maximum-NSTI filter.
- Hidden-state prediction of 16S marker-gene and EC-family copy numbers uses
  maximum parsimony (edge exponent 0.5; PICRUSt2 seed 100). The workflow
  normalizes ASV abundances by predicted 16S rRNA marker copy number, applies
  the default maximum NSTI of 2, and reconstructs unstratified MetaCyc v24
  pathways with MinPath and pathway gap filling.
- Weighted NSTI and the fraction of prokaryotic input reads retained after
  placement and NSTI filtering are reported.

## Analytical unit and dependence

- Pond identity and the nursery-versus-grow-out contrast are descriptive
  because the study contains only one nursery pond.
- Available surface and bottom profiles contribute equally to each pond-month
  mean.
- Season is tested within T2–T4 after controlling for pond, using 9,999
  circular shifts within pond.
- Sampling layer is tested in 75 complete surface–bottom pairs after
  controlling for pond-month, with labels swapped only within pairs.
- Meteorological indicators are tested after controlling for pond and season,
  using circular shifts of the shared sampling-month exposure sequence.
- Event p-values are Benjamini–Hochberg adjusted within analytical cohort.

## Functional modules and interpretation

- Six pathway modules were specified before inspection of functional-test
  results: nitrate reduction; sulfur oxidation and assimilation;
  methanogenesis; C1 assimilation and oxidation; fermentation and short-chain
  products; and compatible-solute metabolism.
- Module definitions are fixed in `input/pathway_module_definitions.csv`.
- PICRUSt2 outputs are predictions of genomic potential from 16S profiles.
  They are not measurements of absolute gene abundance, transcription,
  translation, metabolites, process rates, or pathway activity.
- Taxonomic–predicted-functional turnover concordance is not independent
  functional validation because both representations derive from the same 16S
  dataset.

## Unresolved laboratory metadata

Before submission, the authors must supply the DNA extraction method, filtered
water volume and filter pore size, primer sequences, PCR conditions, PacBio
chemistry and instrument, taxonomic classifier, and reference-database name
and version. These details are not inferable from the supplied ASV table or
analysis scripts and are not fabricated in the revised manuscript.
