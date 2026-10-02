# Release notes

This package contains the scientific analysis inputs, executable scripts,
reference outputs, and validation summaries used to reproduce and check the
reported analyses.

Machine-specific absolute paths in provenance records were replaced with
explicit placeholders. Analytical parameters, hashes of scientific inputs,
study data, model settings, and reference results were not changed.

Public accessibility of all 190 SRA experiments was confirmed without login on
2026-10-02. The package is ready for GitHub release and Zenodo archiving. At
deposit, the uploader must select the applicable license metadata in Zenodo.

This audited repository also separates the strict 168-profile analysis matrix
from nine non-formal `-G`/water columns, records why the workflow has 204 inputs
but the exported feature table has 199 columns, updates all current manuscript
figure/table numbering, and adds an automated manuscript-invariant audit.

The package includes the primary `combined_EC_predicted.tsv.gz`, the
conservative blank-screened PICRUSt2 sensitivity analysis, compact screened
PICRUSt2 and downstream comparison outputs, and the one-command
`run_blank_screened_functional_sensitivity.sh` workflow. The rerun
reproduces the manuscript's pathway-distance correlation (`rho = 0.9999` at
reported precision) and shows zero changes among 22 significance
classifications.
