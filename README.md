# Tainan aquaculture microbiome reproducibility repository

This repository accompanies the manuscript **“Monthly microbiome surveillance
in subtropical aquaculture ponds resolves seasonal community structure but not
weather-driven disturbances.”** It consolidates the sequencing-workflow
provenance, sample mapping, ASV filtering, taxonomic statistics, grouped
cross-validation, PICRUSt2 analysis, and figure-generation code used for the
submitted manuscript.

The repository contains the audited inputs, executable analysis code,
reference results, and validation summaries needed to reproduce and check the
reported analyses. Machine-specific absolute paths in provenance records are
replaced with explicit placeholders without changing analytical settings.

## What is included

- the complete PacBio `HiFi-16S-workflow` v0.9 source snapshot at commit
  `d29598f296b222c85283c49f50a0b34ef0ee4a76`, including modules, scripts,
  environments, and its BSD-3-Clause-Clear license;
- the study-specific Nextflow configuration and recorded parameter log;
- QC summaries for all 204 workflow inputs;
- explicit mappings from 204 workflow inputs to 190 formal biological
  libraries and 168 canonical pond-month-layer profiles;
- the 199-column workflow feature table, a strict 168-profile formal-analysis
  table, the nine separately archived non-formal `-G`/water columns, finalized
  metadata, water-quality data, and climate data;
- scripts for ASV filtering, rarefaction, Shannon diversity, Bray–Curtis PCoA,
  dependence-aware PERMANOVA/PERMDISP, environmental-vector fitting,
  blockwise variation partitioning, turnover, taxonomic summaries, grouped
  cross-validation, model sensitivity, PICRUSt2, and plotting; and
- reference result tables for numerical comparison.

The upstream workflow is also available from the official PacBio repository:
<https://github.com/PacificBiosciences/HiFi-16S-workflow>.

## Verified data flow

1. The sequencing workflow lists **204 inputs**.
2. Fourteen inputs are not part of the formal water-column cohort: eight
   `-G` records, one water blank, three probiotic references, and two seawater
   references.
   The five probiotic/seawater references are listed in the workflow and QC
   records but are absent from the exported feature table; all 190 formal
   libraries are present in that table.
3. The formal analysis therefore begins with **190 biological libraries**.
4. Twenty-two pond-month-layer profiles have two source libraries; arithmetic
   averaging reconstructs **168 canonical profiles** exactly (maximum cell
   difference = 0).
5. The workflow ASV table contains **29,461 ASVs**. Removing 4,242 variants
   absent from all 168 profiles leaves 25,219 positive ASVs. Removing 1,323
   mitochondrial, chloroplast, or explicitly eukaryotic ASVs leaves
   **23,896 prokaryotic ASVs**.

## Directory guide

- `workflow/`: full upstream workflow snapshot, study configuration, parameters,
  sample manifest, and compact QC tables.
- `data/primary/`: compressed feature tables and finalized metadata.
- `data/manifests/`: the audited 204 → 190 → 168 mappings and 29,461 → 23,896
  ASV audit.
- `data/accessions/`: consolidated NCBI SRA submission metadata for all 190
  formal libraries.
- `data/supporting/`: library-depth and submitted rarefaction summaries.
- `code/`: analysis and plotting scripts.
- `input/`: prespecified PICRUSt2 pathway-module definitions.
- `picrust2_primary/`: archived primary PICRUSt2 marker and EC predictions plus
  compact outputs needed for downstream functional analyses without repeating
  phylogenetic placement.
- `picrust2_blank_screened/`: compact rerun outputs after conservatively
  removing 109 blank-enriched ASVs.
- `reference_results/`: submitted/reference taxonomic, functional, and model
  sensitivity results.
- `environment/`: PICRUSt2/SEPP environment and compatibility records.
- `validation/`: audit reports, deterministic rarefaction reruns, and compact
  submitted-versus-rerun comparisons.
- `MANUSCRIPT_TO_CODE_MAP.md`: method, figure, and table crosswalk.
- `REPRODUCIBILITY_NOTES.md`: documented limitations, provenance notes, and
  release metadata guidance.

## Installation

The manuscript used Python 3.12.13. For the downstream analyses:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

PICRUSt2 uses the separately recorded environment under `environment/`.
Nextflow/QIIME 2/DADA2 dependencies are defined by the included upstream
workflow.

## Reproduce the core analysis

Run from the repository root:

```bash
bash run_core_analysis.sh
```

This command rebuilds the sample/ASV audit, repeated-rarefaction sensitivity,
dependence-aware taxonomic results, random-forest reconstruction, six-model
sensitivity analysis, and analysis-generated figures. Restricted tests use
9,999 permutations and may take time.

Before a full rerun, the manuscript-facing counts and archived numerical
results can be checked quickly with:

```bash
python code/00_validate_manuscript_invariants.py
```

The current automated audit status is `PASS`; see
`MANUSCRIPT_REPOSITORY_AUDIT_20261001.md` for the verified scope. The SRA
public-access check and the blank-screened functional-sensitivity rerun are
complete.

## Reproduce the functional analysis

Rerun downstream functional statistics from the archived primary PICRUSt2
outputs. If the core analysis has already been run, its reconstructed turnover
table is used; otherwise, the archived reference turnover table is used:

```bash
bash run_functional_analysis.sh
```

To rerun the complete conservative blank-screened sensitivity path, including
PICRUSt2 metagenome/pathway inference, downstream functional statistics, and
the primary-versus-screened comparison, activate the recorded PICRUSt2
environment and run:

```bash
PICRUST2_PROCESSES=4 bash run_blank_screened_functional_sensitivity.sh
```

The archived rerun retained 23,787 ASVs and produced 555 MetaCyc pathways. Its
pairwise pathway Bray–Curtis distances were nearly identical to the primary
analysis (Spearman rho = 0.999877; `0.9999` at manuscript precision), and none
of the 22 compared test classifications changed at the 0.05 threshold.

To repeat phylogenetic placement and prediction from the ASV sequences, use
`code/03_run_picrust2_primary.sh` with the recorded PICRUSt2 environment and
compatibility helpers. Set `BAC_EC_TABLE` and `ARC_EC_TABLE` to the installed
official PICRUSt2-SC EC reference tables. These large official reference assets
are not duplicated here.

## Re-run the upstream sequencing workflow

Raw FASTQ files are not duplicated in this archive because they are deposited
in SRA. After populating local FASTQ paths (or generating them from the
sample–SRR crosswalk), use the included full workflow snapshot and the
study-specific configuration:

```bash
cd workflow/HiFi-16S-workflow_v0.9_commit_d29598f
cp ../nextflow.config.used_snapshot nextflow.config
nextflow run main.nf --input /path/to/sample_input_with_local_paths.tsv \
  --metadata /path/to/workflow_metadata.tsv --outdir results -profile conda
```

`workflow/sample_input_204.tsv` intentionally contains filenames rather than
machine-specific absolute paths.

## Interpretation limits

- Pond identity and the nursery-versus-grow-out contrast are descriptive.
- Cross-validation evaluates reconstruction within the sampled farm and study
  period, not chronological forecasting or transfer to new ponds.
- PICRUSt2 outputs are predictions of genomic potential derived from 16S data,
  not direct measurements of genes, expression, metabolites, or activity.
- See `REPRODUCIBILITY_NOTES.md` for the completed blank-screened sensitivity
  verification, map-source item, licensing guidance, and numerical provenance
  notes retained transparently with the submitted reference tables.

## Data accession

BioProject: `PRJNA1513996`; SRA study: `SRP727697`. The complete 190-library
sample–BioSample–SRA-run crosswalk is included under `data/manifests/`.
Public accessibility without login was confirmed on 2026-10-02: the BioProject
lists 190 SRA experiments. See
<https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1513996> and
<https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP727697>.
