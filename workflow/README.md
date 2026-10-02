# Sequencing workflow provenance

`HiFi-16S-workflow_v0.9_commit_d29598f/` is a source snapshot of the official
Pacific Biosciences HiFi-16S workflow at commit
`d29598f296b222c85283c49f50a0b34ef0ee4a76`. Upstream demonstration outputs and
tutorial images were omitted; all executable modules, scripts, environment
files, documentation, and the upstream license are retained.

The recorded study run used workflow version 0.9. Its `main.nf` is byte-for-byte
identical to `main.nf.used_snapshot`. `nextflow.config.used_snapshot` is the
study-specific configuration and differs from the upstream default because it
records the database paths, resources, and enabled options used for this run.

`sample_input_204.tsv` contains sample identifiers and filenames only. Local
absolute paths were intentionally removed. Populate a two-column manifest with
`sample-id` and `absolute-filepath`, or generate one from the final SRA
crosswalk, before re-running the workflow.

Machine-specific paths in `parameters_used.txt` were replaced with the
`<HIFI16S_WORKFLOW_DIR>` placeholder; all analytical parameter values are
unchanged.

The study configuration points to Greengenes2 2024.09, GTDB R220, and SILVA
138.2, matching the manuscript. The unmodified upstream
`scripts/dada2_assign_tax.R` contains legacy literal labels `GTDB r207` and
`Silva 138.1` in its optional assignment-source column even though the actual
database files are supplied by the R220/138.2 configuration paths. These are
upstream display labels, not the databases used for classification. The used
workflow snapshot is intentionally preserved byte-for-byte rather than editing
that provenance artifact.
