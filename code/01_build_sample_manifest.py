#!/usr/bin/env python3
"""Build the 204-input, 190-library, and 168-profile sample manifests.

The formal analysis starts from the 190 biological U/D libraries listed in
Table_D_Sample_Metadata. Other workflow inputs are retained in the all-input
manifest for transparency but are not treated as unexplained missing samples.
"""

from __future__ import annotations

import argparse
from pathlib import Path, PurePosixPath

import numpy as np
import pandas as pd


EXPECTED_ALL_INPUTS = 204
EXPECTED_ANALYSIS_LIBRARIES = 190
EXPECTED_PROFILES = 168


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def classify_nonanalysis(sample_id: str) -> tuple[str, str]:
    if sample_id.endswith("-G"):
        return (
            "G_designated_nonanalysis_sample",
            "G-designated sample; not part of the formal U/D water-column cohort",
        )
    if sample_id == "Water":
        return "water_blank_control", "Water blank/control"
    if sample_id.startswith("Probiotic"):
        return "probiotic_reference", "Probiotic reference sample"
    if sample_id.startswith("SeaWater"):
        return "seawater_reference", "Seawater reference sample"
    return "unclassified_nonanalysis_input", "Not listed in formal sample metadata"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workflow-manifest", required=True, type=Path)
    parser.add_argument("--sample-metadata", required=True, type=Path)
    parser.add_argument("--dada2-qc", required=True, type=Path)
    parser.add_argument("--read-qc", required=True, type=Path)
    parser.add_argument("--trim-qc", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--include-source-paths",
        action="store_true",
        help="Retain machine-specific FASTQ paths in exported manifests.",
    )
    args = parser.parse_args()

    outdir = args.output_dir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    workflow = read_tsv(args.workflow_manifest)
    metadata = pd.read_csv(args.sample_metadata)
    dada2 = read_tsv(args.dada2_qc)
    read_qc = read_tsv(args.read_qc)
    trim_qc = read_tsv(args.trim_qc)

    if len(workflow) != EXPECTED_ALL_INPUTS:
        raise ValueError(f"Expected {EXPECTED_ALL_INPUTS} workflow inputs, found {len(workflow)}")
    if len(metadata) != EXPECTED_ANALYSIS_LIBRARIES:
        raise ValueError(
            f"Expected {EXPECTED_ANALYSIS_LIBRARIES} biological libraries, found {len(metadata)}"
        )
    if metadata["SampleID_raw"].duplicated().any():
        raise ValueError("SampleID_raw must be unique in the source-library metadata")
    if metadata["SampleID"].nunique() != EXPECTED_PROFILES:
        raise ValueError(
            f"Expected {EXPECTED_PROFILES} canonical profiles, found {metadata['SampleID'].nunique()}"
        )

    # Remove the QIIME 2 semantic-type row before converting numeric columns.
    dada2 = dada2.loc[~dada2["sample-id"].astype(str).str.startswith("#")].copy()
    dada2_numeric = [c for c in dada2.columns if c != "sample-id"]
    for column in dada2_numeric:
        dada2[column] = pd.to_numeric(dada2[column], errors="raise")

    if set(workflow["sample-id"]) != set(dada2["sample-id"]):
        raise ValueError("Workflow manifest and DADA2 QC sample identifiers do not match")
    if set(workflow["sample-id"]) != set(read_qc["sample"]):
        raise ValueError("Workflow manifest and read-QC sample identifiers do not match")
    if set(workflow["sample-id"]) != set(trim_qc["sample"]):
        raise ValueError("Workflow manifest and trim-QC sample identifiers do not match")

    path_column = (
        "absolute-filepath"
        if "absolute-filepath" in workflow.columns
        else "input-filename"
    )
    if path_column not in workflow.columns:
        raise ValueError(
            "Workflow manifest needs either absolute-filepath or input-filename"
        )
    manifest = workflow.rename(
        columns={"sample-id": "workflow_sample_id", path_column: "input_path"}
    ).copy()
    manifest["input_file"] = manifest["input_path"].map(
        lambda value: PurePosixPath(str(value)).name
    )

    metadata_for_join = metadata.rename(columns={"SampleID_raw": "workflow_sample_id"})
    manifest = manifest.merge(metadata_for_join, on="workflow_sample_id", how="left", validate="one_to_one")

    included = manifest["SampleID"].notna()
    manifest["workflow_class"] = "study_water_column_library"
    manifest["included_in_formal_190_library_cohort"] = included
    manifest["exclusion_reason_from_190"] = ""
    for index in manifest.index[~included]:
        group, reason = classify_nonanalysis(str(manifest.at[index, "workflow_sample_id"]))
        manifest.at[index, "workflow_class"] = group
        manifest.at[index, "exclusion_reason_from_190"] = reason

    manifest["canonical_profile_id"] = manifest["SampleID"]
    profile_sizes = metadata.groupby("SampleID").size()
    manifest["source_library_count_for_profile"] = manifest["canonical_profile_id"].map(profile_sizes)
    manifest["included_in_168_profile_matrix"] = included
    manifest["profile_aggregation_rule"] = np.where(
        included,
        "Arithmetic mean of source-library ASV counts within canonical profile",
        "",
    )
    manifest["biosample_accession"] = ""
    manifest["sra_run_accession"] = ""
    manifest["accession_status_at_audit"] = np.where(
        included,
        "Run 01_attach_sra_accessions.py to attach and validate the archived NCBI accessions",
        "Not part of the formal 190-library analytical cohort",
    )

    dada2 = dada2.rename(
        columns={
            "sample-id": "workflow_sample_id",
            "input": "dada2_input_reads",
            "filtered": "dada2_filtered_reads",
            "percentage of input passed filter": "dada2_percent_passed_filter",
            "denoised": "dada2_denoised_reads",
            "non-chimeric": "dada2_nonchimeric_reads",
            "percentage of input non-chimeric": "dada2_percent_input_nonchimeric",
        }
    )
    read_qc = read_qc[["sample", "num_seqs", "sum_len", "avg_len", "Q20(%)", "Q30(%)"]].rename(
        columns={
            "sample": "workflow_sample_id",
            "num_seqs": "raw_read_count",
            "sum_len": "raw_total_bases",
            "avg_len": "raw_mean_read_length",
            "Q20(%)": "raw_Q20_percent",
            "Q30(%)": "raw_Q30_percent",
        }
    )
    trim_qc = trim_qc.rename(
        columns={
            "sample": "workflow_sample_id",
            "input_reads": "trim_input_reads",
            "demuxed_reads": "primer_matched_reads",
        }
    )
    manifest = manifest.merge(dada2, on="workflow_sample_id", how="left", validate="one_to_one")
    manifest = manifest.merge(read_qc, on="workflow_sample_id", how="left", validate="one_to_one")
    manifest = manifest.merge(trim_qc, on="workflow_sample_id", how="left", validate="one_to_one")

    if manifest["dada2_nonchimeric_reads"].isna().any():
        raise ValueError("DADA2 QC values are missing after the sample join")
    if (manifest["dada2_nonchimeric_reads"] <= 0).any():
        raise ValueError("At least one workflow input has zero non-chimeric reads")

    preferred = [
        "workflow_sample_id",
        "input_file",
        "workflow_class",
        "included_in_formal_190_library_cohort",
        "exclusion_reason_from_190",
        "canonical_profile_id",
        "Pond",
        "Date",
        "Occasion",
        "Layer",
        "HasOccasionLabel",
        "source_library_count_for_profile",
        "included_in_168_profile_matrix",
        "profile_aggregation_rule",
        "biosample_accession",
        "sra_run_accession",
        "accession_status_at_audit",
        "raw_read_count",
        "raw_total_bases",
        "raw_mean_read_length",
        "raw_Q20_percent",
        "raw_Q30_percent",
        "trim_input_reads",
        "primer_matched_reads",
        "dada2_input_reads",
        "dada2_filtered_reads",
        "dada2_percent_passed_filter",
        "dada2_denoised_reads",
        "dada2_nonchimeric_reads",
        "dada2_percent_input_nonchimeric",
        "input_path",
    ]
    manifest = manifest[preferred]

    analysis_190 = manifest.loc[manifest["included_in_formal_190_library_cohort"]].copy()
    if len(analysis_190) != EXPECTED_ANALYSIS_LIBRARIES:
        raise ValueError("The included study-library count is not 190")

    profile_rows: list[dict[str, object]] = []
    for profile_id, group in analysis_190.groupby("canonical_profile_id", sort=False):
        source_ids = group["workflow_sample_id"].astype(str).tolist()
        profile_rows.append(
            {
                "canonical_profile_id": profile_id,
                "Pond": group["Pond"].iloc[0],
                "Date": int(group["Date"].iloc[0]),
                "Layer": group["Layer"].iloc[0],
                "source_library_count": len(group),
                "source_library_ids": ";".join(source_ids),
                "within_month_occasion_averaging_applied": len(group) > 1,
                "profile_aggregation_rule": "Arithmetic mean of source-library ASV counts",
                "total_source_nonchimeric_reads": int(group["dada2_nonchimeric_reads"].sum()),
                "mean_source_nonchimeric_reads": float(group["dada2_nonchimeric_reads"].mean()),
                "minimum_source_nonchimeric_reads": int(group["dada2_nonchimeric_reads"].min()),
                "maximum_source_nonchimeric_reads": int(group["dada2_nonchimeric_reads"].max()),
            }
        )
    profiles_168 = pd.DataFrame(profile_rows)
    if len(profiles_168) != EXPECTED_PROFILES:
        raise ValueError("The canonical-profile count is not 168")

    summary = pd.DataFrame(
        [
            ("workflow_inputs", len(manifest), "All inputs listed by the sequencing workflow"),
            ("formal_biological_libraries", len(analysis_190), "Official analysis starts here"),
            ("excluded_G_designated", int((manifest["workflow_class"] == "G_designated_nonanalysis_sample").sum()), "Not U/D water-column libraries"),
            ("excluded_water_blank", int((manifest["workflow_class"] == "water_blank_control").sum()), "Water blank/control"),
            ("excluded_probiotic_reference", int((manifest["workflow_class"] == "probiotic_reference").sum()), "Reference samples"),
            ("excluded_seawater_reference", int((manifest["workflow_class"] == "seawater_reference").sum()), "Reference samples"),
            ("canonical_profiles", len(profiles_168), "Pond-month-layer profiles"),
            ("profiles_with_one_source_library", int((profiles_168["source_library_count"] == 1).sum()), "No within-month occasion averaging"),
            ("profiles_with_two_source_libraries", int((profiles_168["source_library_count"] == 2).sum()), "Two within-month sampling occasions averaged"),
            ("zero_nonchimeric_workflow_inputs", int((manifest["dada2_nonchimeric_reads"] == 0).sum()), "Expected 0"),
        ],
        columns=["metric", "value", "interpretation"],
    )

    exported_manifest = manifest.copy()
    exported_analysis = analysis_190.copy()
    if not args.include_source_paths:
        exported_manifest = exported_manifest.drop(columns=["input_path"])
        exported_analysis = exported_analysis.drop(columns=["input_path"])
    exported_manifest.to_csv(outdir / "sample_manifest_all_204.tsv", sep="\t", index=False)
    exported_analysis.to_csv(outdir / "sample_manifest_analysis_libraries_190.tsv", sep="\t", index=False)
    profiles_168.to_csv(outdir / "sample_manifest_profiles_168.tsv", sep="\t", index=False)
    summary.to_csv(outdir / "sample_manifest_summary.tsv", sep="\t", index=False)

    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
