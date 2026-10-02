#!/usr/bin/env python3
"""Validate and attach SRA/BioSample accessions to the audited manifests."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


REQUIRED = {
    "accession",
    "study",
    "bioproject_accession",
    "biosample_accession",
    "sample_name",
    "filename",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sra-metadata", nargs="+", required=True, type=Path)
    parser.add_argument("--library-metadata", required=True, type=Path)
    parser.add_argument("--manifest-dir", required=True, type=Path)
    parser.add_argument("--output-metadata", required=True, type=Path)
    args = parser.parse_args()

    pieces = [pd.read_csv(path, sep="\t", dtype=str) for path in args.sra_metadata]
    sra = pd.concat(pieces, ignore_index=True)
    missing_columns = REQUIRED - set(sra.columns)
    if missing_columns:
        raise ValueError(f"SRA metadata lacks columns: {sorted(missing_columns)}")

    libraries = pd.read_csv(args.library_metadata, dtype=str)
    expected = libraries["SampleID_raw"].tolist()
    if len(expected) != 190 or len(set(expected)) != 190:
        raise ValueError("Formal library metadata must contain 190 unique source libraries.")
    if len(sra) != 190 or sra["sample_name"].nunique() != 190:
        raise ValueError("Combined SRA metadata must contain 190 unique sample names.")
    if sra["biosample_accession"].nunique() != 190:
        raise ValueError("BioSample accessions are not one-to-one with formal libraries.")
    if sra["accession"].nunique() != 190:
        raise ValueError("SRA run accessions are not one-to-one with formal libraries.")
    if set(expected) != set(sra["sample_name"]):
        raise ValueError(
            "SRA/formal-library sample mismatch. Missing from SRA: "
            f"{sorted(set(expected) - set(sra['sample_name']))}; extra in SRA: "
            f"{sorted(set(sra['sample_name']) - set(expected))}"
        )
    if set(sra["bioproject_accession"]) != {"PRJNA1513996"}:
        raise ValueError("Unexpected BioProject accession.")
    if set(sra["study"]) != {"SRP727697"}:
        raise ValueError("Unexpected SRA study accession.")
    expected_filenames = sra["sample_name"] + ".gz"
    if not expected_filenames.equals(sra["filename"]):
        raise ValueError("At least one SRA filename does not match sample_name + '.gz'.")

    order = pd.Series(range(len(expected)), index=expected)
    sra = sra.assign(_order=sra["sample_name"].map(order)).sort_values("_order")
    sra = sra.drop(columns="_order").reset_index(drop=True)
    args.output_metadata.parent.mkdir(parents=True, exist_ok=True)
    sra.to_csv(args.output_metadata, sep="\t", index=False)

    crosswalk = sra[
        [
            "sample_name",
            "biosample_accession",
            "accession",
            "study",
            "bioproject_accession",
            "filename",
        ]
    ].rename(
        columns={
            "sample_name": "workflow_sample_id",
            "accession": "sra_run_accession",
            "study": "sra_study_accession",
            "filename": "sra_filename",
        }
    )
    crosswalk.to_csv(
        args.manifest_dir / "sample_SAMN_SRR_crosswalk.tsv", sep="\t", index=False
    )

    accession_columns = crosswalk[
        ["workflow_sample_id", "biosample_accession", "sra_run_accession"]
    ]
    for name in (
        "sample_manifest_all_204.tsv",
        "sample_manifest_analysis_libraries_190.tsv",
    ):
        path = args.manifest_dir / name
        manifest = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        manifest = manifest.drop(
            columns=["biosample_accession", "sra_run_accession"], errors="ignore"
        ).merge(accession_columns, on="workflow_sample_id", how="left", validate="one_to_one")
        manifest["biosample_accession"] = manifest["biosample_accession"].fillna("")
        manifest["sra_run_accession"] = manifest["sra_run_accession"].fillna("")
        formal = manifest["included_in_formal_190_library_cohort"].eq("True")
        manifest.loc[formal, "accession_status_at_audit"] = (
            "Accession assigned; confirm public accessibility without login before resubmission"
        )
        manifest.to_csv(path, sep="\t", index=False)

    profile_path = args.manifest_dir / "sample_manifest_profiles_168.tsv"
    profiles = pd.read_csv(profile_path, sep="\t", dtype=str, keep_default_na=False)
    by_sample = crosswalk.set_index("workflow_sample_id")
    profiles["source_biosample_accessions"] = profiles["source_library_ids"].map(
        lambda values: ";".join(
            by_sample.loc[item, "biosample_accession"] for item in values.split(";")
        )
    )
    profiles["source_sra_run_accessions"] = profiles["source_library_ids"].map(
        lambda values: ";".join(
            by_sample.loc[item, "sra_run_accession"] for item in values.split(";")
        )
    )
    profiles.to_csv(profile_path, sep="\t", index=False)

    print(
        "Attached 190 unique BioSample/SRA run accessions to the audited manifests "
        f"({len(pieces[0])}+{sum(len(piece) for piece in pieces[1:])})."
    )


if __name__ == "__main__":
    main()
