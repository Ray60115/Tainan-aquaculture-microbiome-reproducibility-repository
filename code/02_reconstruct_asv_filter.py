#!/usr/bin/env python3
"""Reconstruct and audit the 29,461 -> 25,219 -> 23,896 ASV flow."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


ANNOTATION_COLUMNS = ("id", "Sequence", "Taxon")
ORGANELLE_EUKARYOTE_PATTERN = re.compile(
    r"mitochond|chloroplast|d__eukary", flags=re.IGNORECASE
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-asv", required=True, type=Path)
    parser.add_argument("--averaged-asv", required=True, type=Path)
    parser.add_argument("--sample-metadata", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    outdir = args.output_dir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # The raw QIIME 2 CSV has a semantic-type row immediately below the header.
    raw = pd.read_csv(args.raw_asv, skiprows=[1], low_memory=False)
    averaged = pd.read_csv(args.averaged_asv, low_memory=False)
    metadata = pd.read_csv(args.sample_metadata)

    for label, table in (("raw", raw), ("averaged", averaged)):
        missing = [column for column in ANNOTATION_COLUMNS if column not in table.columns]
        if missing:
            raise ValueError(f"{label} ASV table lacks annotation columns: {missing}")
        if table["id"].duplicated().any():
            raise ValueError(f"{label} ASV table contains duplicated ASV identifiers")

    if len(raw) != 29461 or len(averaged) != 29461:
        raise ValueError(
            f"Expected 29,461 ASVs in both tables; found {len(raw):,} and {len(averaged):,}"
        )
    if set(raw["id"]) != set(averaged["id"]):
        raise ValueError("Raw and averaged ASV identifier sets are not identical")

    biological_libraries = metadata["SampleID_raw"].astype(str).tolist()
    canonical_profiles = metadata["SampleID"].drop_duplicates().astype(str).tolist()
    if len(biological_libraries) != 190 or len(canonical_profiles) != 168:
        raise ValueError("Expected 190 biological libraries and 168 canonical profiles")

    missing_raw = [sample for sample in biological_libraries if sample not in raw.columns]
    missing_averaged = [sample for sample in canonical_profiles if sample not in averaged.columns]
    if missing_raw or missing_averaged:
        raise ValueError(
            f"Missing raw libraries: {missing_raw}; missing averaged profiles: {missing_averaged}"
        )

    raw_by_id = raw.set_index("id")
    averaged_by_id = averaged.set_index("id")
    ordered_ids = averaged["id"].astype(str).tolist()
    raw_by_id = raw_by_id.loc[ordered_ids]
    averaged_by_id = averaged_by_id.loc[ordered_ids]

    aggregation_rows: list[dict[str, object]] = []
    global_max_difference = 0.0
    for profile_id, group in metadata.groupby("SampleID", sort=False):
        source_ids = group["SampleID_raw"].astype(str).tolist()
        calculated = raw_by_id[source_ids].astype(float).mean(axis=1)
        supplied = averaged_by_id[str(profile_id)].astype(float)
        difference = (calculated - supplied).abs()
        maximum = float(difference.max())
        global_max_difference = max(global_max_difference, maximum)
        aggregation_rows.append(
            {
                "canonical_profile_id": profile_id,
                "source_library_count": len(source_ids),
                "source_library_ids": ";".join(source_ids),
                "maximum_absolute_cell_difference": maximum,
                "number_of_nonmatching_ASV_cells": int((difference > 0).sum()),
                "exact_match": maximum == 0.0,
            }
        )

    aggregation_check = pd.DataFrame(aggregation_rows)
    if global_max_difference != 0.0:
        raise ValueError(
            f"Reconstructed 190-to-168 aggregation differs from supplied table; max difference {global_max_difference}"
        )

    counts = averaged[canonical_profiles].fillna(0).astype(float)
    if (counts < 0).any().any():
        raise ValueError("Negative ASV counts were detected")
    biological_total = counts.sum(axis=1)
    biological_prevalence = (counts > 0).mean(axis=1)
    positive = biological_total > 0

    taxonomy = averaged["Taxon"].fillna("").astype(str)
    sequences = averaged["Sequence"].fillna("").astype(str).str.upper()
    mitochondria = taxonomy.str.contains("mitochond", case=False, regex=False)
    chloroplast = taxonomy.str.contains("chloroplast", case=False, regex=False)
    explicit_eukaryote = taxonomy.str.contains("d__eukary", case=False, regex=False)
    organelle_or_eukaryote = taxonomy.str.contains(ORGANELLE_EUKARYOTE_PATTERN)
    invalid_sequence = ~sequences.str.fullmatch(r"[ACGTN]+")
    retained = positive & ~organelle_or_eukaryote & ~invalid_sequence

    if int(retained.sum()) != 23896:
        raise ValueError(f"Expected 23,896 retained ASVs, found {int(retained.sum()):,}")

    exclusion_reason: list[str] = []
    for index in averaged.index:
        reasons: list[str] = []
        if not bool(positive.iloc[index]):
            reasons.append("absent_from_all_168_biological_profiles")
        else:
            if bool(mitochondria.iloc[index]):
                reasons.append("mitochondrial_assignment")
            if bool(chloroplast.iloc[index]):
                reasons.append("chloroplast_assignment")
            if bool(explicit_eukaryote.iloc[index]):
                reasons.append("explicitly_eukaryotic_assignment")
            if bool(invalid_sequence.iloc[index]):
                reasons.append("invalid_sequence_characters")
        exclusion_reason.append(";".join(reasons))

    domain = taxonomy.str.extract(r"(?:^|;\s*)(d__[^;]+)", expand=False).fillna("unresolved")
    audit = pd.DataFrame(
        {
            "ASV_ID": averaged["id"].astype(str),
            "Sequence_length": sequences.str.len(),
            "Taxon": taxonomy,
            "Domain": domain,
            "Total_abundance_across_168_profiles": biological_total,
            "Prevalence_across_168_profiles": biological_prevalence,
            "Positive_in_at_least_one_profile": positive,
            "Mitochondrial_assignment": mitochondria,
            "Chloroplast_assignment": chloroplast,
            "Explicitly_eukaryotic_assignment": explicit_eukaryote,
            "Invalid_sequence_characters": invalid_sequence,
            "Included_in_23896_prokaryotic_matrix": retained,
            "Exclusion_reason": exclusion_reason,
        }
    )

    retained_counts = counts.loc[retained]
    archaeal = domain.str.casefold().eq("d__archaea")
    retained_total_reads = float(retained_counts.to_numpy().sum())
    retained_archaeal_reads = float(counts.loc[retained & archaeal].to_numpy().sum())
    archaeal_fraction = retained_archaeal_reads / retained_total_reads

    summary_rows = [
        ("raw_ASV_records", len(averaged), "Starting feature count"),
        ("ASVs_absent_from_all_168_profiles", int((~positive).sum()), "Removed before biological analysis"),
        ("ASVs_positive_in_at_least_one_profile", int(positive.sum()), "29,461 - 4,242"),
        ("positive_mitochondrial_ASVs", int((positive & mitochondria).sum()), "Excluded"),
        ("positive_chloroplast_ASVs", int((positive & chloroplast).sum()), "Excluded"),
        ("positive_explicitly_eukaryotic_ASVs", int((positive & explicit_eukaryote).sum()), "Excluded"),
        ("positive_organelle_or_eukaryotic_ASVs", int((positive & organelle_or_eukaryote).sum()), "728 + 590 + 5"),
        ("invalid_sequence_ASVs", int(invalid_sequence.sum()), "Expected 0"),
        ("retained_prokaryotic_ASVs", int(retained.sum()), "25,219 - 1,323"),
        ("retained_archaeal_ASVs", int((retained & archaeal).sum()), "Included in prokaryotic matrix"),
        ("retained_archaeal_read_fraction", archaeal_fraction, "Fraction, not percent"),
        ("retained_archaeal_read_percent", archaeal_fraction * 100.0, "Matches 0.0073% in manuscript"),
        ("source_biological_libraries", len(biological_libraries), "Formal source-library cohort"),
        ("canonical_profiles", len(canonical_profiles), "Pond-month-layer profiles"),
        ("maximum_aggregation_cell_difference", global_max_difference, "Expected 0"),
        ("nonmatching_aggregation_cells", int(aggregation_check["number_of_nonmatching_ASV_cells"].sum()), "Expected 0"),
    ]
    summary = pd.DataFrame(summary_rows, columns=["metric", "value", "interpretation"])

    retained_ids = audit.loc[
        audit["Included_in_23896_prokaryotic_matrix"],
        ["ASV_ID", "Taxon"],
    ]

    audit.to_csv(
        outdir / "asv_filter_audit_29461.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    retained_ids.to_csv(outdir / "retained_prokaryotic_asv_ids_23896.tsv", sep="\t", index=False)
    summary.to_csv(outdir / "asv_filter_summary.tsv", sep="\t", index=False)
    aggregation_check.to_csv(outdir / "profile_aggregation_check_168.tsv", sep="\t", index=False)

    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
