#!/usr/bin/env python3
"""Strict structural and numerical validation of the completed PICRUSt2 run."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_numeric_table(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", compression="infer", index_col=0)
    table.index = table.index.astype(str)
    table.columns = table.columns.astype(str)
    if table.index.duplicated().any() or table.columns.duplicated().any():
        raise ValueError(f"Duplicated row or column identifiers in {path}.")
    numeric = table.apply(pd.to_numeric, errors="raise")
    values = numeric.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError(f"Non-finite values in {path}.")
    if (values < 0).any():
        raise ValueError(f"Negative values in {path}.")
    return numeric


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--picrust", required=True, type=Path)
    parser.add_argument("--abundance", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    marker_path = args.picrust / "combined_marker_predicted_and_nsti.tsv.gz"
    ec_prediction_path = args.picrust / "combined_EC_predicted.tsv.gz"
    normalized_path = args.picrust / "EC_metagenome_out" / "seqtab_norm.tsv.gz"
    weighted_path = args.picrust / "EC_metagenome_out" / "weighted_nsti.tsv.gz"
    metagenome_path = (
        args.picrust / "EC_metagenome_out" / "pred_metagenome_unstrat.tsv.gz"
    )
    pathway_path = args.picrust / "pathways_out" / "path_abun_unstrat.tsv.gz"
    required_paths = [
        marker_path,
        ec_prediction_path,
        normalized_path,
        weighted_path,
        metagenome_path,
        pathway_path,
    ]
    for path in required_paths:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)

    abundance = read_numeric_table(args.abundance)
    marker = pd.read_csv(marker_path, sep="\t", compression="infer", index_col=0)
    marker.index = marker.index.astype(str)
    if marker.index.duplicated().any():
        raise ValueError("Duplicated ASVs in the combined marker table.")
    required_marker_columns = {
        "16S_rRNA_Count",
        "metadata_NSTI",
        "best_domain",
        "closest_reference_genome",
    }
    missing_marker_columns = required_marker_columns - set(marker.columns)
    if missing_marker_columns:
        raise ValueError(
            f"Marker table is missing columns: {sorted(missing_marker_columns)}"
        )
    marker_nsti = pd.to_numeric(marker["metadata_NSTI"], errors="raise")
    if not np.isfinite(marker_nsti).all() or (marker_nsti < 0).any():
        raise ValueError("Invalid NSTI values.")
    marker_copy = pd.to_numeric(marker["16S_rRNA_Count"], errors="raise")
    if not np.isfinite(marker_copy).all() or (marker_copy <= 0).any():
        raise ValueError("Invalid predicted 16S marker copy numbers.")

    sample_ids = abundance.columns.astype(str)
    if len(sample_ids) != 168 or len(set(sample_ids)) != 168:
        raise ValueError("Expected 168 unique formal sample columns.")
    if len(abundance) != 23896:
        raise ValueError(f"Expected 23,896 input ASVs, found {len(abundance):,}.")
    if not set(marker.index).issubset(abundance.index):
        raise ValueError("Marker table contains ASVs absent from the abundance input.")

    normalized = read_numeric_table(normalized_path)
    expected_retained_ids = marker.index[marker_nsti.le(2)]
    if set(normalized.index) != set(expected_retained_ids):
        raise ValueError(
            "Normalized ASV table does not exactly equal marker-assigned ASVs "
            "with NSTI ≤ 2."
        )
    if set(normalized.columns) != set(sample_ids):
        raise ValueError("Normalized ASV sample identifiers do not match the input.")

    weighted = read_numeric_table(weighted_path)
    if set(weighted.index) != set(sample_ids) or weighted.shape[1] != 1:
        raise ValueError("Weighted-NSTI output does not contain exactly 168 samples.")

    ec_metagenome = read_numeric_table(metagenome_path)
    pathways = read_numeric_table(pathway_path)
    for label, table in [
        ("EC metagenome", ec_metagenome),
        ("MetaCyc pathways", pathways),
    ]:
        if set(table.columns) != set(sample_ids):
            raise ValueError(f"{label} sample identifiers do not match the input.")
        if (table.sum(axis=0) <= 0).any():
            raise ValueError(f"{label} contains a zero-sum sample.")

    raw_read_total = abundance.sum(axis=0)
    retained_read_total = abundance.loc[normalized.index].sum(axis=0)
    retained_fraction = retained_read_total / raw_read_total
    if retained_fraction.isna().any() or not retained_fraction.between(0, 1).all():
        raise ValueError("Invalid retained-read fractions.")

    missing_assignment_ids = sorted(set(abundance.index) - set(marker.index))
    summary = {
        "status": "passed",
        "input_ASVs": int(len(abundance)),
        "formal_samples": int(len(sample_ids)),
        "ASVs_with_domain_assignment": int(len(marker)),
        "ASVs_without_domain_assignment": int(len(missing_assignment_ids)),
        "ASV_IDs_without_domain_assignment": missing_assignment_ids,
        "ASVs_above_max_NSTI_2": int(marker_nsti.gt(2).sum()),
        "ASVs_retained_for_metagenome_prediction": int(len(normalized)),
        "predicted_EC_families": int(len(ec_metagenome)),
        "predicted_MetaCyc_pathways": int(len(pathways)),
        "weighted_NSTI": {
            "minimum": float(weighted.iloc[:, 0].min()),
            "q1": float(weighted.iloc[:, 0].quantile(0.25)),
            "median": float(weighted.iloc[:, 0].median()),
            "q3": float(weighted.iloc[:, 0].quantile(0.75)),
            "maximum": float(weighted.iloc[:, 0].max()),
        },
        "retained_input_read_fraction": {
            "minimum": float(retained_fraction.min()),
            "q1": float(retained_fraction.quantile(0.25)),
            "median": float(retained_fraction.median()),
            "q3": float(retained_fraction.quantile(0.75)),
            "maximum": float(retained_fraction.max()),
        },
        "output_sha256": {
            path.name: sha256_file(path) for path in required_paths
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
