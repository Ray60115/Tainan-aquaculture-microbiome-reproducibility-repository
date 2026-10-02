#!/usr/bin/env python3
"""Compare original and prokaryote-filtered dependence-aware results."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


TABLE_SPECS = {
    "shannon_dependence_aware_tests.csv": {
        "keys": ["Analysis"],
        "values": ["n", "Effect", "p_value", "q_value"],
    },
    "permanova_dependence_aware_results.csv": {
        "keys": ["Dataset", "Factor"],
        "values": [
            "n",
            "F",
            "Marginal_R2",
            "Partial_R2",
            "p_value",
            "q_value",
        ],
    },
    "environmental_vector_fits_dependence_aware.csv": {
        "keys": ["Environmental_variable"],
        "values": ["Vector_fit_R2", "p_value", "q_value"],
    },
    "turnover_dependence_aware_test.csv": {
        "keys": ["Analysis"],
        "values": ["n", "F", "partial_R2", "p_value"],
    },
    "blockwise_unique_R2_pond_month.csv": {
        "keys": ["Dataset", "Predictor_block"],
        "values": ["Full_model_R2", "Unique_R2"],
    },
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original", required=True, type=Path)
    parser.add_argument("--filtered", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    outputs = []
    for filename, spec in TABLE_SPECS.items():
        original = pd.read_csv(args.original / filename)
        filtered = pd.read_csv(args.filtered / filename)
        available_values = [
            c for c in spec["values"] if c in original.columns and c in filtered.columns
        ]
        left = original[spec["keys"] + available_values]
        right = filtered[spec["keys"] + available_values]
        merged = left.merge(
            right,
            on=spec["keys"],
            how="outer",
            suffixes=("_original", "_prokaryotic"),
            validate="one_to_one",
        )
        for value in available_values:
            old = pd.to_numeric(merged[f"{value}_original"], errors="coerce")
            new = pd.to_numeric(merged[f"{value}_prokaryotic"], errors="coerce")
            merged[f"{value}_change"] = new - old
        merged.insert(0, "Source_table", filename)
        outputs.append(merged)

    combined = pd.concat(outputs, ignore_index=True, sort=False)
    combined.to_csv(
        args.outdir / "organelle_filter_result_comparison.csv", index=False
    )

    original_quality = json.loads(
        (args.original.parent / "analysis_quality_report.json").read_text(
            encoding="utf-8"
        )
    )
    filtered_quality = json.loads(
        (args.filtered.parent / "analysis_quality_report.json").read_text(
            encoding="utf-8"
        )
    )

    significance_columns = [
        c
        for c in combined.columns
        if c.endswith("_original")
        and (c.startswith("p_value") or c.startswith("q_value"))
    ]
    classification_changes = []
    for original_col in significance_columns:
        base = original_col.removesuffix("_original")
        filtered_col = f"{base}_prokaryotic"
        if filtered_col not in combined.columns:
            continue
        old = pd.to_numeric(combined[original_col], errors="coerce")
        new = pd.to_numeric(combined[filtered_col], errors="coerce")
        changed = old.lt(0.05).ne(new.lt(0.05)) & old.notna() & new.notna()
        if changed.any():
            classification_changes.extend(
                combined.loc[changed, ["Source_table"] + [
                    c for c in ["Analysis", "Dataset", "Factor", "Environmental_variable", "Predictor_block"]
                    if c in combined.columns
                ]].to_dict("records")
            )

    numeric_changes = [
        pd.to_numeric(combined[c], errors="coerce").abs().max()
        for c in combined.columns
        if c.endswith("_change")
    ]
    numeric_changes = [x for x in numeric_changes if np.isfinite(x)]
    summary = {
        "organelle_or_eukaryotic_read_fraction_removed": filtered_quality[
            "feature_filter"
        ]["Organelle_or_eukaryotic_read_fraction_removed"],
        "original_PCoA1_explained_percent": original_quality["pcoa1_explained"],
        "prokaryotic_PCoA1_explained_percent": filtered_quality[
            "pcoa1_explained"
        ],
        "original_PCoA2_explained_percent": original_quality["pcoa2_explained"],
        "prokaryotic_PCoA2_explained_percent": filtered_quality[
            "pcoa2_explained"
        ],
        "statistical_significance_classification_changes_at_0.05": (
            classification_changes
        ),
        "maximum_absolute_numeric_change_across_compared_fields": float(
            max(numeric_changes)
        )
        if numeric_changes
        else None,
        "interpretation": (
            "Organelle removal changed numerical estimates but did not change "
            "the significance classification of the dependence-aware tests "
            "listed in the comparison."
            if not classification_changes
            else "At least one significance classification changed."
        ),
    }
    (args.outdir / "organelle_filter_comparison_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
