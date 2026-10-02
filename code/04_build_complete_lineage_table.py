#!/usr/bin/env python3
"""Build Additional file 1 from the strict 168-profile ASV matrix."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


SEASONS = ("Spring", "Summer", "Autumn", "Winter")
EXCLUDED_TAXA = r"mitochond|chloroplast|d__eukary"


def rank_value(path: str, prefix: str) -> str:
    for item in str(path).split(";"):
        item = item.strip()
        if item.startswith(prefix):
            return item[len(prefix) :].strip()
    return ""


def valid_genus(path: str) -> bool:
    genus = rank_value(path, "g__")
    return bool(genus) and genus.casefold() not in {
        "unclassified", "uncultured", "unknown", "na", "none"
    }


def lineage_key(path: str) -> str:
    return "; ".join(str(path).split(";")[:6]).strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--asv", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()

    asv = pd.read_csv(args.asv, low_memory=False)
    metadata = pd.read_csv(args.metadata)
    samples = metadata["SampleID"].astype(str).tolist()
    if len(samples) != 168 or metadata["PondMonth"].nunique() != 93:
        raise ValueError("Expected 168 profiles and 93 pond-months.")
    sample_columns = [c for c in asv.columns if c not in ("id", "Sequence", "Taxon")]
    if sample_columns != samples:
        raise ValueError("The ASV table must contain exactly the 168 metadata samples in order.")

    taxonomy = asv["Taxon"].fillna("").astype(str)
    positive = asv[samples].fillna(0).sum(axis=1) > 0
    excluded = taxonomy.str.contains(EXCLUDED_TAXA, case=False, regex=True)
    kept = asv.loc[positive & ~excluded].copy()
    if len(kept) != 23896:
        raise ValueError(f"Expected 23,896 retained ASVs; found {len(kept):,}.")

    counts = kept[samples].T.astype(float)
    relative = counts.div(counts.sum(axis=1), axis=0)
    resolved = kept["Taxon"].map(valid_genus)
    relative = relative.loc[:, resolved.to_numpy()]
    relative.columns = kept.loc[resolved, "Taxon"].map(lineage_key).to_numpy()
    relative = relative.T.groupby(level=0).sum().T
    if relative.shape[1] != 2673:
        raise ValueError(f"Expected 2,673 resolved genus lineages; found {relative.shape[1]:,}.")

    meta = metadata.set_index("SampleID").loc[relative.index]
    pond_month = relative.assign(PondMonth=meta["PondMonth"]).groupby("PondMonth").mean()
    pm_meta = (
        metadata[["PondMonth", "Pond", "Date", "Season_verified"]]
        .drop_duplicates("PondMonth")
        .set_index("PondMonth")
        .loc[pond_month.index]
    )
    grow_mask = pm_meta["Pond"].ne("T1") & pm_meta["Season_verified"].notna()
    grow = pond_month.loc[grow_mask]
    grow_season = pm_meta.loc[grow.index, "Season_verified"]
    if len(grow) != 70:
        raise ValueError(f"Expected 70 grow-out pond-months; found {len(grow)}.")

    rows: list[dict[str, object]] = []
    for lineage in pond_month.columns:
        path = str(lineage)
        all_values = pond_month[lineage]
        grow_values = grow[lineage]
        row: dict[str, object] = {
            "Lineage": path,
            "Phylum": rank_value(path, "p__"),
            "Class": rank_value(path, "c__"),
            "Order": rank_value(path, "o__"),
            "Family": rank_value(path, "f__"),
            "Genus": rank_value(path, "g__"),
            "All_pond_month_mean_relative_abundance_percent": 100 * all_values.mean(),
            "All_pond_month_prevalence_count": int((all_values > 0).sum()),
            "All_pond_month_prevalence_percent": 100 * (all_values > 0).mean(),
            "Persistent_at_least_84_of_93": int((all_values > 0).sum()) >= 84,
            "Ubiquitous_93_of_93": int((all_values > 0).sum()) == 93,
            "Grow_out_70_mean_relative_abundance_percent": 100 * grow_values.mean(),
            "Grow_out_70_prevalence_percent": 100 * (grow_values > 0).mean(),
        }
        abundance_values: list[float] = []
        prevalence_values: list[float] = []
        for season in SEASONS:
            values = grow.loc[grow_season.eq(season), lineage]
            abundance = 100 * values.mean()
            prevalence = 100 * (values > 0).mean()
            row[f"{season}_mean_relative_abundance_percent"] = abundance
            row[f"{season}_prevalence_percent"] = prevalence
            abundance_values.append(abundance)
            prevalence_values.append(prevalence)
        row["Seasonal_abundance_range_percentage_points"] = max(abundance_values) - min(abundance_values)
        row["Seasonal_prevalence_range_percentage_points"] = max(prevalence_values) - min(prevalence_values)
        row["Highest_prevalence_season"] = SEASONS[int(np.argmax(prevalence_values))]
        rows.append(row)

    output = pd.DataFrame(rows)
    eligible = output[
        output["Grow_out_70_prevalence_percent"].ge(10)
        & output["Grow_out_70_mean_relative_abundance_percent"].ge(0.01)
    ].sort_values(
        [
            "Seasonal_prevalence_range_percentage_points",
            "Grow_out_70_mean_relative_abundance_percent",
            "Lineage",
        ],
        ascending=[False, False, True],
    )
    selected = set(eligible.head(30)["Lineage"])
    output["Selected_top30_seasonal_prevalence_range"] = output["Lineage"].isin(selected)
    output = output.sort_values(
        ["All_pond_month_prevalence_count", "All_pond_month_mean_relative_abundance_percent", "Lineage"],
        ascending=[False, False, True],
    ).reset_index(drop=True)

    prevalence = output["All_pond_month_prevalence_percent"].to_numpy(float)
    abundance = output["All_pond_month_mean_relative_abundance_percent"].to_numpy(float)
    correlation = spearmanr(np.log10(abundance), prevalence)
    persistent = output["Persistent_at_least_84_of_93"]
    selected_mask = output["Selected_top30_seasonal_prevalence_range"]
    summary = {
        "resolved_genus_level_lineages": int(len(output)),
        "resolved_lineage_mean_relative_abundance_percent": float(abundance.sum()),
        "ubiquitous_lineages_93_of_93": int(output["Ubiquitous_93_of_93"].sum()),
        "persistent_lineages_at_least_84_of_93": int(persistent.sum()),
        "persistent_combined_mean_relative_abundance_percent": float(abundance[persistent].sum()),
        "abundance_prevalence_spearman_rho": float(correlation.statistic),
        "selected_top30_combined_grow_out_mean_relative_abundance_percent": float(
            output.loc[selected_mask, "Grow_out_70_mean_relative_abundance_percent"].sum()
        ),
        "selected_top30_maximum_individual_grow_out_mean_percent": float(
            output.loc[selected_mask, "Grow_out_70_mean_relative_abundance_percent"].max()
        ),
        "selected_top30_median_seasonal_abundance_range_percentage_points": float(
            output.loc[selected_mask, "Seasonal_abundance_range_percentage_points"].median()
        ),
        "selected_top30_median_seasonal_prevalence_range_percentage_points": float(
            output.loc[selected_mask, "Seasonal_prevalence_range_percentage_points"].median()
        ),
        "persistent_median_seasonal_abundance_range_percentage_points": float(
            output.loc[persistent, "Seasonal_abundance_range_percentage_points"].median()
        ),
        "persistent_median_seasonal_prevalence_range_percentage_points": float(
            output.loc[persistent, "Seasonal_prevalence_range_percentage_points"].median()
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)
    args.summary.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
