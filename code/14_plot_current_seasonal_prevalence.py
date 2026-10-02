#!/usr/bin/env python3
"""Rebuild the current seasonal-community and prevalence figure (Fig. 4)."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial.distance import pdist, squareform


ORGANELLE_EUKARYOTE_PATTERN = r"mitochond|chloroplast|d__eukary"
SEASONS = ["Spring", "Summer", "Autumn", "Winter"]


def relative_abundance_rows(counts: np.ndarray) -> np.ndarray:
    totals = counts.sum(axis=1, keepdims=True)
    return np.divide(
        counts,
        totals,
        out=np.zeros_like(counts, dtype=float),
        where=totals > 0,
    )


def pcoa(distance_matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = distance_matrix.shape[0]
    center = np.eye(n) - np.ones((n, n)) / n
    gower = -0.5 * center @ (distance_matrix**2) @ center
    eigenvalues, eigenvectors = np.linalg.eigh(gower)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    positive = eigenvalues > 0
    coordinates = eigenvectors[:, positive] * np.sqrt(eigenvalues[positive])
    explained = 100 * eigenvalues[positive] / eigenvalues[positive].sum()
    return coordinates, explained


def rank_value(path: str, prefix: str) -> str:
    for item in str(path).split(";"):
        item = item.strip()
        if item.startswith(prefix):
            return item[len(prefix) :].strip()
    return ""


def valid_genus(path: str) -> bool:
    genus = rank_value(path, "g__")
    return bool(genus) and genus.casefold() not in {
        "unclassified",
        "uncultured",
        "unknown",
        "na",
        "none",
    }


def lineage_key(path: str) -> str:
    return "; ".join(str(path).split(";")[:6]).strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--asv", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--permanova", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--selected-table", type=Path)
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata)
    asv = pd.read_csv(args.asv, low_memory=False)
    samples = metadata["SampleID"].astype(str).tolist()
    taxonomy = asv["Taxon"].fillna("").astype(str)
    positive = asv[samples].fillna(0).sum(axis=1) > 0
    excluded = taxonomy.str.contains(
        ORGANELLE_EUKARYOTE_PATTERN, case=False, regex=True
    )
    keep = positive & ~excluded
    kept = asv.loc[keep].copy()
    relative = relative_abundance_rows(kept[samples].T.to_numpy(dtype=float))

    sample_relative = pd.DataFrame(relative, index=samples)
    sample_relative["PondMonth"] = metadata.set_index("SampleID").loc[
        samples, "PondMonth"
    ].to_numpy()
    pond_month_relative = sample_relative.groupby("PondMonth", sort=False).mean()
    pond_month_metadata = (
        metadata[["PondMonth", "Pond", "Date", "Season_verified"]]
        .drop_duplicates("PondMonth")
        .set_index("PondMonth")
    )
    cohort = pond_month_metadata["Pond"].isin(["T2", "T3", "T4"]) & pond_month_metadata[
        "Season_verified"
    ].notna()
    cohort_metadata = pond_month_metadata.loc[cohort].copy()
    cohort_relative = pond_month_relative.loc[cohort_metadata.index]
    if len(cohort_metadata) != 70:
        raise ValueError(f"Expected 70 grow-out pond-months; found {len(cohort_metadata)}")

    distance = squareform(pdist(cohort_relative.to_numpy(), metric="braycurtis"))
    coordinates, explained = pcoa(distance)
    ordination = cohort_metadata.reset_index().rename(
        columns={"Season_verified": "Season"}
    )
    ordination["PCoA1"] = coordinates[:, 0]
    ordination["PCoA2"] = coordinates[:, 1]

    resolved = kept["Taxon"].map(valid_genus)
    resolved_relative = pd.DataFrame(
        cohort_relative.loc[:, resolved.to_numpy()].to_numpy(),
        index=cohort_relative.index,
        columns=kept.loc[resolved, "Taxon"].map(lineage_key).to_numpy(),
    )
    resolved_relative = resolved_relative.T.groupby(level=0).sum().T
    seasons = cohort_metadata.loc[resolved_relative.index, "Season_verified"]
    overall_prevalence = 100 * (resolved_relative > 0).mean(axis=0)
    overall_mean = 100 * resolved_relative.mean(axis=0)

    rows: list[dict[str, object]] = []
    for lineage in resolved_relative.columns:
        path = str(lineage)
        row: dict[str, object] = {
            "Lineage": path,
            "Genus": rank_value(path, "g__"),
            "Phylum": rank_value(path, "p__"),
            "Family": rank_value(path, "f__"),
            "Grow_out_prevalence_percent": overall_prevalence[lineage],
            "Grow_out_mean_relative_abundance_percent": overall_mean[lineage],
        }
        prevalence_values = []
        for season in SEASONS:
            values = resolved_relative.loc[seasons.eq(season), lineage]
            abundance = 100 * values.mean()
            prevalence = 100 * (values > 0).mean()
            row[f"{season}_mean_percent"] = abundance
            row[f"{season}_prevalence_percent"] = prevalence
            prevalence_values.append(prevalence)
        row["Seasonal_prevalence_range_percentage_points"] = max(
            prevalence_values
        ) - min(prevalence_values)
        row["Highest_prevalence_season"] = SEASONS[
            int(np.argmax(prevalence_values))
        ]
        rows.append(row)
    seasonal = pd.DataFrame(rows)
    selected = (
        seasonal.loc[
            seasonal["Grow_out_prevalence_percent"].ge(10)
            & seasonal["Grow_out_mean_relative_abundance_percent"].ge(0.01)
        ]
        .sort_values(
            [
                "Seasonal_prevalence_range_percentage_points",
                "Grow_out_mean_relative_abundance_percent",
                "Lineage",
            ],
            ascending=[False, False, True],
        )
        .head(30)
        .reset_index(drop=True)
    )
    if len(selected) != 30:
        raise ValueError(f"Expected 30 selected lineages; found {len(selected)}")
    if args.selected_table:
        args.selected_table.resolve().parent.mkdir(parents=True, exist_ok=True)
        selected.to_csv(args.selected_table, index=False)

    displayed = selected.head(20).sort_values(
        "Seasonal_prevalence_range_percentage_points"
    )
    prevalence_columns = [f"{season}_prevalence_percent" for season in SEASONS]
    raw = displayed[prevalence_columns].to_numpy(dtype=float)
    row_mean = raw.mean(axis=1, keepdims=True)
    row_sd = raw.std(axis=1, ddof=0, keepdims=True)
    standardized = np.divide(
        raw - row_mean,
        row_sd,
        out=np.zeros_like(raw),
        where=row_sd > 0,
    )
    labels = displayed.apply(
        lambda row: f"{row['Genus']} | {row['Phylum']} | {row['Family']}", axis=1
    )

    permanova = pd.read_csv(args.permanova)
    season_test = permanova.loc[
        permanova["Dataset"].eq("Grow-out ponds")
        & permanova["Factor"].eq("Season")
    ].iloc[0]
    palette = {
        "Spring": "#4C78A8",
        "Summer": "#F58518",
        "Autumn": "#54A24B",
        "Winter": "#E45756",
    }

    fig = plt.figure(figsize=(12.0, 13.4))
    grid = fig.add_gridspec(2, 1, height_ratios=[0.92, 1.52], hspace=0.25)
    axis = fig.add_subplot(grid[0])
    for season, group in ordination.groupby("Season"):
        axis.scatter(
            group["PCoA1"],
            group["PCoA2"],
            s=45,
            alpha=0.82,
            color=palette[season],
            label=season,
        )
    axis.axhline(0, color="#AAAAAA", linewidth=0.7)
    axis.axvline(0, color="#AAAAAA", linewidth=0.7)
    axis.set_xlabel(f"PCoA1 ({explained[0]:.1f}%)")
    axis.set_ylabel(f"PCoA2 ({explained[1]:.1f}%)")
    axis.set_title("(a) Seasonal community structure", loc="left", fontweight="bold")
    axis.text(
        0.98,
        0.03,
        f"Marginal R² = {100 * season_test['Marginal_R2']:.2f}%\n"
        f"restricted p = {season_test['p_value']:.4f}; n = {int(season_test['n'])}",
        transform=axis.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
    )
    axis.legend(title="Season", frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1))
    axis.spines[["top", "right"]].set_visible(False)

    subgrid = grid[1].subgridspec(1, 2, width_ratios=[1.0, 0.28], wspace=0.08)
    heat_axis = fig.add_subplot(subgrid[0, 0])
    annotations = np.vectorize(lambda value: f"{value:.1f}")(raw)
    sns.heatmap(
        standardized,
        cmap="RdBu_r",
        center=0,
        vmin=-2,
        vmax=2,
        annot=annotations,
        fmt="",
        annot_kws={"fontsize": 7.3},
        xticklabels=SEASONS,
        yticklabels=labels,
        linewidths=0.3,
        linecolor="white",
        cbar_kws={"label": "Row-standardized prevalence", "shrink": 0.72},
        ax=heat_axis,
    )
    heat_axis.set_xlabel("Season (cell labels: pond-month prevalence, %)")
    heat_axis.set_ylabel("")
    heat_axis.set_title(
        "(b) Descriptive seasonal prevalence profiles",
        loc="left",
        fontweight="bold",
    )

    range_axis = fig.add_subplot(subgrid[0, 1], sharey=heat_axis)
    ranges = displayed["Seasonal_prevalence_range_percentage_points"].to_numpy()
    range_axis.barh(np.arange(len(displayed)), ranges, color="#7186A0", height=0.72)
    for y, value in enumerate(ranges):
        range_axis.text(value + 1.0, y, f"{value:.1f}", va="center", fontsize=7.2)
    range_axis.set_xlabel("Prevalence range\n(percentage points)")
    range_axis.tick_params(axis="y", left=False, labelleft=False)
    range_axis.spines[["top", "right"]].set_visible(False)
    range_axis.set_xlim(0, max(ranges) * 1.24)
    fig.subplots_adjust(left=0.32, right=0.96, top=0.97, bottom=0.06)

    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=400, bbox_inches="tight", transparent=True)
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight", transparent=True)
    plt.close(fig)


if __name__ == "__main__":
    main()
