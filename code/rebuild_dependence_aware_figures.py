#!/usr/bin/env python3
"""Rebuild publication figures from saved dependence-aware result tables."""

from pathlib import Path
import json
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-tainan")

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Ellipse
from scipy.spatial.distance import pdist, squareform

from reanalyze_dependence_aware import CLIMATE_FILE, pcoa, relative_abundance_rows


CODE = Path(__file__).resolve().parent
ROOT = CODE.parent
RESULTS = Path(os.environ.get("TAINAN_RESULTS_DIR", ROOT / "reproduced_results" / "taxonomy"))
TABLES = RESULTS / "tables"
FIGURES = RESULTS / "figures"
DATA = ROOT / "data" / "primary"


def load_community() -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    meta = pd.read_csv(TABLES / "validated_sample_metadata_with_shannon.csv")
    asv = pd.read_csv(DATA / "asv_table_168_profiles.csv.gz")
    sample_columns = [c for c in asv.columns if c in set(meta["SampleID"])]
    counts = asv[sample_columns].T
    counts.index.name = "SampleID"
    counts = counts.loc[meta["SampleID"]]
    rel = relative_abundance_rows(counts.to_numpy(dtype=float))
    dist = squareform(pdist(rel, metric="braycurtis"))
    return meta, rel, dist


def save_figure(fig: plt.Figure, filename: str) -> None:
    fig.savefig(FIGURES / filename, dpi=400, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def remove_top_right_spines(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def add_covariance_ellipse(
    ax: plt.Axes,
    x: pd.Series,
    y: pd.Series,
    color: str,
) -> None:
    """Add a descriptive 95% bivariate covariance ellipse."""
    xy = np.column_stack([x.to_numpy(dtype=float), y.to_numpy(dtype=float)])
    if len(xy) < 3 or not np.isfinite(xy).all():
        return
    covariance = np.cov(xy, rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    order = eigenvalues.argsort()[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    if np.any(eigenvalues <= 0):
        return
    angle = np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))
    # sqrt(chi-square 0.95 quantile with 2 df) = sqrt(5.991)
    width, height = 2 * np.sqrt(eigenvalues * 5.991)
    ellipse = Ellipse(
        xy=xy.mean(axis=0),
        width=width,
        height=height,
        angle=angle,
        facecolor=color,
        edgecolor=color,
        alpha=0.12,
        linewidth=1.2,
        zorder=0,
    )
    ax.add_patch(ellipse)


def format_abundance(value: float) -> str:
    if value == 0:
        return "0"
    if value < 0.01:
        return "<0.01"
    if value < 1:
        return f"{value:.2f}"
    return f"{value:.1f}"


def format_p_value(value: float) -> str:
    """Format p-values consistently with the manuscript and result tables."""
    if value < 0.0001:
        exponent = int(np.floor(np.log10(value)))
        coefficient = value / (10 ** exponent)
        superscript = str(exponent).translate(
            str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
        )
        return f"{coefficient:.2f} × 10{superscript}"
    return f"{value:.4f}"


def rebuild_figure2() -> None:
    daily = pd.read_csv(CLIMATE_FILE)
    daily["FullDate"] = pd.to_datetime(daily["FullDate"])
    pm = pd.read_csv(TABLES / "pond_month_analysis_metadata.csv")
    event_columns = [
        "Heatwave_Tmax_q90_3day",
        "Typhoon",
        "ColdSurge",
        "LowPressure",
    ]
    event_months = (
        pm[["Date", *event_columns]]
        .drop_duplicates("Date")
        .sort_values("Date")
        .copy()
    )
    event_months["Date_dt"] = pd.to_datetime(
        event_months["Date"].astype(str), format="%Y%m"
    )
    quality = json.loads(
        (RESULTS / "analysis_quality_report.json").read_text(encoding="utf-8")
    )
    monthly_low_pressure_threshold = float(
        quality["low_pressure_threshold_hPa"]
    )

    fig, axes = plt.subplots(
        4, 1, figsize=(12.2, 11.0), sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.0, 1.0, 0.82]},
    )
    start = pd.Timestamp("2022-08-15")
    end = pd.Timestamp("2024-12-31")

    axes[0].plot(daily["FullDate"], daily["Tmax"], color="#C65D3B", lw=1.0)
    heat_days = daily["Tmax_gt_q90"].eq(1)
    axes[0].scatter(
        daily.loc[heat_days, "FullDate"],
        daily.loc[heat_days, "Tmax"],
        s=11,
        color="#9E1B1B",
        zorder=3,
        label="Daily Tmax ≥ 34.8 °C",
    )
    axes[0].axhline(34.8, color="#9E1B1B", ls="--", lw=0.9)
    axes[0].set_ylabel("Daily maximum\nair temperature (°C)")
    axes[0].set_title("(a) Maximum air temperature", loc="left", fontweight="bold")
    axes[0].legend(frameon=False, loc="upper left", fontsize=8.5)

    daily_pressure_q10 = float(daily["PressureMin"].quantile(0.10))
    axes[1].plot(
        daily["FullDate"], daily["PressureMin"], color="#3A7D65", lw=1.0
    )
    low_days = daily["PressureMin"].le(daily_pressure_q10)
    axes[1].scatter(
        daily.loc[low_days, "FullDate"],
        daily.loc[low_days, "PressureMin"],
        s=10,
        color="#145A43",
        zorder=3,
        label=f"Lowest daily decile (≤ {daily_pressure_q10:.1f} hPa)",
    )
    axes[1].axhline(daily_pressure_q10, color="#145A43", ls="--", lw=0.9)
    axes[1].set_ylabel("Daily minimum\npressure (hPa)")
    axes[1].set_title("(b) Minimum atmospheric pressure", loc="left", fontweight="bold")
    axes[1].legend(frameon=False, loc="lower left", fontsize=8.5)

    axes[2].bar(
        daily["FullDate"], daily["Rainfall"], width=1.0,
        color="#4C78A8", alpha=0.78,
    )
    axes[2].set_ylabel("Daily rainfall\n(mm)")
    axes[2].set_title("(c) Rainfall", loc="left", fontweight="bold")

    event_specs = [
        ("Heatwave_Tmax_q90_3day", "Heat event", "#D9534F"),
        ("Typhoon", "Typhoon", "#3978B5"),
        ("ColdSurge", "Cold surge", "#7A4AA5"),
        (
            "LowPressure",
            f"Low pressure\n(monthly minimum ≤ {monthly_low_pressure_threshold:.1f} hPa)",
            "#17877C",
        ),
    ]
    for y, (column, label, color) in enumerate(event_specs):
        selected = event_months.loc[event_months[column].eq(1), "Date_dt"]
        axes[3].scatter(
            selected,
            np.full(len(selected), y),
            marker="s",
            s=100,
            color=color,
            edgecolor="white",
            linewidth=0.6,
            zorder=3,
        )
    axes[3].set_yticks(
        np.arange(len(event_specs)), [x[1] for x in event_specs]
    )
    axes[3].invert_yaxis()
    axes[3].set_title(
        "(d) Operational monthly disturbance classifications",
        loc="left",
        fontweight="bold",
    )
    axes[3].set_xlabel("Date")
    axes[3].grid(axis="x", alpha=0.18)

    for ax in axes:
        ax.set_xlim(start, end)
        remove_top_right_spines(ax)
    axes[3].xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    axes[3].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    plt.setp(axes[3].get_xticklabels(), rotation=45, ha="right")
    fig.tight_layout()
    save_figure(fig, "Fig1_meteorological_context_panels_b_to_f.png")


def rebuild_figure3() -> None:
    meta, _, dist = load_community()
    coords, _, explained = pcoa(dist)
    df = meta.copy()
    df["PCoA1"] = coords[:, 0]
    df["PCoA2"] = coords[:, 1]
    df["Culture category"] = np.where(df["Pond"].eq("T1"), "Nursery (T1)", "Grow-out (T2–T4)")

    pond_palette = {"T1": "#8E78B8", "T2": "#E78B93", "T3": "#F2AE3D", "T4": "#557EA6"}
    stage_palette = {"Nursery (T1)": "#8E78B8", "Grow-out (T2–T4)": "#5A9DB5"}
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 11.5))

    for pond, sub in df.groupby("Pond"):
        add_covariance_ellipse(
            axes[0, 0], sub["PCoA1"], sub["PCoA2"], pond_palette[pond]
        )
        axes[0, 0].scatter(sub["PCoA1"], sub["PCoA2"], s=28, alpha=0.78,
                           color=pond_palette[pond], label=pond)
    axes[0, 0].legend(frameon=False, ncol=2)
    axes[0, 0].set_title("(a) Community ordination by pond", loc="left", fontweight="bold")

    for stage, sub in df.groupby("Culture category"):
        add_covariance_ellipse(
            axes[0, 1], sub["PCoA1"], sub["PCoA2"], stage_palette[stage]
        )
        axes[0, 1].scatter(sub["PCoA1"], sub["PCoA2"], s=28, alpha=0.72,
                           color=stage_palette[stage], label=stage)
        axes[0, 1].scatter(sub["PCoA1"].mean(), sub["PCoA2"].mean(), s=120,
                           marker="X", edgecolor="#333333", linewidth=0.8,
                           color=stage_palette[stage])
    axes[0, 1].legend(frameon=False)
    axes[0, 1].set_title("(b) Community ordination by culture category", loc="left",
                         fontweight="bold")

    sns.boxplot(data=df, x="Pond", y="Shannon", hue="Pond",
                order=["T1", "T2", "T3", "T4"], palette=pond_palette,
                width=0.55, showfliers=False, legend=False, ax=axes[1, 0])
    sns.stripplot(data=df, x="Pond", y="Shannon", hue="Pond",
                  order=["T1", "T2", "T3", "T4"], palette=pond_palette,
                  size=3.7, alpha=0.65, legend=False, ax=axes[1, 0])
    axes[1, 0].set_title("(c) Shannon diversity by pond", loc="left", fontweight="bold")
    axes[1, 0].set_xlabel("")

    pairs = pd.read_csv(TABLES / "shannon_complete_layer_pairs.csv")
    layer_test = pd.read_csv(TABLES / "shannon_dependence_aware_tests.csv")
    layer_test = layer_test.loc[
        layer_test["Analysis"].eq("Sampling layer")
    ].iloc[0]
    long = pairs.melt(id_vars=["Pond", "Date", "PondMonth"], value_vars=["U", "D"],
                      var_name="Layer", value_name="Shannon")
    long["Layer"] = long["Layer"].map({"U": "Surface", "D": "Bottom"})
    for _, pair in long.groupby("PondMonth"):
        pair = pair.set_index("Layer").loc[["Surface", "Bottom"]].reset_index()
        axes[1, 1].plot(
            [0, 1], pair["Shannon"],
            color="#8E8E8E", alpha=0.12, lw=0.55, zorder=0,
        )
    sns.boxplot(data=long, x="Layer", y="Shannon", hue="Layer", order=["Surface", "Bottom"],
                palette={"Surface": "#75A8C2", "Bottom": "#B27A63"},
                width=0.48, showfliers=False, legend=False, ax=axes[1, 1])
    sns.stripplot(data=long, x="Layer", y="Shannon", order=["Surface", "Bottom"],
                  color="#3D3D3D", size=2.7, alpha=0.35, ax=axes[1, 1])
    axes[1, 1].text(0.98, 0.03,
                    "Wilcoxon signed-rank\n"
                    f"W = {int(round(layer_test['Statistic']))}; "
                    f"p = {format_p_value(float(layer_test['p_value']))}\n"
                    f"n = {int(layer_test['n'])} paired pond-months",
                    transform=axes[1, 1].transAxes, ha="right", va="bottom", fontsize=9)
    axes[1, 1].set_title("(d) Paired Shannon diversity by sampling layer", loc="left",
                         fontweight="bold")
    axes[1, 1].set_xlabel("")

    for ax in axes.flat[:2]:
        ax.axhline(0, color="#AAAAAA", lw=0.7)
        ax.axvline(0, color="#AAAAAA", lw=0.7)
        ax.set_xlabel(f"PCoA1 ({explained[0]:.1f}%)")
        ax.set_ylabel(f"PCoA2 ({explained[1]:.1f}%)")
        remove_top_right_spines(ax)
    for ax in axes.flat[2:]:
        ax.set_ylabel("Shannon diversity")
        ax.grid(axis="y", alpha=0.2)
        remove_top_right_spines(ax)
    fig.tight_layout()
    save_figure(fig, "Fig2_community_and_alpha_diversity.png")


def rebuild_figure4() -> None:
    stage = pd.read_csv(TABLES / "stage_taxa_descriptive.csv")
    by_pond = pd.read_csv(TABLES / "stage_taxa_by_pond.csv")
    df = stage.merge(by_pond, on="Representative taxon", how="left")
    fold_column = "log₂ fold change (nursery vs grow-out)"
    pond_columns = ["T1 (%)", "T2 (%)", "T3 (%)", "T4 (%)"]
    df["Phylum"] = (
        df["Full taxonomic classification"]
        .str.split(";")
        .str[0]
        .str.strip()
    )
    df = df.sort_values(fold_column).reset_index(drop=True)
    taxa = df["Representative taxon"].tolist()
    phyla = df["Phylum"].drop_duplicates().tolist()
    palette = dict(
        zip(phyla, sns.color_palette("colorblind", n_colors=len(phyla)))
    )

    fig = plt.figure(figsize=(14.0, 10.7))
    gs = fig.add_gridspec(
        1, 2, width_ratios=[1.22, 1.0], wspace=0.50,
    )
    ax = fig.add_subplot(gs[0, 0])
    y = np.arange(len(df))
    for i, row in df.iterrows():
        value = float(row[fold_column])
        color = palette[row["Phylum"]]
        ax.plot([0, value], [i, i], color=color, alpha=0.72, lw=1.5)
        ax.scatter(value, i, color=color, s=35, zorder=3)
    ax.axvline(0, color="#555555", lw=0.8)
    ax.set_yticks(y, taxa)
    ax.set_xlabel("log₂(mean abundance in T1 / mean abundance in T2–T4)")
    ax.set_title(
        "(a) Nursery–grow-out abundance contrasts",
        loc="left",
        fontweight="bold",
    )
    ax.grid(axis="x", alpha=0.16)
    remove_top_right_spines(ax)
    handles = [
        plt.Line2D(
            [0], [0], marker="o", linestyle="", color=palette[phylum],
            label=phylum, markersize=6,
        )
        for phylum in phyla
    ]
    fig.legend(
        handles=handles,
        title="Phylum",
        frameon=False,
        fontsize=8,
        title_fontsize=8.5,
        loc="lower center",
        bbox_to_anchor=(0.37, 0.015),
        ncol=3,
    )

    ax2 = fig.add_subplot(gs[0, 1])
    # Matplotlib places y = 0 at the bottom in panel (a), whereas seaborn places
    # matrix row 0 at the top.  Reverse the matrix rows so that every horizontal
    # position in panel (b) represents the taxon shown at the same height in
    # panel (a).
    heatmap_df = df.iloc[::-1].reset_index(drop=True)
    raw = heatmap_df[pond_columns].to_numpy(dtype=float)
    means = raw.mean(axis=1, keepdims=True)
    standard_deviations = raw.std(axis=1, ddof=0, keepdims=True)
    standardized = np.divide(
        raw - means,
        standard_deviations,
        out=np.zeros_like(raw),
        where=standard_deviations > 0,
    )
    annotation = np.vectorize(format_abundance)(raw)
    sns.heatmap(
        standardized,
        cmap="RdBu_r",
        center=0,
        vmin=-2,
        vmax=2,
        annot=annotation,
        fmt="",
        annot_kws={"fontsize": 7.6},
        xticklabels=["T1\n(nursery)", "T2", "T3", "T4"],
        yticklabels=False,
        linewidths=0.35,
        linecolor="white",
        cbar_kws={
            "label": "Row-standardized mean abundance",
            "shrink": 0.70,
            "pad": 0.03,
        },
        ax=ax2,
    )
    ax2.set_xlabel("Pond (cell labels: mean relative abundance, %)")
    ax2.set_title(
        "(b) Pond-specific abundance patterns",
        loc="left",
        fontweight="bold",
    )
    fig.subplots_adjust(left=0.28, right=0.96, top=0.94, bottom=0.16)
    save_figure(fig, "Fig3_taxonomic_descriptive_contrasts.png")


def rebuild_figure5() -> None:
    meta, rel, _ = load_community()
    community = pd.DataFrame(rel, index=meta["PondMonth"])
    pm_rel = community.groupby(level=0).mean()
    pm_meta = (meta[["PondMonth", "Pond", "Date", "Season_verified"]]
               .drop_duplicates("PondMonth").set_index("PondMonth"))
    keep = pm_meta["Pond"].ne("T1") & pm_meta["Season_verified"].notna()
    pm_rel = pm_rel.loc[pm_meta.index[keep]]
    pm_meta = pm_meta.loc[pm_rel.index]
    dist = squareform(pdist(pm_rel.to_numpy(), metric="braycurtis"))
    scores, _, explained = pcoa(dist)
    coords = pm_meta.reset_index().rename(columns={"Season_verified": "Season"})
    coords["PCoA1"] = scores[:, 0]
    coords["PCoA2"] = scores[:, 1]
    result = pd.read_csv(TABLES / "permanova_dependence_aware_results.csv")
    row = result[(result["Dataset"] == "Grow-out ponds") & (result["Factor"] == "Season")].iloc[0]
    seasonal = pd.read_csv(TABLES / "season_taxa_descriptive.csv")
    mean_columns = [
        "Spring mean (%)", "Summer mean (%)",
        "Autumn mean (%)", "Winter mean (%)",
    ]
    seasonal["Taxon"] = (
        seasonal["Taxon and full taxonomic classification"]
        .str.split(" | ", regex=False)
        .str[0]
    )
    seasonal["Seasonal range"] = (
        seasonal[mean_columns].max(axis=1) - seasonal[mean_columns].min(axis=1)
    )
    seasonal = (
        seasonal.nlargest(20, "Seasonal range")
        .sort_values("Seasonal range")
        .reset_index(drop=True)
    )

    fig = plt.figure(figsize=(12.0, 13.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[0.92, 1.52], hspace=0.24)
    ax = fig.add_subplot(gs[0])
    palette = {"Spring": "#4C78A8", "Summer": "#F58518",
               "Autumn": "#54A24B", "Winter": "#E45756"}
    for season, sub in coords.groupby("Season"):
        ax.scatter(sub["PCoA1"], sub["PCoA2"], s=45, alpha=0.8,
                   color=palette[season], label=season)
    ax.axhline(0, color="#AAAAAA", lw=0.7)
    ax.axvline(0, color="#AAAAAA", lw=0.7)
    ax.set_xlabel(f"PCoA1 ({explained[0]:.1f}%)")
    ax.set_ylabel(f"PCoA2 ({explained[1]:.1f}%)")
    ax.set_title("(a) Seasonal community structure in grow-out ponds\n"
                 "(pond-month means)", loc="left", fontweight="bold")
    ax.text(0.98, 0.03,
            f"Marginal R² = {100 * row['Marginal_R2']:.2f}%\n"
            f"restricted p = {row['p_value']:.4f}; n = {int(row['n'])}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9)
    ax.legend(title="Season", frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")
    remove_top_right_spines(ax)

    subgs = gs[1].subgridspec(1, 2, width_ratios=[1.0, 0.28], wspace=0.08)
    ax2 = fig.add_subplot(subgs[0, 0])
    raw = seasonal[mean_columns].to_numpy(dtype=float)
    means = raw.mean(axis=1, keepdims=True)
    standard_deviations = raw.std(axis=1, ddof=0, keepdims=True)
    standardized = np.divide(
        raw - means,
        standard_deviations,
        out=np.zeros_like(raw),
        where=standard_deviations > 0,
    )
    annotations = np.vectorize(format_abundance)(raw)
    heat = sns.heatmap(
        standardized,
        cmap="RdBu_r",
        center=0,
        vmin=-2,
        vmax=2,
        annot=annotations,
        fmt="",
        annot_kws={"fontsize": 7.6},
        xticklabels=["Spring", "Summer", "Autumn", "Winter"],
        yticklabels=seasonal["Taxon"],
        linewidths=0.3,
        linecolor="white",
        cbar_kws={
            "label": "Row-standardized mean abundance",
            "shrink": 0.72,
            "pad": 0.02,
            "ticks": [-2, -1, 0, 1, 2],
        },
        ax=ax2,
    )
    heat.collections[0].colorbar.ax.tick_params(pad=2)
    ax2.set_xlabel("Season (cell labels: mean relative abundance, %)")
    ax2.set_ylabel("")
    ax2.set_title(
        "(b) Seasonal abundance profiles",
        loc="left",
        fontweight="bold",
        pad=9,
    )

    ax3 = fig.add_subplot(subgs[0, 1], sharey=ax2)
    ax3.barh(
        np.arange(len(seasonal)),
        seasonal["Seasonal range"],
        color="#7186A0",
        alpha=0.9,
        height=0.72,
    )
    for y_value, value in enumerate(seasonal["Seasonal range"]):
        ax3.text(
            value + 0.04, y_value, f"{value:.2f}",
            va="center", fontsize=7.2,
        )
    ax3.set_xlabel("Seasonal range\n(percentage points)")
    ax3.tick_params(axis="y", left=False, labelleft=False)
    ax3.grid(axis="x", alpha=0.16)
    ax3.set_xlim(
        0,
        max(seasonal["Seasonal range"].max() * 1.22, 1),
    )
    remove_top_right_spines(ax3)
    fig.subplots_adjust(left=0.28, right=0.96, top=0.96, bottom=0.055)
    save_figure(fig, "Fig5_season_dependence_aware.png")


def rebuild_figure6() -> None:
    turnover = pd.read_csv(TABLES / "turnover_values.csv")
    turnover["Date_dt"] = pd.to_datetime(
        turnover["Date"].astype(str), format="%Y%m"
    )
    test = pd.read_csv(TABLES / "turnover_dependence_aware_test.csv").iloc[0]
    palette = {"T1": "#8E78B8", "T2": "#E78B93", "T3": "#F2AE3D", "T4": "#557EA6"}
    season_palette = {
        "Spring": "#4C78A8",
        "Summer": "#F58518",
        "Autumn": "#54A24B",
        "Winter": "#E45756",
    }

    fig = plt.figure(figsize=(12.0, 9.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 0.85], hspace=0.28, wspace=0.24)
    ax = fig.add_subplot(gs[0, :])
    for pond, sub in turnover.groupby("Pond"):
        sub = sub.sort_values("Date_dt")
        ax.plot(
            sub["Date_dt"], sub["Turnover"],
            marker="o", ms=3.8, lw=1.1, alpha=0.82,
            color=palette[pond], label=pond,
        )
    ax.set_ylabel("Bray–Curtis turnover")
    ax.set_xlabel("")
    ax.set_title(
        "(a) Turnover between successive available pond-month observations",
        loc="left",
        fontweight="bold",
    )
    ax.legend(title="Pond", frameon=False, ncol=4)
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    ax.grid(axis="y", alpha=0.18)
    remove_top_right_spines(ax)

    ax2 = fig.add_subplot(gs[1, 0])
    sns.boxplot(
        data=turnover, x="Pond", y="Turnover", hue="Pond",
        order=["T1", "T2", "T3", "T4"], palette=palette,
        width=0.55, showfliers=False, legend=False, ax=ax2,
    )
    sns.stripplot(
        data=turnover, x="Pond", y="Turnover",
        order=["T1", "T2", "T3", "T4"],
        color="#3D3D3D", size=3.4, alpha=0.50, ax=ax2,
    )
    ax2.set_title("(b) Turnover by pond", loc="left", fontweight="bold")
    ax2.set_xlabel("")
    ax2.set_ylabel("Bray–Curtis turnover")
    ax2.text(
        0.98, 0.04,
        "Descriptive only\n(no pond-level test)",
        transform=ax2.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
    )
    ax2.grid(axis="y", alpha=0.18)
    remove_top_right_spines(ax2)

    ax3 = fig.add_subplot(gs[1, 1])
    growout = turnover.loc[
        turnover["Pond"].ne("T1") & turnover["Season"].notna()
    ].copy()
    order = ["Spring", "Summer", "Autumn", "Winter"]
    sns.boxplot(
        data=growout, x="Season", y="Turnover", hue="Season",
        order=order, palette=season_palette,
        width=0.58, showfliers=False, legend=False, ax=ax3,
    )
    sns.stripplot(
        data=growout, x="Season", y="Turnover",
        order=order, color="#3D3D3D", size=3.3, alpha=0.50, ax=ax3,
    )
    ax3.set_title(
        "(c) Grow-out turnover by season",
        loc="left",
        fontweight="bold",
    )
    ax3.set_xlabel("")
    ax3.set_ylabel("Bray–Curtis turnover")
    ax3.text(
        0.98, 0.96,
        f"Partial F = {test['F']:.3f}\n"
        f"restricted p = {test['p_value']:.4f}\n"
        f"n = {int(test['n'])} transitions",
        transform=ax3.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78},
    )
    ax3.grid(axis="y", alpha=0.18)
    remove_top_right_spines(ax3)
    fig.tight_layout()
    save_figure(fig, "Fig5_turnover_dependence_aware.png")


def rebuild_figure9() -> None:
    meta, rel, _ = load_community()
    community = pd.DataFrame(rel, index=meta["PondMonth"])
    pm_rel = community.groupby(level=0).mean()
    pm_meta = (
        meta[["PondMonth", "Pond", "Date", "Season_verified"]]
        .drop_duplicates("PondMonth")
        .set_index("PondMonth")
    )
    keep = pm_meta["Pond"].ne("T1") & pm_meta["Season_verified"].notna()
    pm_meta = pm_meta.loc[keep]
    pm_rel = pm_rel.loc[pm_meta.index]
    scores, _, explained = pcoa(
        squareform(pdist(pm_rel.to_numpy(dtype=float), metric="braycurtis"))
    )
    coords = pm_meta.reset_index().rename(columns={"Season_verified": "Season"})
    coords["PCoA1"] = scores[:, 0]
    coords["PCoA2"] = scores[:, 1]
    pm = pd.read_csv(TABLES / "pond_month_analysis_metadata.csv")
    df = coords.merge(pm[["PondMonth", "Typhoon", "ColdSurge",
                          "Heatwave_Tmax_q90_3day", "LowPressure"]],
                      on="PondMonth", how="left")
    results = pd.read_csv(TABLES / "permanova_dependence_aware_results.csv")
    results = results[results["Dataset"].eq("Grow-out ponds")].set_index("Factor")
    explained1 = float(explained[0])
    explained2 = float(explained[1])
    panels = [
        ("Heatwave_Tmax_q90_3day", "Heat event", "#D9534F"),
        ("Typhoon", "Typhoon", "#3978B5"),
        ("ColdSurge", "Cold surge", "#7A4AA5"),
        ("LowPressure", "Low pressure", "#17877C"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 10.2), sharex=True, sharey=True)
    for label, (factor, title, color), ax in zip("abcd", panels, axes.flat):
        event = df[factor].eq(1)
        ax.scatter(df.loc[~event, "PCoA1"], df.loc[~event, "PCoA2"],
                   s=32, color="#D4D4D4", alpha=0.55,
                   label=f"Non-event (n = {(~event).sum()})")
        ax.scatter(df.loc[event, "PCoA1"], df.loc[event, "PCoA2"],
                   s=38, color=color, alpha=0.9,
                   label=f"Event (n = {event.sum()})")
        c0 = df.loc[~event, ["PCoA1", "PCoA2"]].mean()
        c1 = df.loc[event, ["PCoA1", "PCoA2"]].mean()
        ax.annotate("", xy=(c1["PCoA1"], c1["PCoA2"]),
                    xytext=(c0["PCoA1"], c0["PCoA2"]),
                    arrowprops={"arrowstyle": "->", "color": color, "lw": 1.6})
        ax.scatter([c0["PCoA1"], c1["PCoA1"]], [c0["PCoA2"], c1["PCoA2"]],
                   marker="X", s=90, c=["#555555", color],
                   edgecolor="#333333", linewidth=0.6)
        row = results.loc[factor]
        ax.text(0.98, 0.03,
                f"Marginal R² = {100 * row['Marginal_R2']:.2f}%\n"
                f"BH q = {row['q_value']:.3f}",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
                bbox={"facecolor": "white", "edgecolor": "#BBBBBB", "alpha": 0.85})
        ax.axhline(0, color="#AAAAAA", lw=0.7)
        ax.axvline(0, color="#AAAAAA", lw=0.7)
        ax.set_title(f"({label}) {title}", loc="left", fontweight="bold")
        ax.legend(frameon=False, fontsize=8)
        remove_top_right_spines(ax)
    for ax in axes[1]:
        ax.set_xlabel(f"PCoA1 ({explained1:.1f}%)")
    for ax in axes[:, 0]:
        ax.set_ylabel(f"PCoA2 ({explained2:.1f}%)")
    fig.tight_layout()
    save_figure(fig, "Fig8_disturbance_dependence_aware.png")


def rebuild_figure10() -> None:
    per = pd.read_csv(TABLES / "permanova_dependence_aware_results.csv")
    block = pd.read_csv(TABLES / "blockwise_unique_R2_pond_month.csv")
    factors = ["Pond", "Season", "Heatwave_Tmax_q90_3day",
               "Typhoon", "ColdSurge", "LowPressure", "Layer"]
    labels = ["Pond identity", "Season", "Heat event", "Typhoon",
              "Cold surge", "Low pressure", "Sampling layer"]
    datasets = [("All ponds", "#2E86C1"), ("Grow-out ponds", "#F28E2B")]
    values = {}
    annotations = {}
    for dataset, _ in datasets:
        for factor in factors:
            hit = per[(per["Dataset"] == dataset) & (per["Factor"] == factor)]
            if hit.empty:
                values[(dataset, factor)] = np.nan
                continue
            row = hit.iloc[0]
            values[(dataset, factor)] = 100 * row["Marginal_R2"]
            if factor == "Pond":
                annotations[(dataset, factor)] = "descriptive"
            elif np.isfinite(row["q_value"]):
                annotations[(dataset, factor)] = f"q={row['q_value']:.3f}"
            elif np.isfinite(row["p_value"]):
                annotations[(dataset, factor)] = f"p={row['p_value']:.4f}"

    fig, axes = plt.subplots(2, 1, figsize=(11.5, 12.8))
    y = np.arange(len(factors))
    h = 0.35
    for i, (dataset, color) in enumerate(datasets):
        vals = [values[(dataset, f)] for f in factors]
        offset = (i - 0.5) * h
        axes[0].barh(y + offset, np.nan_to_num(vals, nan=0), height=h,
                     color=color, label=dataset)
        for yy, factor, value in zip(y + offset, factors, vals):
            if np.isfinite(value):
                note = annotations.get((dataset, factor), "")
                axes[0].text(value + 0.12, yy, f"{value:.2f}% ({note})",
                             va="center", fontsize=8.3)
    axes[0].set_yticks(y, labels)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Marginal community variation explained, R² (%)")
    axes[0].set_title("(a) Dependence-aware community-structure analyses",
                      loc="left", fontweight="bold")
    axes[0].legend(frameon=False, loc="lower right")
    axes[0].grid(axis="x", alpha=0.2)
    finite_values = [v for v in values.values() if np.isfinite(v)]
    axes[0].set_xlim(0, max(finite_values) + 3.6)
    remove_top_right_spines(axes[0])

    block_order = ["Pond identity", "Season", "Water quality",
                   "Meteorological disturbances"]
    y2 = np.arange(len(block_order))
    for i, (dataset, color) in enumerate(datasets):
        sub = block[block["Dataset"] == dataset].set_index("Predictor_block")
        vals = [100 * sub.loc[x, "Unique_R2"] for x in block_order]
        offset = (i - 0.5) * h
        axes[1].barh(y2 + offset, vals, height=h, color=color, label=dataset)
        for yy, value in zip(y2 + offset, vals):
            axes[1].text(value + 0.10, yy, f"{value:.2f}%", va="center", fontsize=8.5)
    axes[1].set_yticks(y2, block_order)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("Unique community variation explained (%)")
    axes[1].set_title("(b) Unique contributions in the multivariable model",
                      loc="left", fontweight="bold")
    axes[1].grid(axis="x", alpha=0.2)
    max_block = 100 * block["Unique_R2"].max()
    axes[1].set_xlim(0, max_block + 1.8)
    remove_top_right_spines(axes[1])
    fig.tight_layout()
    save_figure(fig, "Fig9_dependence_aware_variance_partition.png")


def rebuild_environmental_vectors() -> None:
    coords = pd.read_csv(TABLES / "pond_month_pcoa_coordinates_T2_T4.csv")
    vectors = pd.read_csv(TABLES / "environmental_vector_fits_dependence_aware.csv")
    quality = pd.read_json(RESULTS / "analysis_quality_report.json", typ="series")
    explained1 = float(quality["pcoa1_explained"])
    explained2 = float(quality["pcoa2_explained"])

    fig, ax = plt.subplots(figsize=(8.2, 6.8))
    palette = {
        "Spring": "#6BAF92",
        "Summer": "#F2A65A",
        "Autumn": "#C67C5C",
        "Winter": "#6C8EBF",
    }
    for season, sub in coords.groupby("Season"):
        ax.scatter(
            sub["PCoA1"], sub["PCoA2"], s=45, alpha=0.8,
            label=season, color=palette[season],
        )
    boundary = coords[coords["Season"].isna()]
    if len(boundary):
        ax.scatter(
            boundary["PCoA1"], boundary["PCoA2"], s=48, alpha=0.85,
            label="Boundary-spanning average", color="#777777", marker="X",
        )

    scale = 0.8 * min(coords["PCoA1"].max() - coords["PCoA1"].min(),
                      coords["PCoA2"].max() - coords["PCoA2"].min())
    offsets = {
        "Temperature": (0.012, 0.004, "left"),
        "Salinity": (0.008, 0.008, "left"),
        "pH": (0.020, 0.010, "left"),
        "DO": (0.008, -0.011, "left"),
        "Turbidity": (-0.004, 0.012, "right"),
        "MeanPressure": (-0.004, -0.026, "center"),
        "MinPressure": (0.010, -0.014, "left"),
        "Tmax": (-0.010, 0.012, "right"),
        "TotalRainfall": (-0.006, 0.021, "right"),
        "MaxDailyRainfall": (0.008, 0.017, "left"),
    }
    display_labels = {
        "Temperature": "Water temperature",
        "Salinity": "Salinity",
        "DO": "Dissolved oxygen",
        "Tmax": "Maximum air temperature",
        "TotalRainfall": "Total rainfall",
        "MaxDailyRainfall": "Maximum daily rainfall",
        "MeanPressure": "Mean pressure",
        "MinPressure": "Minimum pressure",
    }
    for _, row in vectors.iterrows():
        x = row["Correlation_PCoA1"] * scale
        y = row["Correlation_PCoA2"] * scale
        significant = row["q_value"] < 0.05
        color = "#A61B1B" if significant else "#3F3F3F"
        ax.arrow(
            0, 0, x, y, color=color, alpha=0.82, width=0.00035,
            head_width=0.008, length_includes_head=True,
        )
        dx, dy, ha = offsets[row["Environmental_variable"]]
        display_label = display_labels.get(
            row["Environmental_variable"], row["Environmental_variable"]
        )
        ax.text(
            x + dx, y + dy, display_label,
            fontsize=8.2, ha=ha, va="center", color=color,
            fontweight="bold" if significant else "normal",
        )

    ax.axhline(0, color="#AAAAAA", lw=0.8)
    ax.axvline(0, color="#AAAAAA", lw=0.8)
    ax.set_xlim(left=min(coords["PCoA1"].min() - 0.035, -0.335))
    ax.set_xlabel(f"PCoA1 ({explained1:.1f}%)")
    ax.set_ylabel(f"PCoA2 ({explained2:.1f}%)")
    ax.set_title(
        "Environmental alignment in grow-out ponds (pond-month means)",
        fontweight="bold",
    )
    ax.legend(
        title="Season", frameon=False,
        bbox_to_anchor=(1.02, 1), loc="upper left",
    )
    remove_top_right_spines(ax)
    fig.tight_layout()
    save_figure(fig, "Fig7_environmental_vectors_dependence_aware.png")


if __name__ == "__main__":
    rebuild_figure2()
    rebuild_figure3()
    rebuild_figure4()
    rebuild_figure5()
    rebuild_figure6()
    rebuild_figure9()
    rebuild_figure10()
    rebuild_environmental_vectors()
