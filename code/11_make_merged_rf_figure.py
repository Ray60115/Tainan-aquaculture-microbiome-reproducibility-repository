#!/usr/bin/env python3
"""Merge random-forest reconstruction and importance into final Figure 11."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-tainan-rf")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


POND_COLORS = {
    "T1": "#6A51A3",
    "T2": "#D95F5F",
    "T3": "#E69F00",
    "T4": "#2C7FB8",
}


def clean_axis(axis: plt.Axes) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.grid(color="#E6E6E6", linewidth=0.6, alpha=0.75)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    tables = args.results / "tables"
    performance = pd.read_csv(tables / "rf_cross_validated_performance.csv")
    predictions = pd.read_csv(tables / "rf_oof_predictions.csv")
    importance = pd.read_csv(tables / "rf_grouped_permutation_importance.csv")

    sns.set_theme(style="white", context="paper", font_scale=1.02)
    # Panels (b) and (c) are deliberately placed side by side so that the two
    # reconstructed PCoA axes can be compared directly. Full-width panels (a)
    # and (d) preserve legibility of model labels and importance estimates.
    fig = plt.figure(figsize=(13.0, 14.0), constrained_layout=True)
    grid = fig.add_gridspec(
        3,
        2,
        height_ratios=[0.72, 1.0, 1.18],
    )
    ax_a = fig.add_subplot(grid[0, :])
    ax_b = fig.add_subplot(grid[1, 0])
    ax_c = fig.add_subplot(grid[1, 1])
    ax_d = fig.add_subplot(grid[2, :])

    positions = np.arange(len(performance))
    width = 0.34
    bars1 = ax_a.bar(
        positions - width / 2,
        performance["PCoA1_R2"],
        width,
        color="#4C78A8",
        label="PCoA1",
    )
    bars2 = ax_a.bar(
        positions + width / 2,
        performance["PCoA2_R2"],
        width,
        color="#F58518",
        label="PCoA2",
    )
    ax_a.bar_label(bars1, fmt="%.3f", padding=2, fontsize=8)
    ax_a.bar_label(bars2, fmt="%.3f", padding=2, fontsize=8)
    ax_a.set_xticks(
        positions,
        ["Structure", "+ Season", "+ Water\nquality", "+ Disturbances"],
    )
    ax_a.set_ylabel("Grouped cross-validated R²")
    ax_a.set_title(
        "(a) Sequential-model performance", loc="left", fontweight="bold"
    )
    ax_a.set_ylim(
        min(-0.03, float(performance[["PCoA1_R2", "PCoA2_R2"]].min().min()) - 0.03),
        float(performance[["PCoA1_R2", "PCoA2_R2"]].max().max()) + 0.11,
    )
    ax_a.axhline(0, color="#777777", linewidth=0.7)
    ax_a.legend(frameon=False, ncol=2, loc="upper left")
    clean_axis(ax_a)

    plot = importance.sort_values("Mean_increase_MAE_percent")
    ax_d.barh(
        plot["Predictor_group"],
        plot["Mean_increase_MAE_percent"],
        xerr=plot["SD_across_repeats_percent"],
        color="#4C78A8",
        alpha=0.9,
        capsize=2,
    )
    span = float(
        (
            plot["Mean_increase_MAE_percent"]
            + plot["SD_across_repeats_percent"]
        ).max()
    )
    for y_position, (value, error) in enumerate(
        zip(
            plot["Mean_increase_MAE_percent"],
            plot["SD_across_repeats_percent"],
        )
    ):
        ax_d.text(
            float(value + error + max(0.35, span * 0.012)),
            y_position,
            f"{value:.2f}".replace("-", "−"),
            va="center",
            fontsize=7.8,
        )
    ax_d.axvline(0, color="#777777", linewidth=0.7)
    ax_d.set_xlabel("Increase in held-out MAE after permutation (%)")
    ax_d.set_title(
        "(d) Grouped permutation importance", loc="left", fontweight="bold"
    )
    ax_d.set_xlim(
        min(
            -1.0,
            float(
                (
                    plot["Mean_increase_MAE_percent"]
                    - plot["SD_across_repeats_percent"]
                ).min()
            )
            - 1,
        ),
        span + 6,
    )
    clean_axis(ax_d)

    full = performance.iloc[-1]
    for axis_number, axis in [(1, ax_b), (2, ax_c)]:
        observed = predictions[f"PCoA{axis_number}"]
        reconstructed = predictions[f"Predicted_PCoA{axis_number}"]
        for pond in ["T1", "T2", "T3", "T4"]:
            selected = predictions["Pond"].eq(pond)
            axis.scatter(
                observed[selected],
                reconstructed[selected],
                s=29,
                color=POND_COLORS[pond],
                alpha=0.78,
                edgecolor="white",
                linewidth=0.35,
                label=pond,
            )
        lower = float(min(observed.min(), reconstructed.min()))
        upper = float(max(observed.max(), reconstructed.max()))
        axis.plot(
            [lower, upper],
            [lower, upper],
            linestyle="--",
            color="#777777",
            linewidth=1,
        )
        r2 = float(full[f"PCoA{axis_number}_R2"])
        mae = float(full[f"PCoA{axis_number}_MAE"])
        axis.text(
            0.04,
            0.95,
            f"R² = {r2:.3f}\nMAE = {mae:.3f}",
            transform=axis.transAxes,
            va="top",
            ha="left",
            fontsize=9.4,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 2},
        )
        axis.set_xlabel(f"Observed global PCoA{axis_number} score")
        axis.set_ylabel(f"Out-of-fold reconstructed PCoA{axis_number} score")
        axis.set_title(
            f"({'b' if axis_number == 1 else 'c'}) PCoA{axis_number}",
            loc="left",
            fontweight="bold",
        )
        clean_axis(axis)
    ax_c.legend(title="Pond", frameon=False, ncol=2, loc="lower right")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=400, bbox_inches="tight")
    fig.savefig(args.output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
