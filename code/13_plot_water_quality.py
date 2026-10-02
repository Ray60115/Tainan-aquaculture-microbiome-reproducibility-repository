#!/usr/bin/env python3
"""Rebuild the pond-month water-quality time-series figure (current Fig. 6)."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata)
    metadata["Date_plot"] = pd.to_datetime(
        metadata["Date"].astype(str), format="%Y%m"
    )
    variables = [
        ("Temperature", "Water temperature (°C)"),
        ("Salinity", "Salinity (specific gravity)"),
        ("pH", "pH"),
        ("DO", "Dissolved oxygen"),
        ("Turbidity", "Turbidity (NTU)"),
    ]
    pond_month = (
        metadata.groupby(["Pond", "Date_plot"], as_index=False)[
            [name for name, _ in variables]
        ]
        .mean()
        .sort_values(["Pond", "Date_plot"])
    )
    colors = {
        "T1": "#4C78A8",
        "T2": "#F58518",
        "T3": "#54A24B",
        "T4": "#E45756",
    }

    fig, axes = plt.subplots(5, 1, figsize=(10.5, 12.0), sharex=True)
    for panel, (axis, (variable, label)) in enumerate(
        zip(axes, variables, strict=True)
    ):
        for pond, group in pond_month.groupby("Pond", sort=True):
            axis.plot(
                group["Date_plot"],
                group[variable],
                marker="o",
                markersize=3.5,
                linewidth=1.2,
                color=colors.get(pond),
                label=pond,
            )
        axis.set_ylabel(label)
        axis.set_title(f"({'abcde'[panel]})", loc="left", fontweight="bold")
        axis.grid(axis="y", alpha=0.18)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
    axes[0].legend(title="Pond", frameon=False, ncol=4, loc="upper right")
    axes[-1].xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    axes[-1].tick_params(axis="x", rotation=45)
    axes[-1].set_xlabel("Sampling month")
    fig.tight_layout()

    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=400, bbox_inches="tight", transparent=True)
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight", transparent=True)
    plt.close(fig)


if __name__ == "__main__":
    main()
