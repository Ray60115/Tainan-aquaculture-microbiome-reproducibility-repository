#!/usr/bin/env python3
"""Append November–December 2024 C0X120 observations and update metadata.

The existing daily heat-day definitions are retained exactly. In particular,
Tmax_gt_q90 denotes a daily maximum temperature of at least 34.8 °C, matching
the prespecified threshold used throughout the manuscript.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


DAILY_COLUMNS = [
    "FullDate",
    "AirTemp",
    "Tmax",
    "Tmin",
    "Rainfall",
    "WindSpeed",
    "MaxGust",
    "Pressure",
    "PressureMin",
    "PressureMax",
    "RH",
    "Date",
    "SourceFile",
    "Tmax_gt_q85",
    "Tmax_gt_q90",
    "Tmax_gt_q95",
]

RAW_TO_CLEAN = {
    "Temperature": "AirTemp",
    "T Max": "Tmax",
    "T Min": "Tmin",
    "Precp": "Rainfall",
    "WS": "WindSpeed",
    "WSGust": "MaxGust",
    "StnPres": "Pressure",
    "StnPresMin": "PressureMin",
    "StnPresMax": "PressureMax",
    "RH": "RH",
}


def read_raw_month(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, skiprows=[0])
    year_month = path.stem.split("-")[-2:]
    year = int(year_month[0])
    month = int(year_month[1])
    expected_days = pd.Period(f"{year}-{month:02d}").days_in_month
    if len(raw) != expected_days:
        raise ValueError(
            f"{path.name} contains {len(raw)} rows; expected {expected_days}."
        )
    days = pd.to_numeric(raw["ObsTime"], errors="raise").astype(int)
    if days.tolist() != list(range(1, expected_days + 1)):
        raise ValueError(f"{path.name} does not contain a complete daily sequence.")

    clean = pd.DataFrame()
    clean["FullDate"] = pd.to_datetime(
        {"year": year, "month": month, "day": days}
    ).dt.strftime("%Y-%m-%d")
    for raw_name, clean_name in RAW_TO_CLEAN.items():
        clean[clean_name] = pd.to_numeric(raw[raw_name], errors="coerce")
        if clean[clean_name].isna().any():
            missing_days = days.loc[clean[clean_name].isna()].tolist()
            raise ValueError(
                f"{path.name}: nonnumeric or missing {raw_name} on days {missing_days}."
            )
    clean["Date"] = year * 100 + month
    clean["SourceFile"] = path.name
    return clean


def heatwave_months(daily: pd.DataFrame) -> set[int]:
    """Return months touched by a run of at least three q90 heat days."""
    ordered = daily.sort_values("FullDate").reset_index(drop=True)
    dates = pd.to_datetime(ordered["FullDate"])
    flags = ordered["Tmax_gt_q90"].astype(int).to_numpy()
    event_months: set[int] = set()
    start = 0
    while start < len(ordered):
        if flags[start] == 0:
            start += 1
            continue
        end = start + 1
        while (
            end < len(ordered)
            and flags[end] == 1
            and (dates.iloc[end] - dates.iloc[end - 1]).days == 1
        ):
            end += 1
        if end - start >= 3:
            event_months.update(
                (dates.iloc[start:end].dt.year * 100 + dates.iloc[start:end].dt.month)
                .astype(int)
                .tolist()
            )
        start = end
    return event_months


def monthly_summary(daily: pd.DataFrame) -> pd.DataFrame:
    heatwave = heatwave_months(daily)
    summary = daily.groupby("Date", as_index=False).agg(
        HeatEvent_Tmax_q90_any=("Tmax_gt_q90", "max"),
        Days_Tmax_gt_q90=("Tmax_gt_q90", "sum"),
        MeanAirTemp=("AirTemp", "mean"),
        Tmax=("Tmax", "max"),
        Tmin=("Tmin", "min"),
        TotalRainfall=("Rainfall", "sum"),
        MaxDailyRainfall=("Rainfall", "max"),
        RainDays=("Rainfall", lambda x: int((x > 0).sum())),
        MeanWindSpeed=("WindSpeed", "mean"),
        MaxGust=("MaxGust", "max"),
        MeanPressure=("Pressure", "mean"),
        MinPressure=("PressureMin", "min"),
    )
    summary["Heatwave_Tmax_q90_3day"] = (
        summary["Date"].astype(int).isin(heatwave).astype(int)
    )
    return summary


def validate_existing_monthly_values(
    metadata: pd.DataFrame, summary: pd.DataFrame
) -> None:
    columns = [
        "HeatEvent_Tmax_q90_any",
        "Heatwave_Tmax_q90_3day",
        "Days_Tmax_gt_q90",
        "MeanAirTemp",
        "Tmax",
        "Tmin",
        "TotalRainfall",
        "MaxDailyRainfall",
        "RainDays",
        "MeanWindSpeed",
        "MaxGust",
        "MeanPressure",
        "MinPressure",
    ]
    observed = (
        metadata[["Date", *columns]]
        .drop_duplicates("Date")
        .dropna(subset=["MeanAirTemp"])
        .sort_values("Date")
    )
    check = observed.merge(
        summary[["Date", *columns]],
        on="Date",
        suffixes=("_metadata", "_derived"),
        validate="one_to_one",
    )
    for column in columns:
        left = pd.to_numeric(check[f"{column}_metadata"], errors="coerce")
        right = pd.to_numeric(check[f"{column}_derived"], errors="coerce")
        if not np.allclose(left, right, equal_nan=True):
            mismatch = check.loc[
                ~np.isclose(left, right, equal_nan=True),
                ["Date", f"{column}_metadata", f"{column}_derived"],
            ]
            raise ValueError(
                f"Derived {column} does not reproduce existing metadata:\n"
                f"{mismatch.to_string(index=False)}"
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--existing-daily", required=True, type=Path)
    parser.add_argument("--existing-metadata", required=True, type=Path)
    parser.add_argument("--new-raw", nargs=2, required=True, type=Path)
    parser.add_argument("--complete-daily", required=True, type=Path)
    parser.add_argument("--complete-metadata", required=True, type=Path)
    args = parser.parse_args()

    existing = pd.read_csv(args.existing_daily)
    existing["FullDate"] = pd.to_datetime(existing["FullDate"]).dt.strftime("%Y-%m-%d")
    new = pd.concat(
        [read_raw_month(path) for path in args.new_raw],
        ignore_index=True,
    )

    # Preserve the already documented threshold definitions rather than
    # recalculating them after appending two cool-season months.
    thresholds = {}
    for flag in ["Tmax_gt_q85", "Tmax_gt_q90", "Tmax_gt_q95"]:
        thresholds[flag] = float(existing.loc[existing[flag].eq(0), "Tmax"].max())
        new[flag] = (new["Tmax"] > thresholds[flag]).astype(int)
    if thresholds["Tmax_gt_q90"] != 34.7:
        raise ValueError(
            "The existing q90 heat-day rule is not equivalent to Tmax ≥ 34.8 °C."
        )

    complete = pd.concat([existing[DAILY_COLUMNS], new[DAILY_COLUMNS]], ignore_index=True)
    complete = complete.sort_values("FullDate").reset_index(drop=True)
    dates = pd.to_datetime(complete["FullDate"])
    if complete["FullDate"].duplicated().any():
        raise ValueError("Duplicate dates after appending the new observations.")
    expected = pd.date_range("2022-07-01", "2024-12-31", freq="D")
    if not dates.equals(pd.Series(expected)):
        missing = expected.difference(pd.DatetimeIndex(dates))
        extra = pd.DatetimeIndex(dates).difference(expected)
        raise ValueError(f"Daily coverage is not continuous; missing={missing}, extra={extra}.")
    args.complete_daily.parent.mkdir(parents=True, exist_ok=True)
    complete.to_csv(args.complete_daily, index=False)

    summary = monthly_summary(complete)
    metadata = pd.read_csv(args.existing_metadata)
    validate_existing_monthly_values(metadata, summary)

    update_columns = [
        "HeatEvent_Tmax_q90_any",
        "Heatwave_Tmax_q90_3day",
        "Days_Tmax_gt_q90",
        "MeanAirTemp",
        "Tmax",
        "Tmin",
        "TotalRainfall",
        "MaxDailyRainfall",
        "RainDays",
        "MeanWindSpeed",
        "MaxGust",
        "MeanPressure",
        "MinPressure",
    ]
    summary_by_date = summary.set_index("Date")
    target = metadata["Date"].isin([202411, 202412])
    for column in update_columns:
        metadata.loc[target, column] = (
            metadata.loc[target, "Date"].map(summary_by_date[column]).to_numpy()
        )
    if metadata.loc[target, update_columns].isna().any().any():
        raise ValueError("Updated November–December metadata still contain missing climate values.")
    args.complete_metadata.parent.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(args.complete_metadata, index=False)

    report = summary.loc[
        summary["Date"].isin([202411, 202412]),
        ["Date", *update_columns],
    ]
    print(f"Saved {len(complete)} daily observations to {args.complete_daily.name}.")
    print(f"Saved {len(metadata)} metadata rows to {args.complete_metadata.name}.")
    print(report.to_string(index=False))


if __name__ == "__main__":
    main()
