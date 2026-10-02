#!/usr/bin/env python3
"""Apply the user-supplied sampling-occasion season labels to the 168-sample metadata.

The source file distinguishes first and second sampling occasions (``-1`` and
``-2``).  The ASV table used for the main analysis already averages those
occasions within pond, month, and layer.  If the contributing occasions cross
a solar-term season boundary, the averaged microbiome profile cannot be given
one unambiguous season.  Such profiles are retained for non-seasonal analyses
but assigned a missing ``Season`` value for every model that requires season.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


SAMPLE_RE = re.compile(r"^(T[1-4])-(\d{6})(?:-(\d+))?-([UD])$")
VALID_SEASONS = {"Spring", "Summer", "Autumn", "Winter"}


def canonical_sample(value: str) -> tuple[str, str, int, str, str | None]:
    match = SAMPLE_RE.fullmatch(str(value).strip())
    if match is None:
        raise ValueError(f"Unrecognized sample label in season source: {value!r}")
    pond, date, occasion, layer = match.groups()
    return f"{pond}-{date}-{layer}", pond, int(date), layer, occasion


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--seasons", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--excluded", required=True, type=Path)
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata)
    source = pd.read_csv(args.seasons, encoding="utf-8-sig")
    if len(metadata) != 168 or metadata["SampleID"].duplicated().any():
        raise ValueError("Expected 168 unique formal metadata rows.")
    if set(source.columns) != {"Sample", "Group", "condition"}:
        raise ValueError("Season source must contain Sample, Group, and condition.")
    if not set(source["condition"].dropna()).issubset(VALID_SEASONS):
        raise ValueError("Unexpected season value in authoritative source.")

    parsed = source["Sample"].map(canonical_sample)
    source = source.copy()
    source[["CanonicalSample", "Pond", "Date", "Layer", "Occasion"]] = pd.DataFrame(
        parsed.tolist(), index=source.index
    )

    grouped_rows: list[dict] = []
    for sample_id, frame in source.groupby("CanonicalSample", sort=False):
        seasons = frame["condition"].dropna().drop_duplicates().tolist()
        grouped_rows.append(
            {
                "SampleID": sample_id,
                "Source_records": "; ".join(frame["Sample"].astype(str)),
                "Source_occasions": "; ".join(
                    sorted({str(value) for value in frame["Occasion"].dropna()})
                ),
                "Source_seasons": "; ".join(seasons),
                "Season_attachment": seasons[0] if len(seasons) == 1 else pd.NA,
                "Season_status": (
                    "Direct unambiguous assignment"
                    if len(seasons) == 1
                    else "Boundary-spanning averaged profile"
                ),
            }
        )
    grouped = pd.DataFrame(grouped_rows)

    metadata = metadata.copy()
    metadata["Season_original_metadata"] = metadata["Season"]
    metadata = metadata.merge(grouped, on="SampleID", how="left", validate="one_to_one")

    # The attachment contains only a surface entry for three December 2024
    # grow-out pond-months.  The corresponding bottom profiles inherit the
    # same unambiguous season from their same-pond, same-date surface record.
    source_by_pond_date = (
        source.groupby(["Pond", "Date"])["condition"]
        .agg(lambda values: list(pd.unique(values.dropna())))
    )
    source_by_date = source.groupby("Date")["condition"].agg(
        lambda values: list(pd.unique(values.dropna()))
    )
    for index in metadata.index[metadata["Season_status"].isna()]:
        row = metadata.loc[index]
        seasons = source_by_pond_date.get((row["Pond"], int(row["Date"])), [])
        status = "Inherited from same pond and sampling date"
        source_note = "same pond/date record in attachment"
        if len(seasons) != 1:
            seasons = source_by_date.get(int(row["Date"]), [])
            status = "Inherited from unambiguous records in the same sampling month"
            source_note = "same sampling-month records in attachment"
        if len(seasons) != 1:
            raise ValueError(
                f"No unambiguous attachment-based season for {row['SampleID']}."
            )
        metadata.at[index, "Season_attachment"] = seasons[0]
        metadata.at[index, "Season_status"] = status
        metadata.at[index, "Source_seasons"] = seasons[0]
        metadata.at[index, "Source_records"] = source_note

    metadata["Season"] = metadata["Season_attachment"]
    metadata["Season_verified"] = metadata["Season_attachment"]
    metadata["Season_unambiguous"] = metadata["Season_verified"].notna()

    boundary_samples = metadata.loc[
        metadata["Season_status"].eq("Boundary-spanning averaged profile"),
        "SampleID",
    ].tolist()
    expected_boundary = {
        "T1-202302-D",
        "T1-202302-U",
        "T2-202408-D",
        "T2-202408-U",
        "T3-202408-D",
        "T3-202408-U",
        "T4-202408-D",
        "T4-202408-U",
    }
    if set(boundary_samples) != expected_boundary:
        raise ValueError(
            "Unexpected boundary-spanning set: " + ", ".join(boundary_samples)
        )
    if int(metadata["Season_unambiguous"].sum()) != 160:
        raise ValueError("Expected 160 samples with unambiguous season assignment.")

    audit_columns = [
        "SampleID",
        "Pond",
        "Date",
        "Layer",
        "Season_original_metadata",
        "Season_attachment",
        "Season_status",
        "Source_records",
        "Source_occasions",
        "Source_seasons",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(args.output, index=False)
    metadata[audit_columns].to_csv(args.audit, index=False)

    excluded = (
        metadata.loc[~metadata["Season_unambiguous"], ["Pond", "Date"]]
        .drop_duplicates()
        .sort_values(["Pond", "Date"])
    )
    excluded["PondMonth"] = excluded["Pond"] + "-" + excluded["Date"].astype(str)
    excluded["Reason"] = (
        "Averaged first and second sampling occasions fall in different "
        "solar-term seasons in the authoritative attachment"
    )
    excluded[["PondMonth", "Pond", "Date", "Reason"]].to_csv(
        args.excluded, index=False
    )
    print(
        f"Wrote {len(metadata)} rows; {metadata['Season_unambiguous'].sum()} "
        f"have unambiguous seasons and {len(excluded)} pond-months are excluded "
        "from season-dependent inference."
    )


if __name__ == "__main__":
    main()
