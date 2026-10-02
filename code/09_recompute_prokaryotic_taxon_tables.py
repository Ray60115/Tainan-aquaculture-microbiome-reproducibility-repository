#!/usr/bin/env python3
"""Recompute selected-taxon descriptive tables after organelle removal.

The selected display lineages are supplied as a machine-readable CSV rather
than being extracted from a Word supplement.  This keeps the public workflow
self-contained and makes the selection set explicit.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


ORGANELLE_PATTERN = re.compile(
    r"mitochond|chloroplast|d__eukary", flags=re.IGNORECASE
)
RANKS = ["phylum", "class", "order", "family", "genus", "species"]
PREFIX = {
    "phylum": "p__",
    "class": "c__",
    "order": "o__",
    "family": "f__",
    "genus": "g__",
    "species": "s__",
}


def normalize_taxon_name(value: str) -> str:
    return re.sub(r"\s+", " ", str(value).replace("_", " ")).strip()


def split_asv_taxonomy(value: str) -> dict[str, str]:
    result = {}
    for item in str(value).split(";"):
        item = item.strip()
        for rank, prefix in PREFIX.items():
            if item.startswith(prefix):
                result[rank] = normalize_taxon_name(item[len(prefix) :])
    return result


def target_from_unlabelled(classification: str) -> dict[str, str]:
    values = [value.strip() for value in classification.split(";") if value.strip()]
    if len(values) > len(RANKS):
        raise ValueError(f"Unexpected classification depth: {classification}")
    return {
        rank: normalize_taxon_name(value)
        for rank, value in zip(RANKS, values, strict=False)
    }


def target_from_labelled(description: str) -> tuple[str, dict[str, str]]:
    representative, annotation = description.split("|", 1)
    target = {}
    for part in annotation.split(";"):
        if ":" not in part:
            continue
        rank, value = part.split(":", 1)
        rank = rank.strip().lower()
        if rank in RANKS:
            target[rank] = normalize_taxon_name(value)
    return representative.strip(), target


def match_target(
    parsed_taxonomy: pd.Series, target: dict[str, str]
) -> np.ndarray:
    mask = np.ones(len(parsed_taxonomy), dtype=bool)
    for rank, value in target.items():
        mask &= parsed_taxonomy.map(lambda item: item.get(rank, "") == value).to_numpy()
    return mask


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--asv", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--display-definitions", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    asv = pd.read_csv(args.asv)
    meta = pd.read_csv(args.metadata)
    samples = meta["SampleID"].tolist()
    if len(samples) != 168:
        raise ValueError("Expected 168 formal samples.")

    taxonomy = asv["Taxon"].fillna("").astype(str)
    positive = asv[samples].fillna(0).sum(axis=1) > 0
    organelle = taxonomy.str.contains(ORGANELLE_PATTERN)
    kept = asv.loc[positive & ~organelle].copy().reset_index(drop=True)
    counts = kept[samples].T.astype(float)
    counts.columns = kept.index
    relative = counts.div(counts.sum(axis=1), axis=0)
    parsed = kept["Taxon"].map(split_asv_taxonomy)

    meta = meta.set_index("SampleID").loc[samples].reset_index()
    if "Season_verified" not in meta.columns:
        meta["Season_verified"] = meta["Season"]
    if int(meta["Season_verified"].notna().sum()) != 160:
        raise ValueError("Expected 160 unambiguous sample-level season labels.")
    meta_indexed = meta.set_index("SampleID")

    definitions = pd.read_csv(args.display_definitions)
    required = {"Set", "Representative_taxon", "Target_definition", "Target_format"}
    if not required.issubset(definitions.columns):
        raise ValueError(
            "Display definitions lack required columns: "
            f"{sorted(required - set(definitions.columns))}"
        )
    source_stage = definitions.loc[definitions["Set"].eq("stage")]
    source_season = definitions.loc[definitions["Set"].eq("season")]
    if len(source_stage) != 20 or len(source_season) != 30:
        raise ValueError(
            "Expected 20 stage lineages and 30 seasonal lineages; found "
            f"{len(source_stage)} and {len(source_season)}."
        )

    stage_rows = []
    pond_rows = []
    for _, row in source_stage.iterrows():
        representative = row["Representative_taxon"]
        classification = row["Target_definition"]
        if row["Target_format"] != "unlabelled_full_path":
            raise ValueError("Unexpected target format in stage definitions.")
        target = target_from_unlabelled(classification)
        selected = match_target(parsed, target)
        if not selected.any():
            raise ValueError(f"No ASV matched stage taxon: {representative}")
        abundance = relative.loc[:, selected].sum(axis=1)
        present = counts.loc[:, selected].sum(axis=1) > 0
        nursery_ids = meta.loc[meta["Pond"].eq("T1"), "SampleID"]
        grow_ids = meta.loc[meta["Pond"].ne("T1"), "SampleID"]
        nursery_mean = float(abundance.loc[nursery_ids].mean() * 100)
        grow_mean = float(abundance.loc[grow_ids].mean() * 100)
        ratio = (
            float(np.log2(nursery_mean / grow_mean))
            if nursery_mean > 0 and grow_mean > 0
            else np.nan
        )
        stage_rows.append(
            {
                "Representative taxon": representative,
                "Full taxonomic classification": classification,
                "Mean relative abundance, nursery (%)": nursery_mean,
                "Mean relative abundance, grow-out (%)": grow_mean,
                "Mean difference (nursery − grow-out; percentage points)": nursery_mean
                - grow_mean,
                "log₂ fold change (nursery vs grow-out)": ratio,
                "Prevalence, nursery (proportion)": float(
                    present.loc[nursery_ids].mean()
                ),
                "Prevalence, grow-out (proportion)": float(
                    present.loc[grow_ids].mean()
                ),
            }
        )
        pond_row = {"Representative taxon": representative}
        for pond in ["T1", "T2", "T3", "T4"]:
            ids = meta.loc[meta["Pond"].eq(pond), "SampleID"]
            pond_row[f"{pond} (%)"] = float(abundance.loc[ids].mean() * 100)
        pond_rows.append(pond_row)

    seasonal_rows = []
    grow_meta = meta[meta["Pond"].ne("T1")]
    for description in source_season["Target_definition"]:
        representative, target = target_from_labelled(description)
        selected = match_target(parsed, target)
        if not selected.any():
            raise ValueError(f"No ASV matched seasonal taxon: {representative}")
        abundance = relative.loc[:, selected].sum(axis=1)
        present = counts.loc[:, selected].sum(axis=1) > 0
        output = {"Taxon and full taxonomic classification": description}
        for season in ["Spring", "Summer", "Autumn", "Winter"]:
            ids = grow_meta.loc[
                grow_meta["Season_verified"].eq(season), "SampleID"
            ]
            output[f"{season} mean (%)"] = float(abundance.loc[ids].mean() * 100)
            output[f"{season} prevalence"] = float(present.loc[ids].mean())
        seasonal_rows.append(output)

    pd.DataFrame(stage_rows).to_csv(
        args.outdir / "stage_taxa_descriptive.csv", index=False
    )
    pd.DataFrame(pond_rows).to_csv(
        args.outdir / "stage_taxa_by_pond.csv", index=False
    )
    pd.DataFrame(seasonal_rows).to_csv(
        args.outdir / "season_taxa_descriptive.csv", index=False
    )


if __name__ == "__main__":
    main()
