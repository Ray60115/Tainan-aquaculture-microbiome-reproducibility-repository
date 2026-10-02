#!/usr/bin/env python3
"""Repeat dependence-aware taxonomic analyses after organelle removal.

This wrapper deliberately reuses the already-audited dependence-aware
statistical functions and permutation restrictions in
``reanalyze_dependence_aware.py``. It changes only the feature inclusion rule:
mitochondria, chloroplast, and explicit eukaryotic records are removed before
sample-wise normalization, Shannon calculation, Bray-Curtis dissimilarities,
ordination, turnover, and environmental analyses.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ORGANELLE_PATTERN = re.compile(
    r"mitochond|chloroplast|d__eukary", flags=re.IGNORECASE
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    outdir = args.outdir.resolve()
    os.environ["TAINAN_RESULTS_DIR"] = str(outdir)
    sys.path.insert(0, str(project_root))

    import reanalyze_dependence_aware as analysis

    def validate_and_load_prokaryotic():
        meta = pd.read_csv(analysis.META_FILE)
        water = pd.read_csv(analysis.WATER_FILE)
        climate = pd.read_csv(analysis.CLIMATE_FILE)
        asv_all = pd.read_csv(analysis.ASV_FILE)

        if meta["SampleID"].duplicated().any():
            raise ValueError("Duplicate SampleID values in metadata.")
        key = ["Pond", "Date", "Layer"]
        check_m = (
            meta[
                key
                + ["Temperature", "Salinity", "pH", "DO", "Turbidity"]
            ]
            .sort_values(key)
            .reset_index(drop=True)
        )
        check_w = water.sort_values(key).reset_index(drop=True)
        if not check_m[key].equals(check_w[key]):
            raise ValueError("Water and metadata keys do not match.")
        numeric = ["Temperature", "Salinity", "pH", "DO", "Turbidity"]
        if not np.allclose(check_m[numeric], check_w[numeric], equal_nan=True):
            raise ValueError("Water values differ between Table_A and metadata.")

        sample_cols = [
            c for c in asv_all.columns if c in set(meta["SampleID"])
        ]
        if len(sample_cols) != len(meta):
            missing = sorted(set(meta["SampleID"]) - set(sample_cols))
            raise ValueError(f"ASV table is missing metadata samples: {missing}")

        taxonomy = asv_all["Taxon"].fillna("").astype(str)
        organelle = taxonomy.str.contains(ORGANELLE_PATTERN)
        positive = asv_all[sample_cols].fillna(0).sum(axis=1) > 0
        keep = positive & ~organelle
        asv = asv_all.loc[keep].reset_index(drop=True)

        counts = asv[sample_cols].T.to_numpy(dtype=float)
        rel = analysis.relative_abundance_rows(counts)
        rel_df = pd.DataFrame(
            rel,
            index=sample_cols,
            columns=asv["id"].astype(str),
        )

        meta = meta.set_index("SampleID").loc[sample_cols].reset_index()
        meta["PondMonth"] = (
            meta["Pond"].astype(str) + "-" + meta["Date"].astype(str)
        )
        meta["Date_dt"] = pd.to_datetime(
            meta["Date"].astype(str), format="%Y%m"
        )

        if "Season_verified" not in meta.columns:
            meta["Season_verified"] = meta["Season"]
        if int(meta["Season_verified"].notna().sum()) != 160:
            raise ValueError(
                "Expected 160 unambiguous sample-level season assignments."
            )

        unique_pressure = (
            meta[["Date", "MinPressure"]].drop_duplicates().dropna()
        )
        low_threshold = float(
            unique_pressure["MinPressure"].quantile(0.10)
        )
        meta["LowPressure"] = np.where(
            meta["MinPressure"].notna(),
            (meta["MinPressure"] <= low_threshold).astype(float),
            np.nan,
        )
        complete_cols = [
            "Heatwave_Tmax_q90_3day",
            "Typhoon",
            "ColdSurge",
            "LowPressure",
        ]
        meta["DisturbanceComplete"] = meta[complete_cols].notna().all(axis=1)
        meta["Shannon"] = np.apply_along_axis(
            analysis.shannon_from_counts, 1, counts
        )

        validate_and_load_prokaryotic.qc = {
            "ASV_records_source": int(len(asv_all)),
            "ASVs_positive_in_formal_samples_before_organelle_filter": int(
                positive.sum()
            ),
            "Organelle_or_eukaryotic_ASVs_removed": int(
                (positive & organelle).sum()
            ),
            "Prokaryotic_ASVs_analyzed": int(len(asv)),
            "Organelle_or_eukaryotic_read_fraction_removed": float(
                asv_all.loc[positive & organelle, sample_cols]
                .to_numpy(dtype=float)
                .sum()
                / asv_all.loc[positive, sample_cols]
                .to_numpy(dtype=float)
                .sum()
            ),
        }
        return meta, asv, rel_df, climate, low_threshold

    analysis.validate_and_load = validate_and_load_prokaryotic
    analysis.main()

    quality_path = outdir / "analysis_quality_report.json"
    report = json.loads(quality_path.read_text(encoding="utf-8"))
    report["feature_filter"] = validate_and_load_prokaryotic.qc
    report["analysis_label"] = (
        "Dependence-aware prokaryotic-community sensitivity analysis"
    )
    quality_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
