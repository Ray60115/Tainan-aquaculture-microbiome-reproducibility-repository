#!/usr/bin/env python3
"""Dependence-aware analysis of PICRUSt2-inferred MetaCyc pathway profiles.

The analytical unit is matched to the sampling design:

* descriptive functional PCoA: equal-weighted pond-month mean profiles;
* season: grow-out pond-months, controlling pond, with circular shifts within
  pond;
* sampling layer: complete surface-bottom pairs, controlling pond-month, with
  labels swapped only within pairs;
* meteorological indicators: pond-month profiles, controlling pond and season,
  with the shared monthly exposure sequence shifted across dates;
* taxonomic-functional turnover association: successive pond-month
  transitions, with a within-pond circular-shift test.

The pathway modules are prespecified in ``pathway_module_definitions.csv``.
They are summaries of *predicted functional potential*, not measurements of
gene expression, pathway activity, or biogeochemical process rates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.stats import rankdata, spearmanr, wilcoxon


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_PROJECT = SCRIPT_DIR.parent
DEFAULT_HELPERS = (
    SCRIPT_DIR.parent.parent
    / "delivery_complete_climate_20260728_v2"
    / "code_and_results"
    / "scripts"
)
if str(DEFAULT_HELPERS) not in sys.path:
    sys.path.insert(0, str(DEFAULT_HELPERS))

from reanalyze_dependence_aware import (  # noqa: E402
    bh_adjust,
    marginal_r2_distance,
    pcoa,
    restricted_dispersion_test,
    restricted_permanova,
    restricted_univariate_test,
)


SEED = 20260728
N_PERM = 9999
EVENT_COLUMNS = [
    "Typhoon",
    "ColdSurge",
    "Heatwave_Tmax_q90_3day",
    "LowPressure",
]
POND_COLORS = {
    "T1": "#6A51A3",
    "T2": "#D95F5F",
    "T3": "#E69F00",
    "T4": "#2C7FB8",
}
SEASON_ORDER = ["Spring", "Summer", "Autumn", "Winter"]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize_rows(frame: pd.DataFrame) -> pd.DataFrame:
    totals = frame.sum(axis=1)
    if (totals <= 0).any():
        raise ValueError(
            "Zero-sum pathway profiles: "
            + ", ".join(frame.index[totals <= 0].astype(str))
        )
    return frame.div(totals, axis=0)


def prepare_metadata(path: Path) -> pd.DataFrame:
    meta = pd.read_csv(path)
    required = {
        "SampleID",
        "Pond",
        "Date",
        "Layer",
        "Season",
        "Typhoon",
        "ColdSurge",
        "Heatwave_Tmax_q90_3day",
        "MinPressure",
    }
    missing = sorted(required - set(meta.columns))
    if missing:
        raise ValueError(f"Metadata columns missing: {missing}")
    if len(meta) != 168 or meta["SampleID"].duplicated().any():
        raise ValueError("Expected 168 unique formal SampleID rows.")

    meta = meta.copy()
    meta["Date"] = pd.to_numeric(meta["Date"], errors="raise").astype(int)
    meta["Date_dt"] = pd.to_datetime(meta["Date"].astype(str), format="%Y%m")
    meta["PondMonth"] = meta["Pond"].astype(str) + "-" + meta["Date"].astype(str)
    if "Season_verified" not in meta.columns:
        meta["Season_verified"] = meta["Season"]
    # Missing values are intentional for averaged profiles whose contributing
    # occasions straddle a solar-term season boundary.  They are retained for
    # non-seasonal analyses and excluded below only where season is required.
    if int(meta["Season_verified"].notna().sum()) != 160:
        raise ValueError("Expected 160 unambiguous sample-level season labels.")

    monthly_pressure = meta[["Date", "MinPressure"]].drop_duplicates()
    if monthly_pressure.groupby("Date")["MinPressure"].nunique(dropna=False).max() != 1:
        raise ValueError("MinPressure is not shared within sampling month.")
    pressure_q10 = float(monthly_pressure["MinPressure"].quantile(0.10))
    meta["LowPressure"] = (meta["MinPressure"] <= pressure_q10).astype(int)
    if not set(meta["Layer"].astype(str)).issubset({"U", "D"}):
        raise ValueError("Only U and D layers are permitted in the formal cohort.")
    return meta


def read_pathway_profiles(
    path: Path, sample_ids: list[str]
) -> tuple[pd.DataFrame, pd.Series]:
    raw = pd.read_csv(path, sep="\t", compression="infer", index_col=0)
    raw.index = raw.index.astype(str)
    if raw.index.duplicated().any():
        raise ValueError("Duplicated pathway identifiers in PICRUSt2 output.")
    missing = sorted(set(sample_ids) - set(raw.columns))
    extra = sorted(set(raw.columns) - set(sample_ids))
    if missing or extra:
        raise ValueError(
            f"Pathway/sample mismatch: missing={missing[:8]}, extra={extra[:8]}"
        )
    raw = raw.loc[:, sample_ids].apply(pd.to_numeric, errors="raise")
    if (raw < 0).any().any():
        raise ValueError("Negative inferred pathway abundances detected.")
    raw = raw.loc[raw.sum(axis=1) > 0]
    profiles = normalize_rows(raw.T)
    prevalence = (raw > 0).mean(axis=1)
    return profiles, prevalence


def aggregate_pond_month(
    meta: pd.DataFrame, profiles: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    indexed = meta.set_index("SampleID")
    joined = profiles.join(indexed[["PondMonth"]])
    pm_profiles = joined.groupby("PondMonth").mean(numeric_only=True)
    pm_profiles = normalize_rows(pm_profiles)
    pm = (
        meta.groupby("PondMonth", as_index=False)
        .agg(
            Pond=("Pond", "first"),
            Date=("Date", "first"),
            Date_dt=("Date_dt", "first"),
            Season=("Season_verified", "first"),
            n_layers=("Layer", "nunique"),
            Typhoon=("Typhoon", "first"),
            ColdSurge=("ColdSurge", "first"),
            Heatwave_Tmax_q90_3day=("Heatwave_Tmax_q90_3day", "first"),
            LowPressure=("LowPressure", "first"),
        )
        .set_index("PondMonth")
        .loc[pm_profiles.index]
        .reset_index()
    )
    if len(pm) != 93:
        raise ValueError(f"Expected 93 pond-months, found {len(pm)}.")
    return pm, pm_profiles.loc[pm["PondMonth"]]


def subset_distance(frame: pd.DataFrame, ids: list[str]) -> np.ndarray:
    return squareform(
        pdist(frame.loc[ids].to_numpy(dtype=float), metric="braycurtis")
    )


def run_profile_tests(
    meta: pd.DataFrame,
    profiles: pd.DataFrame,
    pm: pd.DataFrame,
    pm_profiles: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict] = []
    dispersion_rows: list[dict] = []

    for dataset, frame in [
        ("All ponds", pm.copy()),
        ("Grow-out ponds", pm[pm["Pond"].ne("T1")].copy()),
    ]:
        frame = frame.reset_index(drop=True)
        distance = subset_distance(pm_profiles, frame["PondMonth"].tolist())
        rows.append(
            {
                "Response": "Predicted MetaCyc pathway profile",
                "Dataset": dataset,
                "Factor": "Pond identity",
                "Unit": "Pond-month mean",
                "n": len(frame),
                "Pseudo_F": np.nan,
                "Marginal_R2": marginal_r2_distance(distance, frame["Pond"]),
                "Partial_R2": np.nan,
                "p_value": np.nan,
                "q_value": np.nan,
                "Permutations": 0,
                "Restriction": "Descriptive only; pond identity was not permuted",
            }
        )

    grow = (
        pm[pm["Pond"].ne("T1") & pm["Season"].notna()]
        .sort_values(["Pond", "Date"])
        .reset_index(drop=True)
    )
    grow_distance = subset_distance(pm_profiles, grow["PondMonth"].tolist())
    result = restricted_permanova(
        grow,
        grow_distance,
        "Season",
        ["Pond"],
        "within_pond_circular",
        n_perm=N_PERM,
    )
    rows.append(
        {
            "Response": "Predicted MetaCyc pathway profile",
            "Dataset": "Grow-out ponds",
            "Factor": "Season",
            "Unit": "Pond-month mean",
            "n": result["n"],
            "Pseudo_F": result["F"],
            "Marginal_R2": result["marginal_R2"],
            "Partial_R2": result["partial_R2"],
            "p_value": result["p_value"],
            "q_value": np.nan,
            "Permutations": result["permutations"],
            "Restriction": "Circular shifts of season labels within pond",
        }
    )
    dispersion = restricted_dispersion_test(
        grow,
        grow_distance,
        "Season",
        "within_pond_circular",
        n_perm=N_PERM,
    )
    dispersion_rows.append(
        {
            "Dataset": "Grow-out ponds",
            "Factor": "Season",
            "n": dispersion["n"],
            "F": dispersion["F"],
            "p_value": dispersion["p_value"],
            "Permutations": dispersion["permutations"],
            "Restriction": "Circular shifts of season labels within pond",
        }
    )

    pair_index = meta.groupby("PondMonth")["Layer"].nunique()
    pair_ids = pair_index[pair_index.eq(2)].index
    layer_meta = (
        meta[meta["PondMonth"].isin(pair_ids)]
        .sort_values(["PondMonth", "Layer"])
        .reset_index(drop=True)
    )
    layer_distance = squareform(
        pdist(
            profiles.loc[layer_meta["SampleID"]].to_numpy(dtype=float),
            metric="braycurtis",
        )
    )
    result = restricted_permanova(
        layer_meta,
        layer_distance,
        "Layer",
        ["PondMonth"],
        "paired_layer_swap",
        n_perm=N_PERM,
    )
    rows.append(
        {
            "Response": "Predicted MetaCyc pathway profile",
            "Dataset": "All ponds",
            "Factor": "Sampling layer",
            "Unit": "75 complete surface-bottom pairs",
            "n": result["n"],
            "Pseudo_F": result["F"],
            "Marginal_R2": result["marginal_R2"],
            "Partial_R2": result["partial_R2"],
            "p_value": result["p_value"],
            "q_value": np.nan,
            "Permutations": result["permutations"],
            "Restriction": "Surface/bottom labels swapped only within pond-month",
        }
    )
    dispersion = restricted_dispersion_test(
        layer_meta,
        layer_distance,
        "Layer",
        "paired_layer_swap",
        n_perm=N_PERM,
    )
    dispersion_rows.append(
        {
            "Dataset": "All ponds",
            "Factor": "Sampling layer",
            "n": dispersion["n"],
            "F": dispersion["F"],
            "p_value": dispersion["p_value"],
            "Permutations": dispersion["permutations"],
            "Restriction": "Surface/bottom labels swapped only within pond-month",
        }
    )

    event_row_indices: dict[str, list[int]] = {"All ponds": [], "Grow-out ponds": []}
    for dataset, frame in [
        ("All ponds", pm.copy()),
        ("Grow-out ponds", pm[pm["Pond"].ne("T1")].copy()),
    ]:
        frame = frame[frame["Season"].notna()].sort_values(
            ["Pond", "Date"]
        ).reset_index(drop=True)
        distance = subset_distance(pm_profiles, frame["PondMonth"].tolist())
        for factor in EVENT_COLUMNS:
            result = restricted_permanova(
                frame,
                distance,
                factor,
                ["Pond", "Season"],
                "date_circular",
                n_perm=N_PERM,
            )
            n_event = int(pd.to_numeric(frame[factor], errors="raise").sum())
            n_non_event = int(len(frame) - n_event)
            rows.append(
                {
                    "Response": "Predicted MetaCyc pathway profile",
                    "Dataset": dataset,
                    "Factor": factor,
                    "Unit": "Pond-month mean; sampling month permuted as one block",
                    "n": result["n"],
                    "n_event": n_event,
                    "n_non_event": n_non_event,
                    "Pseudo_F": result["F"],
                    "Marginal_R2": result["marginal_R2"],
                    "Partial_R2": result["partial_R2"],
                    "p_value": result["p_value"],
                    "q_value": np.nan,
                    "Permutations": result["permutations"],
                    "Restriction": "Circular shifts of the shared monthly exposure sequence",
                }
            )
            event_row_indices[dataset].append(len(rows) - 1)

    for dataset, indices in event_row_indices.items():
        adjusted = bh_adjust([rows[i]["p_value"] for i in indices])
        for index, q_value in zip(indices, adjusted, strict=True):
            rows[index]["q_value"] = q_value

    return pd.DataFrame(rows), pd.DataFrame(dispersion_rows)


def build_module_scores(
    definitions: pd.DataFrame,
    profiles: pd.DataFrame,
    prevalence: pd.Series,
    pathway_names: pd.Series,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    definitions = definitions.copy()
    definitions["Observed_MetaCyc_name"] = definitions["Pathway_ID"].map(pathway_names)
    definitions["Detected_in_primary_output"] = definitions["Pathway_ID"].isin(
        profiles.columns
    )
    definitions["Sample_prevalence"] = (
        definitions["Pathway_ID"].map(prevalence).fillna(0.0)
    )
    definitions["Name_matches_expected"] = (
        definitions["Observed_MetaCyc_name"]
        .fillna("")
        .str.casefold()
        .eq(definitions["Expected_MetaCyc_name"].str.casefold())
    )

    transformed = np.log1p(profiles * 1_000_000.0)
    means = transformed.mean(axis=0)
    standard_deviation = transformed.std(axis=0, ddof=0)
    stable = standard_deviation > 0
    z = transformed.loc[:, stable].sub(means[stable], axis=1).div(
        standard_deviation[stable], axis=1
    )

    scores = pd.DataFrame(index=profiles.index)
    included_flag = []
    for module, group in definitions.groupby("Module", sort=False):
        pathways = [
            pathway
            for pathway in group["Pathway_ID"]
            if pathway in z.columns
        ]
        if len(pathways) < 1:
            raise ValueError(f"No detected, variable pathways for module {module!r}.")
        scores[module] = z[pathways].mean(axis=1)
        included_flag.extend(
            [
                (module, pathway, pathway in pathways)
                for pathway in group["Pathway_ID"]
            ]
        )
    flag = pd.DataFrame(
        included_flag,
        columns=["Module", "Pathway_ID", "Included_in_module_score"],
    )
    definitions = definitions.merge(
        flag, on=["Module", "Pathway_ID"], how="left", validate="one_to_one"
    )
    return scores, definitions


def run_module_tests(
    meta: pd.DataFrame, module_scores: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sample = meta.set_index("SampleID").join(module_scores)
    modules = module_scores.columns.tolist()
    pond_month_scores = sample.groupby("PondMonth")[modules].mean()
    pm_meta = (
        meta.groupby("PondMonth", as_index=False)
        .agg(
            Pond=("Pond", "first"),
            Date=("Date", "first"),
            Season=("Season_verified", "first"),
        )
        .set_index("PondMonth")
        .join(pond_month_scores)
        .reset_index()
    )

    season_rows = []
    grow = (
        pm_meta[pm_meta["Pond"].ne("T1") & pm_meta["Season"].notna()]
        .sort_values(["Pond", "Date"])
        .reset_index(drop=True)
    )
    for module in modules:
        result = restricted_univariate_test(
            grow,
            module,
            "Season",
            ["Pond"],
            "within_pond_circular",
            n_perm=N_PERM,
        )
        season_rows.append(
            {
                "Response": module,
                "Factor": "Season",
                "Dataset": "Grow-out ponds",
                "Unit": "Pond-month mean module score",
                "n": result["n"],
                "F": result["F"],
                "Partial_R2": result["partial_R2"],
                "p_value": result["p_value"],
                "q_value": np.nan,
                "Permutations": result["permutations"],
                "Restriction": "Circular shifts of season labels within pond",
            }
        )
    q_values = bh_adjust([row["p_value"] for row in season_rows])
    for row, q_value in zip(season_rows, q_values, strict=True):
        row["q_value"] = q_value

    pair_ids = (
        meta.groupby("PondMonth")["Layer"]
        .nunique()
        .loc[lambda x: x.eq(2)]
        .index
    )
    paired = sample[sample["PondMonth"].isin(pair_ids)].copy()
    layer_rows = []
    for module in modules:
        wide = paired.pivot(
            index="PondMonth", columns="Layer", values=module
        ).dropna(subset=["D", "U"])
        test = wilcoxon(
            wide["D"], wide["U"], alternative="two-sided", method="auto"
        )
        layer_rows.append(
            {
                "Response": module,
                "Factor": "Sampling layer",
                "Dataset": "All ponds",
                "Unit": "Complete surface-bottom pond-month pair",
                "n_pairs": len(wide),
                "Median_paired_difference_bottom_minus_surface": float(
                    np.median(wide["D"] - wide["U"])
                ),
                "Wilcoxon_W": float(test.statistic),
                "p_value": float(test.pvalue),
                "q_value": np.nan,
                "Restriction": "Explicit within-pond-month pairing",
            }
        )
    q_values = bh_adjust([row["p_value"] for row in layer_rows])
    for row, q_value in zip(layer_rows, q_values, strict=True):
        row["q_value"] = q_value

    return (
        pd.DataFrame(season_rows),
        pd.DataFrame(layer_rows),
        pm_meta,
    )


def calculate_turnover(
    pm: pd.DataFrame, pm_profiles: pd.DataFrame
) -> pd.DataFrame:
    rows = []
    for pond, group in pm.sort_values("Date").groupby("Pond", sort=False):
        ordered = group.sort_values("Date")
        for position in range(1, len(ordered)):
            previous = ordered.iloc[position - 1]
            current = ordered.iloc[position]
            a = pm_profiles.loc[previous["PondMonth"]].to_numpy(dtype=float)
            b = pm_profiles.loc[current["PondMonth"]].to_numpy(dtype=float)
            rows.append(
                {
                    "Pond": pond,
                    "PreviousDate": int(previous["Date"]),
                    "Date": int(current["Date"]),
                    "Season": current["Season"],
                    "Functional_turnover": float(
                        np.abs(a - b).sum() / (a + b).sum()
                    ),
                }
            )
    return pd.DataFrame(rows)


def residual_rank_correlation(frame: pd.DataFrame) -> float:
    x = rankdata(frame["Taxonomic_turnover"].to_numpy(dtype=float))
    y = rankdata(frame["Functional_turnover"].to_numpy(dtype=float))
    pond = pd.get_dummies(frame["Pond"], drop_first=True, dtype=float)
    design = np.column_stack([np.ones(len(frame)), pond.to_numpy(dtype=float)])
    x_resid = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    y_resid = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    return float(np.corrcoef(x_resid, y_resid)[0, 1])


def turnover_association(
    functional: pd.DataFrame, taxonomic_path: Path
) -> tuple[pd.DataFrame, dict]:
    taxonomic = pd.read_csv(taxonomic_path).rename(
        columns={"Turnover": "Taxonomic_turnover"}
    )
    keys = ["Pond", "PreviousDate", "Date"]
    merged = taxonomic.merge(
        functional, on=keys, how="inner", validate="one_to_one"
    )
    if len(merged) != 89:
        raise ValueError(f"Expected 89 matched turnover transitions, found {len(merged)}.")

    raw = spearmanr(
        merged["Taxonomic_turnover"], merged["Functional_turnover"]
    )
    observed = residual_rank_correlation(merged)
    rng = np.random.default_rng(SEED)
    groups = []
    for _, indices in merged.groupby("Pond", sort=False).groups.items():
        idx = np.asarray(list(indices), dtype=int)
        idx = idx[np.argsort(merged.loc[idx, "Date"].to_numpy())]
        groups.append(idx)
    exceed = 0
    for _ in range(N_PERM):
        permuted = merged.copy()
        values = merged["Functional_turnover"].to_numpy(copy=True)
        for indices in groups:
            shift = int(rng.integers(0, len(indices)))
            values[indices] = np.roll(
                merged.loc[indices, "Functional_turnover"].to_numpy(), shift
            )
        permuted["Functional_turnover"] = values
        correlation = residual_rank_correlation(permuted)
        exceed += int(abs(correlation) >= abs(observed) - 1e-12)

    summary = {
        "n_transitions": len(merged),
        "n_ponds": int(merged["Pond"].nunique()),
        "raw_spearman_rho": float(raw.statistic),
        "raw_spearman_p_unrestricted_descriptive": float(raw.pvalue),
        "pond_adjusted_rank_correlation": observed,
        "restricted_two_sided_p": (exceed + 1) / (N_PERM + 1),
        "permutations": N_PERM,
        "restriction": "Independent circular shifts of functional-turnover values within each pond",
        "interpretation_guardrail": (
            "The predicted functional profiles are derived from the same 16S "
            "composition and therefore do not provide independent validation "
            "of functional activity."
        ),
    }
    return merged, summary


def load_pathway_names(path: Path) -> pd.Series:
    if path.suffix.lower() == ".csv":
        audited = pd.read_csv(path)
        required = {"Pathway_ID", "Observed_MetaCyc_name"}
        if not required.issubset(audited.columns):
            raise ValueError(
                "CSV pathway-name map must contain Pathway_ID and "
                "Observed_MetaCyc_name."
            )
        names = audited.drop_duplicates("Pathway_ID").set_index("Pathway_ID")[
            "Observed_MetaCyc_name"
        ]
        names.index = names.index.astype(str)
        return names
    names = pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        header=None,
        index_col=0,
        names=["Pathway_ID", "MetaCyc_name"],
    )["MetaCyc_name"]
    names.index = names.index.astype(str)
    return names


def add_confidence_ellipse(
    ax: plt.Axes, x: np.ndarray, y: np.ndarray, color: str
) -> None:
    if len(x) < 3:
        return
    covariance = np.cov(x, y)
    values, vectors = np.linalg.eigh(covariance)
    order = values.argsort()[::-1]
    values = values[order]
    vectors = vectors[:, order]
    angle = math.degrees(math.atan2(vectors[1, 0], vectors[0, 0]))
    # 80% bivariate-normal contour: sqrt(chi2.ppf(0.80, 2)).
    scale = 1.7941225779941015
    width, height = 2 * scale * np.sqrt(np.maximum(values, 0))
    ellipse = Ellipse(
        (float(np.mean(x)), float(np.mean(y))),
        width,
        height,
        angle=angle,
        facecolor="none",
        edgecolor=color,
        linewidth=1.4,
        alpha=0.75,
    )
    ax.add_patch(ellipse)


def make_main_figure(
    pm: pd.DataFrame,
    pm_profiles: pd.DataFrame,
    turnover: pd.DataFrame,
    turnover_summary: dict,
    pm_modules: pd.DataFrame,
    output: Path,
) -> pd.DataFrame:
    distance = squareform(
        pdist(pm_profiles.to_numpy(dtype=float), metric="braycurtis")
    )
    coordinates, _, explained = pcoa(distance)
    coord = pm[["PondMonth", "Pond", "Date", "Season"]].copy()
    coord["PCoA1"] = coordinates[:, 0]
    coord["PCoA2"] = coordinates[:, 1]

    grow_modules = pm_modules[pm_modules["Pond"].ne("T1")]
    modules = [
        column
        for column in pm_modules.columns
        if column not in {"PondMonth", "Pond", "Date", "Season"}
    ]
    heatmap = (
        grow_modules.groupby("Season")[modules]
        .mean()
        .reindex(SEASON_ORDER)
        .T
    )

    sns.set_theme(style="white", context="paper", font_scale=1.05)
    fig = plt.figure(figsize=(16.2, 5.2), constrained_layout=True)
    grid = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.0, 1.38])
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[0, 2])

    for pond in ["T1", "T2", "T3", "T4"]:
        sub = coord[coord["Pond"].eq(pond)].sort_values("Date")
        ax_a.plot(
            sub["PCoA1"],
            sub["PCoA2"],
            color=POND_COLORS[pond],
            linewidth=0.65,
            alpha=0.22,
            zorder=1,
        )
        ax_a.scatter(
            sub["PCoA1"],
            sub["PCoA2"],
            s=28,
            color=POND_COLORS[pond],
            alpha=0.82,
            edgecolor="white",
            linewidth=0.35,
            label=pond,
            zorder=2,
        )
        add_confidence_ellipse(
            ax_a,
            sub["PCoA1"].to_numpy(),
            sub["PCoA2"].to_numpy(),
            POND_COLORS[pond],
        )
    ax_a.set_title("(a) Predicted functional-profile PCoA", loc="left", fontweight="bold")
    ax_a.set_xlabel(f"PCoA1 ({explained[0]:.1f}%)")
    ax_a.set_ylabel(f"PCoA2 ({explained[1]:.1f}%)")
    ax_a.legend(title="Pond", frameon=False, ncol=2, loc="best")

    for pond in ["T1", "T2", "T3", "T4"]:
        sub = turnover[turnover["Pond"].eq(pond)]
        ax_b.scatter(
            sub["Taxonomic_turnover"],
            sub["Functional_turnover"],
            s=32,
            color=POND_COLORS[pond],
            alpha=0.78,
            edgecolor="white",
            linewidth=0.35,
            label=pond,
        )
    x = turnover["Taxonomic_turnover"].to_numpy(dtype=float)
    y = turnover["Functional_turnover"].to_numpy(dtype=float)
    fit = np.polyfit(x, y, 1)
    x_line = np.linspace(x.min(), x.max(), 100)
    ax_b.plot(x_line, np.polyval(fit, x_line), color="#333333", linewidth=1.2)
    ax_b.text(
        0.04,
        0.96,
        (
            "Pond-adjusted rank correlation\n"
            f"$r$ = {turnover_summary['pond_adjusted_rank_correlation']:.3f}; "
            f"restricted $p$ = {turnover_summary['restricted_two_sided_p']:.4f}"
        ),
        transform=ax_b.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.86, "pad": 2},
    )
    ax_b.set_title("(b) Consecutive pond-month turnover", loc="left", fontweight="bold")
    ax_b.set_xlabel("Taxonomic Bray–Curtis dissimilarity")
    ax_b.set_ylabel("Predicted-functional Bray–Curtis dissimilarity")

    labels = [label.replace(" and ", " & ") for label in heatmap.index]
    sns.heatmap(
        heatmap,
        cmap="vlag",
        center=0,
        annot=True,
        fmt=".2f",
        linewidths=0.6,
        linecolor="white",
        cbar=False,
        ax=ax_c,
        yticklabels=labels,
    )
    ax_c.set_title("(c) Grow-out seasonal module profiles", loc="left", fontweight="bold")
    ax_c.set_xlabel("Season")
    ax_c.set_ylabel("")
    ax_c.tick_params(axis="y", labelrotation=0)

    for axis in [ax_a, ax_b]:
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.grid(color="#E6E6E6", linewidth=0.6, alpha=0.8)
    fig.savefig(output, dpi=400, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    coord.attrs["PCoA1_explained"] = float(explained[0])
    coord.attrs["PCoA2_explained"] = float(explained[1])
    return coord


def collect_picrust_qc(
    picrust_dir: Path,
    abundance_input: Path,
    sample_ids: list[str],
) -> tuple[pd.DataFrame, dict]:
    weighted_path = picrust_dir / "EC_metagenome_out" / "weighted_nsti.tsv.gz"
    marker_path = picrust_dir / "combined_marker_predicted_and_nsti.tsv.gz"
    normalized_path = picrust_dir / "EC_metagenome_out" / "seqtab_norm.tsv.gz"
    required = [weighted_path, marker_path, normalized_path]
    if not all(path.exists() for path in required):
        return pd.DataFrame({"SampleID": sample_ids}), {
            "status": "not repeated for derived blank-screened sensitivity output",
            "missing_files": [path.name for path in required if not path.exists()],
            "note": (
                "Placement and NSTI quality metrics are reported for the primary "
                "analysis; the blank-screened run is used only to compare pathway "
                "profiles and test classifications."
            ),
        }
    weighted = pd.read_csv(weighted_path, sep="\t", index_col=0)
    marker = pd.read_csv(marker_path, sep="\t", index_col=0)
    normalized = pd.read_csv(normalized_path, sep="\t", index_col=0)
    abundance = pd.read_csv(abundance_input, sep="\t", index_col=0)
    abundance.index = abundance.index.astype(str)
    marker.index = marker.index.astype(str)
    marker_for_input = marker.loc[marker.index.intersection(abundance.index)]
    retained_ids = abundance.index.intersection(normalized.index.astype(str))
    raw_total = abundance[sample_ids].sum(axis=0)
    retained_total = abundance.loc[retained_ids, sample_ids].sum(axis=0)

    weighted_series = weighted.iloc[:, 0].reindex(sample_ids)
    if weighted_series.isna().any():
        raise ValueError("Weighted NSTI output does not contain all 168 samples.")
    sample_qc = pd.DataFrame(
        {
            "SampleID": sample_ids,
            "Input_prokaryotic_reads": raw_total.reindex(sample_ids).to_numpy(),
            "Reads_from_ASVs_retained_after_placement_and_NSTI": retained_total.reindex(
                sample_ids
            ).to_numpy(),
            "Retained_read_fraction": (
                retained_total.reindex(sample_ids) / raw_total.reindex(sample_ids)
            ).to_numpy(),
            "Weighted_NSTI": weighted_series.to_numpy(dtype=float),
        }
    )
    summary = {
        "input_ASVs": int(len(abundance)),
        "placed_ASVs_with_domain_assignment": int(len(marker_for_input)),
        "ASVs_retained_in_metagenome_step": int(len(retained_ids)),
        "best_domain_counts": marker_for_input[
            "best_domain"
        ].value_counts().to_dict(),
        "ASVs_above_default_NSTI_2": int(
            (pd.to_numeric(marker_for_input["metadata_NSTI"]) > 2).sum()
        ),
        "weighted_NSTI": {
            "minimum": float(sample_qc["Weighted_NSTI"].min()),
            "q1": float(sample_qc["Weighted_NSTI"].quantile(0.25)),
            "median": float(sample_qc["Weighted_NSTI"].median()),
            "q3": float(sample_qc["Weighted_NSTI"].quantile(0.75)),
            "maximum": float(sample_qc["Weighted_NSTI"].max()),
        },
        "retained_read_fraction": {
            "minimum": float(sample_qc["Retained_read_fraction"].min()),
            "q1": float(sample_qc["Retained_read_fraction"].quantile(0.25)),
            "median": float(sample_qc["Retained_read_fraction"].median()),
            "q3": float(sample_qc["Retained_read_fraction"].quantile(0.75)),
            "maximum": float(sample_qc["Retained_read_fraction"].max()),
        },
    }
    return sample_qc, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, default=DEFAULT_PROJECT)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--picrust", required=True, type=Path)
    parser.add_argument("--abundance-input", required=True, type=Path)
    parser.add_argument("--taxonomy-turnover", required=True, type=Path)
    parser.add_argument("--module-definitions", required=True, type=Path)
    parser.add_argument("--pathway-name-map", required=True, type=Path)
    args = parser.parse_args()

    project = args.project.resolve()
    results = project / "results"
    tables = project / "tables"
    figures = project / "figures"
    qc_dir = project / "qc"
    for directory in [results, tables, figures, qc_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    meta = prepare_metadata(args.metadata)
    sample_ids = meta["SampleID"].tolist()
    pathway_path = args.picrust / "pathways_out" / "path_abun_unstrat.tsv.gz"
    profiles, prevalence = read_pathway_profiles(pathway_path, sample_ids)
    pathway_names = load_pathway_names(args.pathway_name_map)
    pm, pm_profiles = aggregate_pond_month(meta, profiles)

    profile_tests, dispersion_tests = run_profile_tests(
        meta, profiles, pm, pm_profiles
    )
    definitions = pd.read_csv(args.module_definitions)
    module_scores, audited_definitions = build_module_scores(
        definitions, profiles, prevalence, pathway_names
    )
    season_modules, layer_modules, pm_modules = run_module_tests(
        meta, module_scores
    )

    functional_turnover = calculate_turnover(pm, pm_profiles)
    turnover, turnover_summary = turnover_association(
        functional_turnover, args.taxonomy_turnover
    )
    coord = make_main_figure(
        pm,
        pm_profiles,
        turnover,
        turnover_summary,
        pm_modules,
        figures / "Fig10_predicted_functional_potential.png",
    )
    sample_qc, picrust_qc = collect_picrust_qc(
        args.picrust, args.abundance_input, sample_ids
    )

    profile_tests.to_csv(tables / "functional_profile_dependence_aware_tests.csv", index=False)
    dispersion_tests.to_csv(
        tables / "functional_profile_permdisp_tests.csv", index=False
    )
    audited_definitions.to_csv(
        tables / "Table_S11_pathway_module_definitions.csv", index=False
    )
    season_modules.to_csv(
        tables / "functional_module_season_tests.csv", index=False
    )
    layer_modules.to_csv(
        tables / "functional_module_layer_tests.csv", index=False
    )
    pm_modules.to_csv(
        tables / "functional_module_scores_pond_month.csv", index=False
    )
    module_scores.to_csv(
        tables / "functional_module_scores_sample.csv", index_label="SampleID"
    )
    turnover.to_csv(
        tables / "taxonomic_functional_turnover_comparison.csv", index=False
    )
    coord.to_csv(
        tables / "functional_pcoa_coordinates_pond_month.csv", index=False
    )
    sample_qc.to_csv(qc_dir / "picrust2_sample_qc.csv", index=False)

    combined_s14 = pd.concat(
        [
            profile_tests.assign(Test_family="Whole predicted pathway profile"),
            season_modules.rename(
                columns={"F": "Pseudo_F", "Partial_R2": "Partial_R2"}
            ).assign(
                Test_family="Prespecified module: season",
                Marginal_R2=np.nan,
            ),
            layer_modules.rename(
                columns={
                    "n_pairs": "n",
                    "Wilcoxon_W": "Pseudo_F",
                    "Median_paired_difference_bottom_minus_surface": "Partial_R2",
                }
            ).assign(
                Test_family="Prespecified module: layer",
                Marginal_R2=np.nan,
                Permutations=0,
            ),
        ],
        ignore_index=True,
        sort=False,
    )
    combined_s14.to_csv(
        tables / "Table_S10_dependence_aware_functional_tests.csv", index=False
    )

    summary = {
        "formal_samples": len(meta),
        "pond_months": len(pm),
        "growout_pond_months_total": int(pm["Pond"].ne("T1").sum()),
        "growout_pond_months_in_season_models": int(
            (pm["Pond"].ne("T1") & pm["Season"].notna()).sum()
        ),
        "predicted_pathways": int(profiles.shape[1]),
        "prespecified_modules": int(definitions["Module"].nunique()),
        "profile_pcoa_positive_axis_variance_percent": {
            "PCoA1": float(coord.attrs["PCoA1_explained"]),
            "PCoA2": float(coord.attrs["PCoA2_explained"]),
        },
        "turnover_association": turnover_summary,
        "picrust2_qc": picrust_qc,
        "input_sha256": {
            pathway_path.name: sha256_file(pathway_path),
            args.module_definitions.name: sha256_file(args.module_definitions),
            args.metadata.name: sha256_file(args.metadata),
        },
        "interpretation": (
            "All functional results represent PICRUSt2-inferred genomic "
            "potential from 16S profiles, not measured genes, transcripts, "
            "proteins, metabolites, process rates, or pathway activity."
        ),
    }
    (results / "functional_analysis_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
