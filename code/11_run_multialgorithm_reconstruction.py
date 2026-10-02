#!/usr/bin/env python3
"""Grouped-CV sensitivity across six regression algorithms and two cohorts.

This script implements Supplemental Table S12 of the submitted manuscript.
The PCoA coordinate system is always derived from all 168 filtered profiles.
Model fitting is then performed for either the 160-sample all-pond cohort or
the 122-sample T2--T4 cohort, using identical fivefold GroupKFold splits by
pond-month within a cohort.

XGBoost and pyGAM are optional imports at program start but are required when
their corresponding algorithms are selected.  See ``requirements.txt`` for
the manuscript environment.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR


SEED = 42
N_FOLDS = 5
ORGANELLE_EUKARYOTE_PATTERN = r"mitochond|chloroplast|d__eukary"


FEATURE_SETS = {
    "Structure": ["Pond", "Layer"],
    "+ Season": [
        "Pond",
        "Layer",
        "Season_verified",
        "Month_sin",
        "Month_cos",
        "Year",
    ],
    "+ Water quality": [
        "Pond",
        "Layer",
        "Season_verified",
        "Month_sin",
        "Month_cos",
        "Year",
        "Temperature",
        "Salinity",
        "pH",
        "DO",
        "Turbidity",
    ],
    "+ Disturbances": [
        "Pond",
        "Layer",
        "Season_verified",
        "Month_sin",
        "Month_cos",
        "Year",
        "Temperature",
        "Salinity",
        "pH",
        "DO",
        "Turbidity",
        "Heatwave_Tmax_q90_3day",
        "Typhoon",
        "ColdSurge",
        "LowPressure",
    ],
}


def relative_abundance_rows(counts: np.ndarray) -> np.ndarray:
    totals = counts.sum(axis=1, keepdims=True)
    return np.divide(
        counts,
        totals,
        out=np.zeros_like(counts, dtype=float),
        where=totals > 0,
    )


def pcoa(distance_matrix: np.ndarray) -> np.ndarray:
    n = distance_matrix.shape[0]
    centering = np.eye(n) - np.ones((n, n)) / n
    gower = -0.5 * centering @ (distance_matrix**2) @ centering
    values, vectors = np.linalg.eigh(gower)
    order = np.argsort(values)[::-1]
    values = values[order]
    vectors = vectors[:, order]
    positive = values > 0
    return vectors[:, positive] * np.sqrt(values[positive])


def load_analysis_data(asv_path: Path, metadata_path: Path) -> pd.DataFrame:
    metadata = pd.read_csv(metadata_path)
    asv = pd.read_csv(asv_path, low_memory=False)
    samples = metadata["SampleID"].astype(str).tolist()
    missing = [sample for sample in samples if sample not in asv.columns]
    if missing:
        raise ValueError(f"ASV table is missing sample profiles: {missing}")

    taxonomy = asv["Taxon"].fillna("").astype(str)
    positive = asv[samples].fillna(0).sum(axis=1) > 0
    excluded = taxonomy.str.contains(
        ORGANELLE_EUKARYOTE_PATTERN, case=False, regex=True
    )
    retained = asv.loc[positive & ~excluded]
    if len(retained) != 23896:
        raise ValueError(f"Expected 23,896 retained ASVs; found {len(retained):,}")

    counts = retained[samples].T.to_numpy(dtype=float)
    distance = squareform(pdist(relative_abundance_rows(counts), metric="braycurtis"))
    coordinates = pcoa(distance)
    metadata = metadata.copy()
    metadata["PCoA1"] = coordinates[:, 0]
    metadata["PCoA2"] = coordinates[:, 1]
    metadata = metadata.loc[
        metadata["DisturbanceComplete"].astype(bool)
        & metadata["Season_verified"].notna()
    ].copy()
    month = metadata["Date"].astype(str).str[4:].astype(int)
    metadata["Month"] = month
    metadata["Month_sin"] = np.sin(2 * np.pi * month / 12)
    metadata["Month_cos"] = np.cos(2 * np.pi * month / 12)
    metadata["Year"] = metadata["Date"].astype(str).str[:4].astype(int)
    return metadata.reset_index(drop=True)


def encoded_features(data: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    categorical = [
        x for x in ("Pond", "Layer", "Season_verified") if x in columns
    ]
    return pd.get_dummies(
        data[columns], columns=categorical, drop_first=False, dtype=float
    )


def mean_only_predictions(
    data: pd.DataFrame, splits: list[tuple[np.ndarray, np.ndarray]]
) -> np.ndarray:
    y = data[["PCoA1", "PCoA2"]].to_numpy(dtype=float)
    predictions = np.full_like(y, np.nan)
    for train, test in splits:
        predictions[test] = y[train].mean(axis=0)
    return predictions


def standard_model_predictions(
    data: pd.DataFrame,
    columns: list[str],
    splits: list[tuple[np.ndarray, np.ndarray]],
    factory: Callable[[], object],
    *,
    scale_x: bool = False,
    scale_y: bool = False,
) -> np.ndarray:
    x = encoded_features(data, columns).to_numpy(dtype=float)
    y = data[["PCoA1", "PCoA2"]].to_numpy(dtype=float)
    predictions = np.full_like(y, np.nan)
    for train, test in splits:
        x_train, x_test = x[train], x[test]
        if scale_x:
            x_scaler = StandardScaler().fit(x_train)
            x_train = x_scaler.transform(x_train)
            x_test = x_scaler.transform(x_test)
        for axis in range(2):
            y_train = y[train, axis]
            if scale_y:
                y_scaler = StandardScaler().fit(y_train[:, None])
                fitted_y = y_scaler.transform(y_train[:, None]).ravel()
            else:
                y_scaler = None
                fitted_y = y_train
            model = factory()
            model.fit(x_train, fitted_y)
            estimate = np.asarray(model.predict(x_test)).reshape(-1)
            if y_scaler is not None:
                estimate = y_scaler.inverse_transform(estimate[:, None]).ravel()
            predictions[test, axis] = estimate
    return predictions


@dataclass(frozen=True)
class GamDesign:
    matrix: np.ndarray
    factor_indices: tuple[int, ...]
    month_index: int | None
    year_index: int | None
    smooth_indices: tuple[int, ...]


def gam_design(data: pd.DataFrame, columns: list[str]) -> GamDesign:
    arrays: list[np.ndarray] = []
    factors: list[int] = []
    smooths: list[int] = []
    month_index: int | None = None
    year_index: int | None = None

    def add(values: np.ndarray) -> int:
        arrays.append(np.asarray(values, dtype=float).reshape(-1))
        return len(arrays) - 1

    for name in ("Pond", "Layer", "Season_verified"):
        if name in columns:
            categories = sorted(data[name].astype(str).unique())
            mapping = {value: index for index, value in enumerate(categories)}
            factors.append(add(data[name].astype(str).map(mapping).to_numpy()))
    if "Month_sin" in columns:
        month_index = add(data["Month"].to_numpy(dtype=float))
    if "Year" in columns:
        year_index = add(data["Year"].to_numpy(dtype=float))
    for name in ("Temperature", "Salinity", "pH", "DO", "Turbidity"):
        if name in columns:
            smooths.append(add(data[name].to_numpy(dtype=float)))
    for name in (
        "Heatwave_Tmax_q90_3day",
        "Typhoon",
        "ColdSurge",
        "LowPressure",
    ):
        if name in columns:
            factors.append(add(data[name].to_numpy(dtype=float)))
    return GamDesign(
        matrix=np.column_stack(arrays),
        factor_indices=tuple(factors),
        month_index=month_index,
        year_index=year_index,
        smooth_indices=tuple(smooths),
    )


def gam_predictions(
    data: pd.DataFrame,
    columns: list[str],
    splits: list[tuple[np.ndarray, np.ndarray]],
) -> np.ndarray:
    try:
        from pygam import LinearGAM, f, l, s
    except ImportError as error:
        raise RuntimeError(
            "pyGAM is required for the generalized additive model"
        ) from error

    design = gam_design(data, columns)
    terms = None
    for index in design.factor_indices:
        term = f(index)
        terms = term if terms is None else terms + term
    if design.month_index is not None:
        # Half-month edge knots make December and January adjacent while
        # preserving integer month centres from 1 through 12.
        term = s(
            design.month_index,
            n_splines=6,
            basis="cp",
            edge_knots=[0.5, 12.5],
        )
        terms = term if terms is None else terms + term
    if design.year_index is not None:
        term = l(design.year_index)
        terms = term if terms is None else terms + term
    for index in design.smooth_indices:
        term = s(index, n_splines=5)
        terms = term if terms is None else terms + term
    if terms is None:
        raise ValueError("GAM design has no terms")

    y = data[["PCoA1", "PCoA2"]].to_numpy(dtype=float)
    predictions = np.full_like(y, np.nan)
    for train, test in splits:
        for axis in range(2):
            model = LinearGAM(terms=terms, lam=0.6)
            model.fit(design.matrix[train], y[train, axis])
            predictions[test, axis] = model.predict(design.matrix[test])
    return predictions


def evaluate(
    observed: np.ndarray, predicted: np.ndarray, algorithm: str, predictor_set: str
) -> dict[str, float | str]:
    r2 = [r2_score(observed[:, axis], predicted[:, axis]) for axis in range(2)]
    mae = [
        mean_absolute_error(observed[:, axis], predicted[:, axis])
        for axis in range(2)
    ]
    return {
        "Algorithm": algorithm,
        "Predictor_set": predictor_set,
        "PCoA1_R2": r2[0],
        "PCoA2_R2": r2[1],
        "Mean_R2": float(np.mean(r2)),
        "PCoA1_MAE": mae[0],
        "PCoA2_MAE": mae[1],
        "Mean_MAE": float(np.mean(mae)),
    }


def run_cohort(data: pd.DataFrame, cohort: str, algorithms: list[str]) -> pd.DataFrame:
    if cohort == "all_ponds":
        work = data.copy()
        expected = (160, 89)
    elif cohort == "grow_out_only":
        work = data.loc[data["Pond"].isin(["T2", "T3", "T4"])].copy()
        expected = (122, 70)
    else:
        raise ValueError(cohort)
    work = work.reset_index(drop=True)
    observed = work[["PCoA1", "PCoA2"]].to_numpy(dtype=float)
    splitter = GroupKFold(n_splits=N_FOLDS)
    split_input = np.zeros((len(work), 1))
    splits = list(splitter.split(split_input, groups=work["PondMonth"]))
    if (len(work), work["PondMonth"].nunique()) != expected:
        raise ValueError(
            f"Unexpected cohort size for {cohort}: "
            f"{len(work)} samples and {work['PondMonth'].nunique()} groups"
        )

    rows = [
        evaluate(
            observed,
            mean_only_predictions(work, splits),
            "Mean-only baseline",
            "None",
        )
    ]

    factories: dict[str, tuple[str, Callable[[], object], bool, bool]] = {
        "ridge": ("Ridge", lambda: Ridge(alpha=1.0), True, False),
        "random_forest": (
            "Random forest",
            lambda: RandomForestRegressor(
                n_estimators=500,
                min_samples_leaf=3,
                max_features=0.70,
                random_state=SEED,
                n_jobs=-1,
            ),
            False,
            False,
        ),
        "gradient_boosting": (
            "Gradient boosting",
            lambda: GradientBoostingRegressor(
                n_estimators=200,
                learning_rate=0.05,
                max_depth=2,
                min_samples_leaf=3,
                random_state=SEED,
            ),
            False,
            False,
        ),
        "support_vector": (
            "Support vector regression",
            lambda: SVR(kernel="rbf", C=1.0, epsilon=0.1, gamma="scale"),
            True,
            True,
        ),
    }

    if "xgboost" in algorithms:
        try:
            from xgboost import XGBRegressor
        except ImportError as error:
            raise RuntimeError("xgboost is required for the XGBoost model") from error
        factories["xgboost"] = (
            "XGBoost",
            lambda: XGBRegressor(
                n_estimators=300,
                learning_rate=0.05,
                max_depth=2,
                min_child_weight=3,
                reg_lambda=1.0,
                reg_alpha=0.0,
                gamma=0.0,
                subsample=1.0,
                colsample_bytree=1.0,
                objective="reg:squarederror",
                random_state=SEED,
                n_jobs=1,
                verbosity=0,
            ),
            False,
            False,
        )

    for algorithm in algorithms:
        for predictor_set, columns in FEATURE_SETS.items():
            if algorithm == "gam":
                predictions = gam_predictions(work, columns, splits)
                label = "Generalized additive model"
            else:
                label, factory, scale_x, scale_y = factories[algorithm]
                predictions = standard_model_predictions(
                    work,
                    columns,
                    splits,
                    factory,
                    scale_x=scale_x,
                    scale_y=scale_y,
                )
            rows.append(evaluate(observed, predictions, label, predictor_set))

    result = pd.DataFrame(rows)
    result["Delta_mean_R2"] = np.nan
    for algorithm, indices in result.groupby("Algorithm", sort=False).groups.items():
        if algorithm == "Mean-only baseline":
            continue
        positions = list(indices)
        result.loc[positions, "Delta_mean_R2"] = result.loc[
            positions, "Mean_R2"
        ].diff()
    result.insert(0, "Cohort", cohort)
    result.insert(1, "Samples", len(work))
    result.insert(2, "Pond_month_groups", work["PondMonth"].nunique())
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--asv", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--algorithms",
        nargs="+",
        choices=[
            "ridge",
            "random_forest",
            "gradient_boosting",
            "xgboost",
            "support_vector",
            "gam",
        ],
        default=[
            "ridge",
            "random_forest",
            "gradient_boosting",
            "xgboost",
            "support_vector",
            "gam",
        ],
    )
    args = parser.parse_args()

    outdir = args.output_dir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    data = load_analysis_data(args.asv, args.metadata)
    all_ponds = run_cohort(data, "all_ponds", args.algorithms)
    grow_out = run_cohort(data, "grow_out_only", args.algorithms)
    all_ponds.to_csv(outdir / "multialgorithm_all_ponds.csv", index=False)
    grow_out.to_csv(outdir / "multialgorithm_grow_out_only.csv", index=False)
    pd.concat([all_ponds, grow_out], ignore_index=True).to_csv(
        outdir / "multialgorithm_both_cohorts.csv", index=False
    )
    print(pd.concat([all_ponds, grow_out], ignore_index=True).to_string(index=False))


if __name__ == "__main__":
    main()
