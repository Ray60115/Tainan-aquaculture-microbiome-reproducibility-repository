#!/usr/bin/env python3
"""Reproduce grouped random-forest reconstruction and permutation importance."""

from __future__ import annotations

from pathlib import Path
import hashlib
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-tainan")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import GroupKFold

from reanalyze_dependence_aware import pcoa, relative_abundance_rows


CODE = Path(__file__).resolve().parent
ROOT = CODE.parent
DATA = ROOT / "data" / "primary"
OUT = Path(os.environ.get("TAINAN_RESULTS_DIR", ROOT / "reproduced_results" / "taxonomy"))
TABLES = OUT / "tables"
FIGURES = OUT / "figures"

SEED = 42
N_TREES = 500
MIN_LEAF = 3
MAX_FEATURES = 0.70
N_FOLDS = 5
N_IMPORTANCE_REPEATS = 50


def independent_rng(*parts: object) -> np.random.Generator:
    key = "|".join(map(str, parts))
    seed = int.from_bytes(hashlib.sha256(key.encode()).digest()[:8], "little")
    return np.random.default_rng(seed)


def load_analysis_data() -> pd.DataFrame:
    meta = pd.read_csv(TABLES / "validated_sample_metadata_with_shannon.csv")
    asv = pd.read_csv(DATA / "asv_table_168_profiles.csv.gz")
    counts = asv[meta["SampleID"].tolist()].T.to_numpy(dtype=float)
    relative = relative_abundance_rows(counts)
    distance = squareform(pdist(relative, metric="braycurtis"))
    scores, _, _ = pcoa(distance)
    meta["PCoA1"] = scores[:, 0]
    meta["PCoA2"] = scores[:, 1]
    work = meta.loc[meta["DisturbanceComplete"]].copy().reset_index(drop=True)
    month = work["Date"].astype(str).str[4:].astype(int)
    work["Month_sin"] = np.sin(2 * np.pi * month / 12)
    work["Month_cos"] = np.cos(2 * np.pi * month / 12)
    work["Year"] = work["Date"].astype(str).str[:4].astype(int)
    return work


def encoded_features(data: pd.DataFrame, raw_columns: list[str]) -> pd.DataFrame:
    categorical = [x for x in ["Pond", "Layer", "Season_verified"] if x in raw_columns]
    return pd.get_dummies(
        data[raw_columns], columns=categorical, drop_first=False, dtype=float
    )


def fit_oof(
    data: pd.DataFrame, feature_columns: list[str]
) -> tuple[pd.DataFrame, list[tuple[np.ndarray, np.ndarray, list[RandomForestRegressor]]], pd.DataFrame]:
    x = encoded_features(data, feature_columns)
    y = data[["PCoA1", "PCoA2"]].to_numpy()
    predictions = np.full_like(y, np.nan)
    fitted = []
    splitter = GroupKFold(n_splits=N_FOLDS)
    for fold, (train, test) in enumerate(splitter.split(x, groups=data["PondMonth"]), start=1):
        models = []
        for axis in range(2):
            model = RandomForestRegressor(
                n_estimators=N_TREES,
                min_samples_leaf=MIN_LEAF,
                max_features=MAX_FEATURES,
                random_state=SEED,
                n_jobs=-1,
            )
            model.fit(x.iloc[train], y[train, axis])
            predictions[test, axis] = model.predict(x.iloc[test])
            models.append(model)
        fitted.append((train, test, models))
    pred = data[["SampleID", "Pond", "Date", "PondMonth", "PCoA1", "PCoA2"]].copy()
    pred["Predicted_PCoA1"] = predictions[:, 0]
    pred["Predicted_PCoA2"] = predictions[:, 1]
    return x, fitted, pred


def evaluate(pred: pd.DataFrame, name: str) -> dict[str, float | str]:
    r2_1 = r2_score(pred["PCoA1"], pred["Predicted_PCoA1"])
    r2_2 = r2_score(pred["PCoA2"], pred["Predicted_PCoA2"])
    mae_1 = mean_absolute_error(pred["PCoA1"], pred["Predicted_PCoA1"])
    mae_2 = mean_absolute_error(pred["PCoA2"], pred["Predicted_PCoA2"])
    return {
        "Model": name,
        "PCoA1_R2": r2_1,
        "PCoA2_R2": r2_2,
        "Mean_R2": np.mean([r2_1, r2_2]),
        "PCoA1_MAE": mae_1,
        "PCoA2_MAE": mae_2,
        "Mean_MAE": np.mean([mae_1, mae_2]),
    }


def grouped_importance(
    data: pd.DataFrame,
    x: pd.DataFrame,
    fitted: list[tuple[np.ndarray, np.ndarray, list[RandomForestRegressor]]],
    baseline_pred: pd.DataFrame,
) -> pd.DataFrame:
    # Keep the original seed keys separate from publication-facing labels.
    # The initial implementation derived each permutation stream from the
    # group label; changing only "Temperature" to "Water temperature" would
    # otherwise change the numerical importance result.
    groups = {
        "Pond identity": ([c for c in x if c.startswith("Pond_")], "Pond identity"),
        "Seasonal timing": (
            [c for c in x if c.startswith("Season_verified_")] + ["Month_sin", "Month_cos"],
            "Seasonal timing",
        ),
        "Year": (["Year"], "Year"),
        "Water temperature": (["Temperature"], "Temperature"),
        "Salinity": (["Salinity"], "Salinity-related variable"),
        "pH": (["pH"], "pH"),
        "Dissolved oxygen": (["DO"], "Dissolved oxygen"),
        "Turbidity": (["Turbidity"], "Turbidity"),
        "Heat event": (["Heatwave_Tmax_q90_3day"], "Heat event"),
        "Typhoon": (["Typhoon"], "Typhoon"),
        "Layer": ([c for c in x if c.startswith("Layer_")], "Layer"),
        "Cold surge": (["ColdSurge"], "Cold surge"),
        "Low pressure": (["LowPressure"], "Low pressure"),
    }
    y = data[["PCoA1", "PCoA2"]].to_numpy()
    baseline = np.array([
        mean_absolute_error(y[:, 0], baseline_pred["Predicted_PCoA1"]),
        mean_absolute_error(y[:, 1], baseline_pred["Predicted_PCoA2"]),
    ])
    rows = []
    for group, (columns, seed_key) in groups.items():
        repeated_predictions = np.full(
            (N_IMPORTANCE_REPEATS, len(y), 2), np.nan, dtype=float
        )
        column_positions = [x.columns.get_loc(column) for column in columns]
        for fold, (_, test, models) in enumerate(fitted):
            base = x.iloc[test].to_numpy(dtype=float)
            stacked = np.tile(base, (N_IMPORTANCE_REPEATS, 1))
            for repeat in range(N_IMPORTANCE_REPEATS):
                start = repeat * len(test)
                stop = start + len(test)
                rng = independent_rng(SEED, seed_key, repeat, fold)
                order = rng.permutation(len(test))
                stacked[start:stop, column_positions] = base[
                    order[:, None], column_positions
                ]
            stacked_frame = pd.DataFrame(stacked, columns=x.columns)
            for axis, model in enumerate(models):
                predicted = model.predict(stacked_frame).reshape(
                    N_IMPORTANCE_REPEATS, len(test)
                )
                repeated_predictions[:, test, axis] = predicted
        repeated_mae = np.mean(
            np.abs(y[None, :, :] - repeated_predictions), axis=1
        )
        values = 100 * (repeated_mae / baseline[None, :] - 1)
        overall = values.mean(axis=1)
        rows.append({
            "Predictor_group": group,
            "Mean_increase_MAE_percent": overall.mean(),
            "SD_across_repeats_percent": overall.std(ddof=1),
            "PCoA1_mean_increase_MAE_percent": values[:, 0].mean(),
            "PCoA2_mean_increase_MAE_percent": values[:, 1].mean(),
            "Repeats": N_IMPORTANCE_REPEATS,
        })
    return pd.DataFrame(rows).sort_values("Mean_increase_MAE_percent", ascending=False)


def make_figures(performance: pd.DataFrame, predictions: pd.DataFrame, importance: pd.DataFrame) -> None:
    colors = {"PCoA1": "#4C78A8", "PCoA2": "#F58518"}
    fig = plt.figure(figsize=(11.5, 9.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=0.34, wspace=0.24)
    top = fig.add_subplot(gs[0, :])
    lower_axes = [fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
    x = np.arange(len(performance))
    width = 0.34
    bars1 = top.bar(
        x - width / 2, performance["PCoA1_R2"], width,
        label="PCoA1", color=colors["PCoA1"],
    )
    bars2 = top.bar(
        x + width / 2, performance["PCoA2_R2"], width,
        label="PCoA2", color=colors["PCoA2"],
    )
    top.bar_label(bars1, fmt="%.3f", padding=3, fontsize=8.2)
    top.bar_label(bars2, fmt="%.3f", padding=3, fontsize=8.2)
    top.set_xticks(
        x,
        ["Structure", "+ Season", "+ Water\nquality", "+ Disturbances"],
    )
    top.set_ylabel("Grouped cross-validated R²")
    top.set_title("(a) Cross-validated reconstruction of global PCoA scores", loc="left", fontweight="bold")
    top.legend(frameon=False, ncol=2)
    top.axhline(0, color="#999999", lw=0.7)
    top.set_ylim(-0.03, max(performance[["PCoA1_R2", "PCoA2_R2"]].max()) + 0.12)
    top.spines["top"].set_visible(False)
    top.spines["right"].set_visible(False)
    for axis, ax2 in enumerate(lower_axes):
        obs = predictions[f"PCoA{axis + 1}"]
        pred = predictions[f"Predicted_PCoA{axis + 1}"]
        for pond, sub in predictions.assign(obs=obs, pred=pred).groupby("Pond"):
            ax2.scatter(sub["obs"], sub["pred"], s=26, alpha=0.72, label=pond)
        limits = [min(obs.min(), pred.min()), max(obs.max(), pred.max())]
        ax2.plot(limits, limits, ls="--", color="#777777", lw=1)
        ax2.set_xlabel(f"Observed PCoA{axis + 1}")
        ax2.set_ylabel(f"Out-of-fold reconstructed PCoA{axis + 1}")
        ax2.set_title(f"({'bc'[axis]}) PCoA{axis + 1}", loc="left", fontweight="bold")
        ax2.text(0.04, 0.94,
                 f"R² = {r2_score(obs, pred):.3f}\nMAE = {mean_absolute_error(obs, pred):.3f}",
                 transform=ax2.transAxes, va="top")
        ax2.spines["top"].set_visible(False)
        ax2.spines["right"].set_visible(False)
    lower_axes[1].legend(title="Pond", frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURES / "Fig11_panels_a_to_c_reconstruction.png", dpi=400, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.3, 6.5))
    plot = importance.sort_values("Mean_increase_MAE_percent")
    ax.barh(plot["Predictor_group"], plot["Mean_increase_MAE_percent"],
            xerr=plot["SD_across_repeats_percent"], color="#4C78A8", alpha=0.9,
            capsize=2)
    for y, (value, error) in enumerate(
        zip(
            plot["Mean_increase_MAE_percent"],
            plot["SD_across_repeats_percent"],
        )
    ):
        ax.text(
            value + error + 0.65,
            y,
            f"{value:.2f}".replace("-", "−"),
            va="center",
            fontsize=8.1,
        )
    ax.axvline(0, color="#777777", lw=0.8)
    ax.set_xlabel("Increase in grouped-CV MAE after permutation (%)")
    ax.set_title("Grouped permutation importance in held-out pond-months", fontweight="bold")
    right_limit = (
        plot["Mean_increase_MAE_percent"]
        + plot["SD_across_repeats_percent"]
    ).max() + 6.0
    left_limit = min(
        -1.0,
        (
            plot["Mean_increase_MAE_percent"]
            - plot["SD_across_repeats_percent"]
        ).min() - 1.0,
    )
    ax.set_xlim(left_limit, right_limit)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIGURES / "Fig11_panel_d_grouped_importance.png", dpi=400, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    data = load_analysis_data()
    feature_sets = {
        "Structure": ["Pond", "Layer"],
        "+ Season": ["Pond", "Layer", "Season_verified", "Month_sin", "Month_cos", "Year"],
        "+ Water quality": [
            "Pond", "Layer", "Season_verified", "Month_sin", "Month_cos", "Year",
            "Temperature", "Salinity", "pH", "DO", "Turbidity",
        ],
        "+ Disturbances": [
            "Pond", "Layer", "Season_verified", "Month_sin", "Month_cos", "Year",
            "Temperature", "Salinity", "pH", "DO", "Turbidity",
            "Heatwave_Tmax_q90_3day", "Typhoon", "ColdSurge", "LowPressure",
        ],
    }
    performance_rows = []
    full = None
    for name, columns in feature_sets.items():
        x, fitted, pred = fit_oof(data, columns)
        performance_rows.append(evaluate(pred, name))
        if name == "+ Disturbances":
            full = (x, fitted, pred)
    performance = pd.DataFrame(performance_rows)
    performance["Incremental_mean_R2"] = performance["Mean_R2"].diff()
    performance.to_csv(TABLES / "rf_cross_validated_performance.csv", index=False)
    assert full is not None
    x, fitted, predictions = full
    predictions.to_csv(TABLES / "rf_oof_predictions.csv", index=False)
    importance = grouped_importance(data, x, fitted, predictions)
    importance.to_csv(TABLES / "rf_grouped_permutation_importance.csv", index=False)
    make_figures(performance, predictions, importance)
    print("Random-forest reconstruction completed.")


if __name__ == "__main__":
    main()
