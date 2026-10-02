#!/usr/bin/env python3
"""
Dependence-aware reanalysis of the Tainan aquaculture microbiome dataset.

This script replaces unrestricted sample-level inference with analyses whose
permutation or pairing units match the sampling design:

* surface and bottom Shannon diversity: paired within pond-month;
* seasonal Shannon diversity and turnover: pond-month/transition observations
  with circular label shifts within pond;
* meteorological comparisons: sampling month as the exposure unit, with
  circular shifts of the complete monthly exposure sequence;
* community season/event PERMANOVA: pond-month-averaged communities, with
  restricted circular permutations;
* layer PERMANOVA: label swaps within complete surface-bottom pairs;
* environmental vector fitting: pond-month ordination, circular permutations
  within pond for water variables and across sampling months for shared
  meteorological variables.

Pond and T1-versus-T2-T4 comparisons are descriptive because the four ponds,
including only one nursery pond, do not provide independent replication of
production category.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
import math
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-tainan")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.stats import wilcoxon


ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "primary"
OUT = Path(os.environ.get("TAINAN_RESULTS_DIR", ROOT / "dependence_aware_results"))
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
for directory in (OUT, TABLES, FIGURES):
    directory.mkdir(parents=True, exist_ok=True)

ASV_FILE = DATA / "asv_table_168_profiles.csv.gz"
META_FILE = DATA / "sample_metadata_168.csv"
WATER_FILE = DATA / "water_quality_168.csv"
CLIMATE_FILE = DATA / "climate_daily_complete_2024.csv"

SEED = 20260728
N_PERM = 9999
RNG = np.random.default_rng(SEED)


def make_rng(key: str) -> np.random.Generator:
    digest = hashlib.sha256(f"{SEED}|{key}".encode("utf-8")).digest()
    seed = int.from_bytes(digest[:8], "little", signed=False)
    return np.random.default_rng(seed)


def bh_adjust(p_values: np.ndarray | list[float]) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    q = np.full_like(p, np.nan)
    valid = np.isfinite(p)
    pv = p[valid]
    if len(pv) == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    q[valid] = restored
    return q


def shannon_from_counts(counts: np.ndarray) -> float:
    x = np.asarray(counts, dtype=float)
    total = x.sum()
    if total <= 0:
        return np.nan
    p = x[x > 0] / total
    return float(-(p * np.log(p)).sum())


def relative_abundance_rows(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    totals = x.sum(axis=1, keepdims=True)
    if np.any(totals <= 0):
        raise ValueError("At least one microbial sample has a zero total.")
    return x / totals


def pcoa(distance_matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    d = np.asarray(distance_matrix, dtype=float)
    n = len(d)
    j = np.eye(n) - np.ones((n, n)) / n
    g = -0.5 * j @ (d ** 2) @ j
    eigvals, eigvecs = np.linalg.eigh(g)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    positive = eigvals > max(1e-12, eigvals[0] * 1e-10)
    vals = eigvals[positive]
    coords = eigvecs[:, positive] * np.sqrt(vals)
    explained = vals / vals.sum() * 100
    return coords, vals, explained


def design_matrix(df: pd.DataFrame, categorical=(), continuous=()) -> np.ndarray:
    parts = [np.ones((len(df), 1), dtype=float)]
    for col in categorical:
        dummies = pd.get_dummies(df[col].astype("category"), drop_first=True, dtype=float)
        if dummies.shape[1]:
            parts.append(dummies.to_numpy(dtype=float))
    for col in continuous:
        z = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
        sd = np.nanstd(z, ddof=0)
        z = (z - np.nanmean(z)) / sd if sd > 0 else np.zeros_like(z)
        parts.append(z[:, None])
    return np.column_stack(parts)


def partial_f_univariate(y: np.ndarray, x_reduced: np.ndarray, x_full: np.ndarray) -> tuple[float, float]:
    y = np.asarray(y, dtype=float)
    br = np.linalg.lstsq(x_reduced, y, rcond=None)[0]
    bf = np.linalg.lstsq(x_full, y, rcond=None)[0]
    sse_r = float(np.sum((y - x_reduced @ br) ** 2))
    sse_f = float(np.sum((y - x_full @ bf) ** 2))
    df1 = np.linalg.matrix_rank(x_full) - np.linalg.matrix_rank(x_reduced)
    df2 = len(y) - np.linalg.matrix_rank(x_full)
    ss = max(0.0, sse_r - sse_f)
    f = (ss / df1) / (sse_f / df2) if df1 > 0 and df2 > 0 and sse_f > 0 else np.nan
    partial_r2 = ss / sse_r if sse_r > 0 else np.nan
    return float(f), float(partial_r2)


def factor_dummy_matrix(labels: np.ndarray | pd.Series) -> np.ndarray:
    """Return treatment-coded factor columns without rebuilding a DataFrame."""
    values = np.asarray(labels)
    levels = np.unique(values)
    if len(levels) <= 1:
        return np.empty((len(values), 0), dtype=float)
    return np.column_stack([(values == level).astype(float) for level in levels[1:]])


def orthonormal_basis(x: np.ndarray) -> tuple[np.ndarray, int]:
    """Return an orthonormal basis for the column space of x."""
    q, r = np.linalg.qr(np.asarray(x, dtype=float), mode="reduced")
    rank = int(np.linalg.matrix_rank(r))
    return q[:, :rank], rank


def residual_factor_basis(q_reduced: np.ndarray, factor_columns: np.ndarray) -> tuple[np.ndarray, int]:
    """Orthogonalize factor columns against a fixed reduced-model basis."""
    z = np.asarray(factor_columns, dtype=float)
    if z.shape[1] == 0:
        return np.empty((len(z), 0), dtype=float), 0
    residual = z - q_reduced @ (q_reduced.T @ z)
    return orthonormal_basis(residual)


def make_label_permuter(
    df: pd.DataFrame,
    factor: str,
    scheme: str,
    rng: np.random.Generator,
):
    """Precompute permutation groups so each iteration uses only NumPy operations."""
    base = df[factor].to_numpy(copy=True)
    if scheme == "within_pond_circular":
        groups = []
        for _, indices in df.groupby("Pond", sort=False).groups.items():
            idx = np.asarray(list(indices), dtype=int)
            ordered = idx[np.argsort(df.loc[idx, "Date"].to_numpy())]
            groups.append((ordered, base[ordered].copy()))

        def permute():
            out = base.copy()
            for ordered, values in groups:
                shift = int(rng.integers(0, len(values)))
                out[ordered] = np.roll(values, shift)
            return out

        return permute

    if scheme == "date_circular":
        date_factor = df[["Date", factor]].drop_duplicates().sort_values("Date")
        if date_factor.groupby("Date")[factor].nunique(dropna=False).max() != 1:
            raise ValueError(f"{factor} is not unique within Date.")
        dates = date_factor["Date"].to_numpy()
        values = date_factor[factor].to_numpy()
        date_codes = pd.Categorical(df["Date"], categories=dates, ordered=True).codes
        if np.any(date_codes < 0):
            raise ValueError("A sampling date could not be encoded for permutation.")

        def permute():
            shift = int(rng.integers(0, len(values)))
            return np.roll(values, shift)[date_codes]

        return permute

    if scheme == "paired_layer_swap":
        pairs = []
        for _, indices in df.groupby("PondMonth", sort=False).groups.items():
            idx = np.asarray(list(indices), dtype=int)
            if len(idx) != 2:
                raise ValueError("Layer permutation requires exactly two rows per pond-month.")
            pairs.append(idx)
        pair_array = np.asarray(pairs, dtype=int)

        def permute():
            out = base.copy()
            flips = rng.random(len(pair_array)) < 0.5
            selected = pair_array[flips]
            if len(selected):
                out[selected[:, 0]] = base[selected[:, 1]]
                out[selected[:, 1]] = base[selected[:, 0]]
            return out

        return permute

    raise ValueError(scheme)


def gower_matrix(distance_matrix: np.ndarray) -> np.ndarray:
    d = np.asarray(distance_matrix, dtype=float)
    n = len(d)
    j = np.eye(n) - np.ones((n, n)) / n
    return -0.5 * j @ (d ** 2) @ j


def hat(x: np.ndarray) -> np.ndarray:
    return x @ np.linalg.pinv(x)


def model_ss(g: np.ndarray, x: np.ndarray) -> tuple[float, int]:
    """Return trace(HG) without forming the n-by-n hat matrix."""
    q, r = np.linalg.qr(x, mode="reduced")
    rank = int(np.linalg.matrix_rank(r))
    q = q[:, :rank]
    return float(np.trace(q.T @ g @ q)), rank


def partial_f_distance(
    distance_matrix: np.ndarray, x_reduced: np.ndarray, x_full: np.ndarray
) -> tuple[float, float, float]:
    g = gower_matrix(distance_matrix)
    ss_total = float(np.trace(g))
    ss_reduced, rank_reduced = model_ss(g, x_reduced)
    ss_full, rank_full = model_ss(g, x_full)
    ss_factor = ss_full - ss_reduced
    ss_resid = ss_total - ss_full
    df1 = rank_full - rank_reduced
    df2 = len(g) - rank_full
    f = (ss_factor / df1) / (ss_resid / df2) if df1 > 0 and df2 > 0 and ss_resid > 0 else np.nan
    partial_r2 = ss_factor / (ss_factor + ss_resid) if ss_factor + ss_resid > 0 else np.nan
    incremental_total_r2 = ss_factor / ss_total if ss_total > 0 else np.nan
    return float(f), float(partial_r2), float(incremental_total_r2)


def marginal_r2_distance(distance_matrix: np.ndarray, labels: pd.Series) -> float:
    temp = pd.DataFrame({"factor": labels.astype(str).to_numpy()})
    x = design_matrix(temp, categorical=["factor"])
    g = gower_matrix(distance_matrix)
    ss, _ = model_ss(g, x)
    total = float(np.trace(g))
    return ss / total


def circular_shift_labels_within_pond(
    df: pd.DataFrame, factor: str, rng: np.random.Generator
) -> np.ndarray:
    out = df[factor].astype(object).to_numpy(copy=True)
    for _, indices in df.groupby("Pond", sort=False).groups.items():
        idx = np.asarray(list(indices), dtype=int)
        ordered = idx[np.argsort(df.loc[idx, "Date"].to_numpy())]
        values = df.loc[ordered, factor].to_numpy()
        shift = int(rng.integers(0, len(values)))
        out[ordered] = np.roll(values, shift)
    return out


def circular_shift_labels_by_date(
    df: pd.DataFrame, factor: str, rng: np.random.Generator
) -> np.ndarray:
    date_factor = (
        df[["Date", factor]]
        .drop_duplicates()
        .sort_values("Date")
    )
    if date_factor.groupby("Date")[factor].nunique(dropna=False).max() != 1:
        raise ValueError(f"{factor} is not unique within Date.")
    dates = date_factor["Date"].to_numpy()
    values = date_factor[factor].to_numpy()
    shift = int(rng.integers(0, len(values)))
    mapping = dict(zip(dates, np.roll(values, shift)))
    return df["Date"].map(mapping).to_numpy()


def paired_layer_permutation(df: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    labels = df["Layer"].to_numpy(copy=True)
    for _, indices in df.groupby("PondMonth", sort=False).groups.items():
        idx = np.asarray(list(indices), dtype=int)
        if len(idx) != 2:
            raise ValueError("Layer permutation requires exactly two rows per pond-month.")
        if rng.random() < 0.5:
            labels[idx] = labels[idx[::-1]]
    return labels


def restricted_univariate_test(
    df: pd.DataFrame,
    response: str,
    factor: str,
    reduced_categorical: list[str],
    scheme: str,
    n_perm: int = N_PERM,
) -> dict:
    work = df.reset_index(drop=True).copy()
    xr = design_matrix(work, categorical=reduced_categorical)
    q_reduced, rank_reduced = orthonormal_basis(xr)
    y = work[response].to_numpy(dtype=float)
    reduced_residual = y - q_reduced @ (q_reduced.T @ y)
    sse_reduced = float(reduced_residual @ reduced_residual)

    def statistic(labels: np.ndarray) -> tuple[float, float]:
        q_factor, rank_factor = residual_factor_basis(
            q_reduced, factor_dummy_matrix(labels)
        )
        if rank_factor == 0:
            return np.nan, np.nan
        ss_factor = float(np.sum((q_factor.T @ y) ** 2))
        sse_full = max(0.0, sse_reduced - ss_factor)
        df2 = len(y) - rank_reduced - rank_factor
        f_value = (
            (ss_factor / rank_factor) / (sse_full / df2)
            if df2 > 0 and sse_full > 0
            else np.nan
        )
        partial_r2 = ss_factor / sse_reduced if sse_reduced > 0 else np.nan
        return float(f_value), float(partial_r2)

    observed, partial_r2 = statistic(work[factor].to_numpy())
    rng = make_rng(f"univariate|{response}|{factor}|{scheme}|{len(work)}|{n_perm}")
    permute = make_label_permuter(work, factor, scheme, rng)
    exceed = 0
    for _ in range(n_perm):
        f_perm, _ = statistic(permute())
        exceed += int(np.isfinite(f_perm) and f_perm >= observed - 1e-12)
    return {
        "F": observed,
        "partial_R2": partial_r2,
        "p_value": (exceed + 1) / (n_perm + 1),
        "permutations": n_perm,
        "scheme": scheme,
        "n": len(work),
    }


def restricted_permanova(
    df: pd.DataFrame,
    distance_matrix: np.ndarray,
    factor: str,
    reduced_categorical: list[str],
    scheme: str,
    n_perm: int = N_PERM,
) -> dict:
    work = df.reset_index(drop=True).copy()
    xr = design_matrix(work, categorical=reduced_categorical)
    g = gower_matrix(distance_matrix)
    total = float(np.trace(g))
    q_reduced, rank_reduced = orthonormal_basis(xr)
    ss_reduced = float(np.trace(q_reduced.T @ g @ q_reduced))

    def statistic(labels: np.ndarray) -> tuple[float, float, float]:
        q_factor, rank_factor = residual_factor_basis(
            q_reduced, factor_dummy_matrix(labels)
        )
        if rank_factor == 0:
            return np.nan, np.nan, np.nan
        ss_factor = float(np.trace(q_factor.T @ g @ q_factor))
        ss_resid = total - ss_reduced - ss_factor
        df2 = len(work) - rank_reduced - rank_factor
        f_value = (
            (ss_factor / rank_factor) / (ss_resid / df2)
            if df2 > 0 and ss_resid > 0
            else np.nan
        )
        partial_r2 = (
            ss_factor / (ss_factor + ss_resid)
            if ss_factor + ss_resid > 0
            else np.nan
        )
        return float(f_value), float(partial_r2), float(ss_factor / total)

    observed, partial_r2, incremental = statistic(work[factor].to_numpy())
    marginal = marginal_r2_distance(distance_matrix, work[factor])
    rng = make_rng(f"permanova|{factor}|{scheme}|{len(work)}|{n_perm}|{'-'.join(reduced_categorical)}")
    permute = make_label_permuter(work, factor, scheme, rng)
    exceed = 0
    for _ in range(n_perm):
        f_perm, _, _ = statistic(permute())
        exceed += int(np.isfinite(f_perm) and f_perm >= observed - 1e-12)
    return {
        "F": observed,
        "marginal_R2": marginal,
        "partial_R2": partial_r2,
        "incremental_total_R2": incremental,
        "p_value": (exceed + 1) / (n_perm + 1),
        "permutations": n_perm,
        "scheme": scheme,
        "n": len(work),
    }


def vector_fit(coords2: np.ndarray, values: np.ndarray) -> tuple[float, float, float]:
    z = np.asarray(values, dtype=float)
    r1 = float(np.corrcoef(z, coords2[:, 0])[0, 1])
    r2 = float(np.corrcoef(z, coords2[:, 1])[0, 1])
    x = np.column_stack([np.ones(len(z)), coords2])
    fitted = x @ np.linalg.lstsq(x, z, rcond=None)[0]
    sst = float(np.sum((z - z.mean()) ** 2))
    ssr = float(np.sum((fitted - z.mean()) ** 2))
    fit_r2 = ssr / sst if sst > 0 else np.nan
    return r1, r2, fit_r2


def make_vector_statistic(coords2: np.ndarray):
    """Precompute fixed ordination terms for repeated vector-fit permutations."""
    centered_coords = np.asarray(coords2, dtype=float) - np.mean(coords2, axis=0)
    coord_ss = np.sum(centered_coords ** 2, axis=0)
    q_coords, _ = orthonormal_basis(centered_coords)

    def statistic(values: np.ndarray) -> tuple[float, float, float]:
        centered_values = np.asarray(values, dtype=float)
        centered_values = centered_values - centered_values.mean()
        value_ss = float(centered_values @ centered_values)
        correlations = (
            centered_coords.T @ centered_values
            / np.sqrt(coord_ss * value_ss)
        )
        fit_r2 = (
            float(np.sum((q_coords.T @ centered_values) ** 2) / value_ss)
            if value_ss > 0
            else np.nan
        )
        return float(correlations[0]), float(correlations[1]), fit_r2

    return statistic


def dispersion_f(coords: np.ndarray, labels: np.ndarray) -> float:
    labels = np.asarray(labels)
    distances = np.empty(len(labels), dtype=float)
    for group in np.unique(labels):
        idx = np.where(labels == group)[0]
        centroid = coords[idx].mean(axis=0)
        distances[idx] = np.sqrt(np.sum((coords[idx] - centroid) ** 2, axis=1))
    grand = distances.mean()
    groups = np.unique(labels)
    ss_between = sum(
        np.sum(labels == group) * (distances[labels == group].mean() - grand) ** 2
        for group in groups
    )
    ss_within = sum(
        np.sum((distances[labels == group] - distances[labels == group].mean()) ** 2)
        for group in groups
    )
    df_between = len(groups) - 1
    df_within = len(labels) - len(groups)
    return float((ss_between / df_between) / (ss_within / df_within))


def restricted_dispersion_test(
    df: pd.DataFrame,
    distance_matrix: np.ndarray,
    factor: str,
    scheme: str,
    n_perm: int = N_PERM,
) -> dict:
    work = df.reset_index(drop=True).copy()
    coords, _, _ = pcoa(distance_matrix)
    observed = dispersion_f(coords, work[factor].astype(str).to_numpy())
    rng = make_rng(f"dispersion|{factor}|{scheme}|{len(work)}|{n_perm}")
    permute = make_label_permuter(work, factor, scheme, rng)
    exceed = 0
    for _ in range(n_perm):
        fp = dispersion_f(coords, np.asarray(permute(), dtype=str))
        exceed += int(fp >= observed - 1e-12)
    return {
        "F": observed,
        "p_value": (exceed + 1) / (n_perm + 1),
        "n": len(work),
        "permutations": n_perm,
        "scheme": scheme,
    }


def permute_continuous_within_pond(
    df: pd.DataFrame, variable: str, rng: np.random.Generator
) -> np.ndarray:
    out = df[variable].to_numpy(copy=True)
    for _, indices in df.groupby("Pond", sort=False).groups.items():
        idx = np.asarray(list(indices), dtype=int)
        ordered = idx[np.argsort(df.loc[idx, "Date"].to_numpy())]
        values = df.loc[ordered, variable].to_numpy()
        shift = int(rng.integers(0, len(values)))
        out[ordered] = np.roll(values, shift)
    return out


def permute_continuous_by_date(
    df: pd.DataFrame, variable: str, rng: np.random.Generator
) -> np.ndarray:
    date_value = df[["Date", variable]].drop_duplicates().sort_values("Date")
    if date_value.groupby("Date")[variable].nunique(dropna=False).max() != 1:
        raise ValueError(f"{variable} is not unique within Date.")
    dates = date_value["Date"].to_numpy()
    values = date_value[variable].to_numpy()
    shift = int(rng.integers(0, len(values)))
    mapping = dict(zip(dates, np.roll(values, shift)))
    return df["Date"].map(mapping).to_numpy(dtype=float)


def make_continuous_permuter(
    df: pd.DataFrame,
    variable: str,
    scheme: str,
    rng: np.random.Generator,
):
    """Precompute indices for circular permutation of a continuous variable."""
    base = df[variable].to_numpy(dtype=float, copy=True)
    if scheme == "within_pond_circular":
        groups = []
        for _, indices in df.groupby("Pond", sort=False).groups.items():
            idx = np.asarray(list(indices), dtype=int)
            ordered = idx[np.argsort(df.loc[idx, "Date"].to_numpy())]
            groups.append((ordered, base[ordered].copy()))

        def permute():
            out = base.copy()
            for ordered, values in groups:
                shift = int(rng.integers(0, len(values)))
                out[ordered] = np.roll(values, shift)
            return out

        return permute

    if scheme == "date_circular":
        date_value = df[["Date", variable]].drop_duplicates().sort_values("Date")
        if date_value.groupby("Date")[variable].nunique(dropna=False).max() != 1:
            raise ValueError(f"{variable} is not unique within Date.")
        dates = date_value["Date"].to_numpy()
        values = date_value[variable].to_numpy(dtype=float)
        date_codes = pd.Categorical(df["Date"], categories=dates, ordered=True).codes
        if np.any(date_codes < 0):
            raise ValueError("A sampling date could not be encoded for permutation.")

        def permute():
            shift = int(rng.integers(0, len(values)))
            return np.roll(values, shift)[date_codes]

        return permute

    raise ValueError(scheme)


def blockwise_r2(
    df: pd.DataFrame,
    distance_matrix: np.ndarray,
    blocks: dict[str, tuple[list[str], list[str]]],
) -> tuple[float, dict[str, float]]:
    categorical = []
    continuous = []
    for cat, cont in blocks.values():
        categorical.extend(cat)
        continuous.extend(cont)
    x_full = design_matrix(df, categorical=categorical, continuous=continuous)
    g = gower_matrix(distance_matrix)
    total = float(np.trace(g))
    full_ss, _ = model_ss(g, x_full)
    full_r2 = float(full_ss / total)
    unique = {}
    for name in blocks:
        cats = []
        conts = []
        for other, (cat, cont) in blocks.items():
            if other == name:
                continue
            cats.extend(cat)
            conts.extend(cont)
        x_red = design_matrix(df, categorical=cats, continuous=conts)
        reduced_ss, _ = model_ss(g, x_red)
        reduced_r2 = float(reduced_ss / total)
        unique[name] = full_r2 - reduced_r2
    return full_r2, unique


def validate_and_load():
    meta = pd.read_csv(META_FILE)
    water = pd.read_csv(WATER_FILE)
    climate = pd.read_csv(CLIMATE_FILE)
    asv = pd.read_csv(ASV_FILE)

    if meta["SampleID"].duplicated().any():
        raise ValueError("Duplicate SampleID values in metadata.")
    key = ["Pond", "Date", "Layer"]
    check_m = meta[key + ["Temperature", "Salinity", "pH", "DO", "Turbidity"]].sort_values(key).reset_index(drop=True)
    check_w = water.sort_values(key).reset_index(drop=True)
    if not check_m[key].equals(check_w[key]):
        raise ValueError("Water and metadata keys do not match.")
    numeric = ["Temperature", "Salinity", "pH", "DO", "Turbidity"]
    if not np.allclose(check_m[numeric], check_w[numeric], equal_nan=True):
        raise ValueError("Water values differ between Table_A and metadata.")

    sample_cols = [c for c in asv.columns if c in set(meta["SampleID"])]
    if len(sample_cols) != len(meta):
        missing = sorted(set(meta["SampleID"]) - set(sample_cols))
        raise ValueError(f"ASV table is missing metadata samples: {missing}")
    counts = asv[sample_cols].T.to_numpy(dtype=float)
    rel = relative_abundance_rows(counts)
    rel_df = pd.DataFrame(rel, index=sample_cols, columns=asv["id"].astype(str))

    meta = meta.set_index("SampleID").loc[sample_cols].reset_index()
    meta["PondMonth"] = meta["Pond"].astype(str) + "-" + meta["Date"].astype(str)
    meta["Date_dt"] = pd.to_datetime(meta["Date"].astype(str), format="%Y%m")

    # The season field has already been reconciled against the authoritative
    # sampling-occasion file.  A missing value is intentional when the ASV
    # profile averages occasions on opposite sides of a solar-term boundary.
    # These rows remain available to analyses that do not use season.
    if "Season_verified" not in meta.columns:
        meta["Season_verified"] = meta["Season"]
    allowed = {"Spring", "Summer", "Autumn", "Winter"}
    observed = set(meta["Season_verified"].dropna().astype(str))
    if not observed.issubset(allowed):
        raise ValueError(f"Unexpected season labels: {sorted(observed - allowed)}")

    unique_pressure = meta[["Date", "MinPressure"]].drop_duplicates().dropna()
    low_threshold = float(unique_pressure["MinPressure"].quantile(0.10))
    meta["LowPressure"] = np.where(
        meta["MinPressure"].notna(),
        (meta["MinPressure"] <= low_threshold).astype(float),
        np.nan,
    )
    complete_cols = ["Heatwave_Tmax_q90_3day", "Typhoon", "ColdSurge", "LowPressure"]
    meta["DisturbanceComplete"] = meta[complete_cols].notna().all(axis=1)

    shannon_values = np.apply_along_axis(shannon_from_counts, 1, counts)
    meta["Shannon"] = shannon_values
    return meta, asv, rel_df, climate, low_threshold


def aggregate_pond_month(meta: pd.DataFrame, rel_df: pd.DataFrame):
    meta_indexed = meta.set_index("SampleID")
    rel_with_key = rel_df.join(meta_indexed[["PondMonth"]])
    pm_rel = rel_with_key.groupby("PondMonth").mean(numeric_only=True)

    numeric_mean = [
        "Temperature", "Salinity", "pH", "DO", "Turbidity", "Shannon",
        "MeanAirTemp", "Tmax", "Tmin", "TotalRainfall", "MaxDailyRainfall",
        "RainDays", "MeanWindSpeed", "MaxGust", "MeanPressure", "MinPressure",
        "Typhoon", "ColdSurge", "Heatwave_Tmax_q90_3day", "LowPressure",
    ]
    pm = meta.groupby("PondMonth", as_index=False).agg(
        Pond=("Pond", "first"),
        Date=("Date", "first"),
        Date_dt=("Date_dt", "first"),
        Season=("Season_verified", "first"),
        n_layers=("Layer", "nunique"),
        **{c: (c, "mean") for c in numeric_mean},
    )
    pm = pm.set_index("PondMonth").loc[pm_rel.index].reset_index()
    return pm, pm_rel.loc[pm["PondMonth"]]


def run_shannon(meta: pd.DataFrame, pm: pd.DataFrame):
    descriptives = []
    for grouping, frame, col in [
        ("Pond", meta, "Pond"),
        ("Layer", meta, "Layer"),
        ("Production category", meta.assign(Production=np.where(meta["Pond"].eq("T1"), "Nursery (T1)", "Grow-out (T2–T4)")), "Production"),
        ("Season, grow-out pond-month", pm[pm["Pond"].ne("T1")], "Season"),
    ]:
        for group, sub in frame.groupby(col, dropna=False):
            x = sub["Shannon"].dropna()
            descriptives.append({
                "Grouping": grouping,
                "Category": group,
                "n": len(x),
                "mean": x.mean(),
                "sd": x.std(ddof=1),
                "median": x.median(),
                "q1": x.quantile(0.25),
                "q3": x.quantile(0.75),
                "minimum": x.min(),
                "maximum": x.max(),
            })
    pd.DataFrame(descriptives).to_csv(TABLES / "shannon_descriptive_summary.csv", index=False)

    paired = (
        meta.pivot_table(index=["Pond", "Date", "PondMonth"], columns="Layer", values="Shannon", aggfunc="first")
        .dropna(subset=["D", "U"])
        .reset_index()
    )
    w = wilcoxon(paired["D"], paired["U"], alternative="two-sided", method="auto")
    tests = [{
        "Analysis": "Sampling layer",
        "Unit": "Complete surface-bottom pond-month pair",
        "n": len(paired),
        "Statistic": float(w.statistic),
        "Effect": float(np.median(paired["D"] - paired["U"])),
        "Effect_definition": "Median paired difference (bottom - surface)",
        "p_value": float(w.pvalue),
        "q_value": np.nan,
        "Method": "Two-sided Wilcoxon signed-rank test",
        "Permutation_scheme": "Not applicable; explicit pairing",
    }]

    grow = pm[pm["Pond"].ne("T1") & pm["Season"].notna()].sort_values(["Pond", "Date"]).reset_index(drop=True)
    season_result = restricted_univariate_test(
        grow, "Shannon", "Season", ["Pond"], "within_pond_circular"
    )
    tests.append({
        "Analysis": "Season within grow-out ponds",
        "Unit": "Pond-month mean",
        "n": season_result["n"],
        "Statistic": season_result["F"],
        "Effect": season_result["partial_R2"],
        "Effect_definition": "Partial R2 after pond identity",
        "p_value": season_result["p_value"],
        "q_value": np.nan,
        "Method": "Partial F test",
        "Permutation_scheme": "9,999 circular shifts of season labels within pond",
    })

    complete = pm[
        pm[["Typhoon", "ColdSurge", "Heatwave_Tmax_q90_3day", "LowPressure"]]
        .notna().all(axis=1)
    ].copy()
    # Pond-center first, then form one system-level value for each sampling month.
    complete["Shannon_centered"] = complete["Shannon"] - complete.groupby("Pond")["Shannon"].transform("mean")
    monthly = complete.groupby("Date", as_index=False).agg(
        Shannon_centered=("Shannon_centered", "mean"),
        Typhoon=("Typhoon", "first"),
        ColdSurge=("ColdSurge", "first"),
        Heatwave_Tmax_q90_3day=("Heatwave_Tmax_q90_3day", "first"),
        LowPressure=("LowPressure", "first"),
        n_ponds=("Pond", "nunique"),
    ).sort_values("Date").reset_index(drop=True)
    event_rows = []
    for factor in ["Typhoon", "ColdSurge", "Heatwave_Tmax_q90_3day", "LowPressure"]:
        x = monthly["Shannon_centered"].to_numpy()
        labels = monthly[factor].to_numpy(dtype=int)
        observed = float(x[labels == 1].mean() - x[labels == 0].mean())
        perm_stats = []
        for shift in range(len(labels)):
            lp = np.roll(labels, shift)
            perm_stats.append(float(x[lp == 1].mean() - x[lp == 0].mean()))
        perm_stats = np.asarray(perm_stats)
        p = float(np.mean(np.abs(perm_stats) >= abs(observed) - 1e-12))
        event_rows.append({
            "Analysis": factor,
            "Unit": "Sampling month after pond centering",
            "n_event_months": int((labels == 1).sum()),
            "n_non_event_months": int((labels == 0).sum()),
            "Effect": observed,
            "Effect_definition": "Difference in mean pond-centered Shannon (event - non-event)",
            "p_value": p,
            "Method": "Two-sided exact circular-shift test",
            "Permutation_scheme": f"All {len(labels)} circular shifts of the monthly exposure sequence",
        })
    q = bh_adjust([r["p_value"] for r in event_rows])
    for row, qv in zip(event_rows, q):
        row["q_value"] = qv
        tests.append({
            "Analysis": row["Analysis"],
            "Unit": row["Unit"],
            "n": row["n_event_months"] + row["n_non_event_months"],
            "Statistic": np.nan,
            "Effect": row["Effect"],
            "Effect_definition": row["Effect_definition"],
            "p_value": row["p_value"],
            "q_value": row["q_value"],
            "Method": row["Method"],
            "Permutation_scheme": row["Permutation_scheme"],
        })

    pd.DataFrame(tests).to_csv(TABLES / "shannon_dependence_aware_tests.csv", index=False)
    pd.DataFrame(event_rows).to_csv(TABLES / "shannon_disturbance_month_tests.csv", index=False)
    paired.to_csv(TABLES / "shannon_complete_layer_pairs.csv", index=False)
    monthly.to_csv(TABLES / "shannon_month_level_analysis_data.csv", index=False)


def run_turnover(pm: pd.DataFrame, pm_rel: pd.DataFrame):
    rows = []
    for pond, group in pm.sort_values("Date").groupby("Pond"):
        group = group.sort_values("Date")
        for i in range(1, len(group)):
            previous = group.iloc[i - 1]
            current = group.iloc[i]
            a = pm_rel.loc[previous["PondMonth"]].to_numpy(dtype=float)
            b = pm_rel.loc[current["PondMonth"]].to_numpy(dtype=float)
            denom = np.sum(a + b)
            bc = np.sum(np.abs(a - b)) / denom if denom > 0 else np.nan
            rows.append({
                "Pond": pond,
                "PreviousDate": int(previous["Date"]),
                "Date": int(current["Date"]),
                "Season": current["Season"],
                "Turnover": bc,
            })
    turnover = pd.DataFrame(rows)
    grow = turnover[
        turnover["Pond"].ne("T1") & turnover["Season"].notna()
    ].sort_values(["Pond", "Date"]).reset_index(drop=True)
    result = restricted_univariate_test(
        grow, "Turnover", "Season", ["Pond"], "within_pond_circular"
    )
    pd.DataFrame([{
        "Analysis": "Seasonal difference in grow-out turnover",
        "Unit": "Successive available pond-month transition",
        "n": len(grow),
        "F": result["F"],
        "partial_R2": result["partial_R2"],
        "p_value": result["p_value"],
        "permutations": result["permutations"],
        "scheme": result["scheme"],
    }]).to_csv(TABLES / "turnover_dependence_aware_test.csv", index=False)
    turnover.to_csv(TABLES / "turnover_values.csv", index=False)

    sns.set_theme(style="whitegrid", font_scale=0.95)
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 9), gridspec_kw={"height_ratios": [1.35, 1]})
    ax = axes[0, :]
    axes[0, 0].remove()
    axes[0, 1].remove()
    top = fig.add_subplot(2, 1, 1)
    palette = {"T1": "#7B61A8", "T2": "#E76F7A", "T3": "#F39C12", "T4": "#355F8A"}
    for pond, sub in turnover.groupby("Pond"):
        top.plot(pd.to_datetime(sub["Date"].astype(str), format="%Y%m"), sub["Turnover"], marker="o", lw=1.5, label=pond, color=palette[pond])
    top.set_title("(a) Consecutive-sampling community turnover", loc="left", fontweight="bold")
    top.set_ylabel("Bray-Curtis dissimilarity")
    top.legend(title="Pond", frameon=False, bbox_to_anchor=(1.02, 0.5), loc="center left")
    sns.boxplot(
        data=turnover, x="Pond", y="Turnover", hue="Pond",
        order=["T1", "T2", "T3", "T4"], palette=palette,
        ax=axes[1, 0], width=0.6, legend=False,
    )
    axes[1, 0].set_title("(b) Turnover by pond", loc="left", fontweight="bold")
    axes[1, 0].set_ylabel("Bray-Curtis dissimilarity")
    sns.boxplot(data=grow, x="Season", y="Turnover", order=["Spring", "Summer", "Autumn", "Winter"], ax=axes[1, 1], width=0.6)
    axes[1, 1].set_title("(c) Grow-out turnover by season", loc="left", fontweight="bold")
    axes[1, 1].set_ylabel("Bray-Curtis dissimilarity")
    axes[1, 1].text(
        0.97, 0.95,
        f"Restricted circular permutation\n$p$ = {result['p_value']:.3f}",
        ha="right", va="top", transform=axes[1, 1].transAxes,
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 2},
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "Fig5_turnover_dependence_aware.png", dpi=400, bbox_inches="tight")
    plt.close(fig)


def subset_distance(rel: pd.DataFrame, ids: list[str]) -> np.ndarray:
    return squareform(pdist(rel.loc[ids].to_numpy(dtype=float), metric="braycurtis"))


def run_permanova(meta: pd.DataFrame, rel_df: pd.DataFrame, pm: pd.DataFrame, pm_rel: pd.DataFrame):
    results = []

    # Pond effect: descriptive marginal R2 only.
    for label, frame, rel in [
        ("All ponds", pm, pm_rel),
        ("Grow-out ponds", pm[pm["Pond"].ne("T1")].reset_index(drop=True), pm_rel),
    ]:
        ids = frame["PondMonth"].tolist()
        d = subset_distance(rel, ids)
        results.append({
            "Dataset": label,
            "Factor": "Pond",
            "Unit": "Pond-month mean community",
            "n": len(frame),
            "F": np.nan,
            "Marginal_R2": marginal_r2_distance(d, frame["Pond"]),
            "Partial_R2": np.nan,
            "Incremental_total_R2": np.nan,
            "p_value": np.nan,
            "q_value": np.nan,
            "Permutations": 0,
            "Restriction": "Descriptive only; pond identity was not permuted",
        })

    grow = pm[pm["Pond"].ne("T1") & pm["Season"].notna()].sort_values(["Pond", "Date"]).reset_index(drop=True)
    dg = subset_distance(pm_rel, grow["PondMonth"].tolist())
    season = restricted_permanova(
        grow, dg, "Season", ["Pond"], "within_pond_circular"
    )
    season_dispersion = restricted_dispersion_test(
        grow, dg, "Season", "within_pond_circular"
    )
    results.append({
        "Dataset": "Grow-out ponds",
        "Factor": "Season",
        "Unit": "Pond-month mean community",
        "n": len(grow),
        "F": season["F"],
        "Marginal_R2": season["marginal_R2"],
        "Partial_R2": season["partial_R2"],
        "Incremental_total_R2": season["incremental_total_R2"],
        "p_value": season["p_value"],
        "q_value": np.nan,
        "Permutations": season["permutations"],
        "Restriction": "Circular shifts of season labels within pond",
    })

    # Paired layer test.
    counts = meta.groupby("PondMonth")["Layer"].nunique()
    pair_ids = counts[counts.eq(2)].index
    layer_meta = meta[meta["PondMonth"].isin(pair_ids)].sort_values(["PondMonth", "Layer"]).reset_index(drop=True)
    dl = subset_distance(rel_df, layer_meta["SampleID"].tolist())
    layer = restricted_permanova(
        layer_meta, dl, "Layer", ["PondMonth"], "paired_layer_swap"
    )
    layer_dispersion = restricted_dispersion_test(
        layer_meta, dl, "Layer", "paired_layer_swap"
    )
    results.append({
        "Dataset": "All ponds",
        "Factor": "Layer",
        "Unit": "Complete surface-bottom pair",
        "n": len(layer_meta),
        "F": layer["F"],
        "Marginal_R2": layer["marginal_R2"],
        "Partial_R2": layer["partial_R2"],
        "Incremental_total_R2": layer["incremental_total_R2"],
        "p_value": layer["p_value"],
        "q_value": np.nan,
        "Permutations": layer["permutations"],
        "Restriction": "Surface/bottom labels swapped only within pond-month",
    })

    event_cols = ["Typhoon", "ColdSurge", "Heatwave_Tmax_q90_3day", "LowPressure"]
    complete = pm[
        pm[event_cols].notna().all(axis=1) & pm["Season"].notna()
    ].copy()
    event_indices = []
    for dataset, frame in [
        ("All ponds", complete),
        ("Grow-out ponds", complete[complete["Pond"].ne("T1")]),
    ]:
        frame = frame.sort_values(["Pond", "Date"]).reset_index(drop=True)
        d = subset_distance(pm_rel, frame["PondMonth"].tolist())
        for factor in event_cols:
            r = restricted_permanova(
                frame, d, factor, ["Pond", "Season"], "date_circular"
            )
            results.append({
                "Dataset": dataset,
                "Factor": factor,
                "Unit": "Pond-month mean community; sampling month permuted as a block",
                "n": len(frame),
                "F": r["F"],
                "Marginal_R2": r["marginal_R2"],
                "Partial_R2": r["partial_R2"],
                "Incremental_total_R2": r["incremental_total_R2"],
                "p_value": r["p_value"],
                "q_value": np.nan,
                "Permutations": r["permutations"],
                "Restriction": "Circular shifts of the complete monthly exposure sequence",
            })
            event_indices.append(len(results) - 1)

    # Control multiplicity separately within each four-indicator dataset.
    for dataset in ["All ponds", "Grow-out ponds"]:
        idx = [i for i in event_indices if results[i]["Dataset"] == dataset]
        q = bh_adjust([results[i]["p_value"] for i in idx])
        for i, qv in zip(idx, q):
            results[i]["q_value"] = qv

    pd.DataFrame(results).to_csv(TABLES / "permanova_dependence_aware_results.csv", index=False)
    pd.DataFrame([
        {
            "Dataset": "Grow-out ponds",
            "Factor": "Season",
            "n": season_dispersion["n"],
            "F": season_dispersion["F"],
            "p_value": season_dispersion["p_value"],
            "permutations": season_dispersion["permutations"],
            "restriction": "Circular shifts of season labels within pond",
        },
        {
            "Dataset": "All ponds",
            "Factor": "Layer",
            "n": layer_dispersion["n"],
            "F": layer_dispersion["F"],
            "p_value": layer_dispersion["p_value"],
            "permutations": layer_dispersion["permutations"],
            "restriction": "Surface/bottom labels swapped only within pond-month",
        },
    ]).to_csv(TABLES / "permdisp_dependence_aware_results.csv", index=False)


def run_vectors_and_blocks(pm: pd.DataFrame, pm_rel: pd.DataFrame):
    env_vars_water = ["Temperature", "Salinity", "pH", "DO", "Turbidity"]
    env_vars_climate = ["MeanPressure", "MinPressure", "Tmax", "TotalRainfall", "MaxDailyRainfall"]
    event_cols = ["Typhoon", "ColdSurge", "Heatwave_Tmax_q90_3day", "LowPressure"]
    vector_needed = env_vars_water + env_vars_climate
    model_needed = vector_needed + event_cols + ["Season"]

    # Vector fitting itself does not condition on season; retaining the three
    # boundary-spanning grow-out profiles preserves the stated 73-profile
    # environmental cohort.  Season-dependent models below use only the 70
    # unambiguous grow-out pond-months.
    vector_df = (
        pm[pm["Pond"].ne("T1")]
        .dropna(subset=vector_needed)
        .sort_values(["Pond", "Date"])
        .reset_index(drop=True)
    )
    d = subset_distance(pm_rel, vector_df["PondMonth"].tolist())
    coords, eigvals, explained = pcoa(d)
    coords2 = coords[:, :2]
    coord_df = vector_df[["PondMonth", "Pond", "Date", "Season"]].copy()
    coord_df["PCoA1"] = coords2[:, 0]
    coord_df["PCoA2"] = coords2[:, 1]
    coord_df.to_csv(TABLES / "pond_month_pcoa_coordinates_T2_T4.csv", index=False)

    vector_statistic = make_vector_statistic(coords2)
    rows = []
    for variable in env_vars_water + env_vars_climate:
        r1, r2, observed = vector_statistic(vector_df[variable].to_numpy())
        rng = make_rng(f"vector|{variable}|{len(vector_df)}|{N_PERM}")
        if variable in env_vars_water:
            scheme = "within_pond_circular"
            restriction = "Circular shifts within pond"
        else:
            scheme = "date_circular"
            restriction = "Circular shifts of shared sampling-month values"
        permute = make_continuous_permuter(vector_df, variable, scheme, rng)
        exceed = 0
        for _ in range(N_PERM):
            _, _, rp = vector_statistic(permute())
            exceed += int(rp >= observed - 1e-12)
        rows.append({
            "Environmental_variable": variable,
            "Correlation_PCoA1": r1,
            "Correlation_PCoA2": r2,
            "Vector_fit_R2": observed,
            "p_value": (exceed + 1) / (N_PERM + 1),
            "q_value": np.nan,
            "n_pond_months": len(vector_df),
            "permutations": N_PERM,
            "restriction": restriction,
        })
    q = bh_adjust([r["p_value"] for r in rows])
    for row, qv in zip(rows, q):
        row["q_value"] = qv
    vector_results = pd.DataFrame(rows)
    vector_results.to_csv(TABLES / "environmental_vector_fits_dependence_aware.csv", index=False)

    # Plot pond-month PCoA and fitted environmental arrows.
    fig, ax = plt.subplots(figsize=(8.2, 6.8))
    palette = {"Spring": "#6BAF92", "Summer": "#F2A65A", "Autumn": "#C67C5C", "Winter": "#6C8EBF"}
    for season, sub in coord_df.groupby("Season"):
        ax.scatter(sub["PCoA1"], sub["PCoA2"], s=45, alpha=0.8, label=season, color=palette.get(season))
    boundary = coord_df[coord_df["Season"].isna()]
    if len(boundary):
        ax.scatter(
            boundary["PCoA1"], boundary["PCoA2"], s=48, alpha=0.85,
            label="Boundary-spanning average", color="#777777", marker="X",
        )
    scale = 0.8 * min(np.ptp(coords2[:, 0]), np.ptp(coords2[:, 1]))
    label_offsets = {
        "Temperature": (-0.006, 0.004, "right"),
        "Salinity": (0.008, 0.006, "left"),
        "pH": (0.020, 0.010, "left"),
        "DO": (0.008, -0.010, "left"),
        "Turbidity": (-0.004, 0.010, "right"),
        "MeanPressure": (-0.004, -0.026, "center"),
        "MinPressure": (0.010, -0.014, "left"),
        "Tmax": (-0.010, 0.010, "right"),
        "TotalRainfall": (-0.006, 0.018, "right"),
        "MaxDailyRainfall": (0.008, 0.016, "left"),
    }
    for _, row in vector_results.iterrows():
        x = row["Correlation_PCoA1"] * scale
        y = row["Correlation_PCoA2"] * scale
        significant = row["q_value"] < 0.05
        color = "#A61B1B" if significant else "#3F3F3F"
        ax.arrow(0, 0, x, y, color=color, alpha=0.82, width=0.00035, head_width=0.008, length_includes_head=True)
        dx, dy, ha = label_offsets[row["Environmental_variable"]]
        display_label = {
            "Temperature": "Water temperature",
            "Salinity": "Salinity",
            "DO": "Dissolved oxygen",
            "Tmax": "Maximum air temperature",
            "TotalRainfall": "Total rainfall",
            "MaxDailyRainfall": "Maximum daily rainfall",
            "MeanPressure": "Mean pressure",
            "MinPressure": "Minimum pressure",
        }.get(row["Environmental_variable"], row["Environmental_variable"])
        ax.text(
            x + dx, y + dy, display_label,
            fontsize=8.2, ha=ha, va="center", color=color,
            fontweight="bold" if significant else "normal",
        )
    ax.axhline(0, color="#AAAAAA", lw=0.8)
    ax.axvline(0, color="#AAAAAA", lw=0.8)
    ax.set_xlabel(f"PCoA1 ({explained[0]:.1f}%)")
    ax.set_ylabel(f"PCoA2 ({explained[1]:.1f}%)")
    ax.set_title("Environmental alignment in grow-out ponds (pond-month means)", fontweight="bold")
    ax.legend(title="Season", frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(FIGURES / "Fig7_environmental_vectors_dependence_aware.png", dpi=400, bbox_inches="tight")
    plt.close(fig)

    blocks = {
        "Pond identity": (["Pond"], []),
        "Season": (["Season"], []),
        "Water quality": ([], env_vars_water),
        "Meteorological disturbances": (event_cols, []),
    }
    block_rows = []
    for dataset, frame in [
        ("All ponds", pm.dropna(subset=model_needed)),
        (
            "Grow-out ponds",
            pm[pm["Pond"].ne("T1")].dropna(subset=model_needed),
        ),
    ]:
        frame = frame.sort_values(["Pond", "Date"]).reset_index(drop=True)
        dd = subset_distance(pm_rel, frame["PondMonth"].tolist())
        full, unique = blockwise_r2(frame, dd, blocks)
        for block, value in unique.items():
            block_rows.append({
                "Dataset": dataset,
                "n_pond_months": len(frame),
                "Full_model_R2": full,
                "Predictor_block": block,
                "Unique_R2": value,
            })
    pd.DataFrame(block_rows).to_csv(TABLES / "blockwise_unique_R2_pond_month.csv", index=False)

    return {
        "vector_n": len(vector_df),
        "pcoa1_explained": explained[0],
        "pcoa2_explained": explained[1],
    }


def run_season_sensitivity(meta: pd.DataFrame, pm: pd.DataFrame, pm_rel: pd.DataFrame):
    """Record the primary season cohort after boundary-spanning averages are excluded."""
    affected = set(meta.loc[meta["Season"].isna(), "PondMonth"])
    grow = (
        pm[pm["Pond"].ne("T1") & ~pm["PondMonth"].isin(affected)]
        .sort_values(["Pond", "Date"])
        .reset_index(drop=True)
    )
    shannon = restricted_univariate_test(
        grow, "Shannon", "Season", ["Pond"], "within_pond_circular"
    )
    d = subset_distance(pm_rel, grow["PondMonth"].tolist())
    community = restricted_permanova(
        grow, d, "Season", ["Pond"], "within_pond_circular"
    )

    turnover = pd.read_csv(TABLES / "turnover_values.csv")
    affected_keys = {
        (row.Pond, int(str(row.PondMonth).split("-")[1]))
        for row in meta.loc[meta["Season"].isna(), ["Pond", "PondMonth"]].drop_duplicates().itertuples()
    }
    turn = turnover[
        turnover["Pond"].ne("T1")
        & ~turnover.apply(lambda x: (x["Pond"], int(x["Date"])) in affected_keys, axis=1)
    ].sort_values(["Pond", "Date"]).reset_index(drop=True)
    turnover_result = restricted_univariate_test(
        turn, "Turnover", "Season", ["Pond"], "within_pond_circular"
    )

    pd.DataFrame([
        {
            "Response": "Shannon diversity",
            "n": len(grow),
            "F": shannon["F"],
            "partial_R2": shannon["partial_R2"],
            "p_value": shannon["p_value"],
        },
        {
            "Response": "Community composition",
            "n": len(grow),
            "F": community["F"],
            "partial_R2": community["partial_R2"],
            "p_value": community["p_value"],
        },
        {
            "Response": "Consecutive-sampling turnover",
            "n": len(turn),
            "F": turnover_result["F"],
            "partial_R2": turnover_result["partial_R2"],
            "p_value": turnover_result["p_value"],
        },
    ]).assign(
        Sensitivity_definition=(
            "Primary season cohort excludes pond-month profiles that average "
            "sampling occasions assigned to different solar-term seasons"
        ),
        Excluded_pond_months="; ".join(sorted(affected)),
        Permutations=N_PERM,
        Restriction="Circular shifts within pond",
    ).to_csv(TABLES / "season_boundary_exclusion_audit.csv", index=False)


def write_quality_report(meta, pm, low_threshold, extra):
    report = {
        "input_files": {
            "ASV": ASV_FILE.name,
            "metadata": META_FILE.name,
            "water": WATER_FILE.name,
            "climate": CLIMATE_FILE.name,
        },
        "sample_rows": int(len(meta)),
        "pond_months": int(len(pm)),
        "complete_surface_bottom_pairs": int((meta.groupby("PondMonth")["Layer"].nunique() == 2).sum()),
        "pond_counts": {k: int(v) for k, v in meta["Pond"].value_counts().items()},
        "missing_season_original_metadata": int(
            meta.get("Season_original_metadata", meta["Season"]).isna().sum()
        ),
        "boundary_spanning_sample_rows_excluded_from_season_models": int(
            meta["Season_verified"].isna().sum()
        ),
        "unambiguous_season_sample_rows": int(meta["Season_verified"].notna().sum()),
        "complete_disturbance_sample_rows": int(meta["DisturbanceComplete"].sum()),
        "complete_disturbance_pond_months": int(
            pm[["Typhoon", "ColdSurge", "Heatwave_Tmax_q90_3day", "LowPressure"]]
            .notna().all(axis=1).sum()
        ),
        "low_pressure_threshold_hPa": low_threshold,
        "random_seed": SEED,
        "permutations": N_PERM,
        **extra,
    }
    (OUT / "analysis_quality_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")


def main():
    meta, asv, rel_df, climate, low_threshold = validate_and_load()
    pm, pm_rel = aggregate_pond_month(meta, rel_df)
    meta.to_csv(TABLES / "validated_sample_metadata_with_shannon.csv", index=False)
    pm.to_csv(TABLES / "pond_month_analysis_metadata.csv", index=False)

    run_shannon(meta, pm)
    run_turnover(pm, pm_rel)
    run_permanova(meta, rel_df, pm, pm_rel)
    extra = run_vectors_and_blocks(pm, pm_rel)
    run_season_sensitivity(meta, pm, pm_rel)
    write_quality_report(meta, pm, low_threshold, extra)
    print(f"Analysis completed. Outputs: {OUT}")


if __name__ == "__main__":
    main()
