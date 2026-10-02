#!/usr/bin/env python3
"""Compare primary and conservative blank-screened PICRUSt2 results.

The single Water blank is not used as a biological sample and is insufficient
for formal prevalence-based contaminant calling.  This script therefore treats
the conservative blank-enriched ASV screen only as a sensitivity analysis.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from scipy.stats import pearsonr, spearmanr


def read_profiles(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", compression="infer", index_col=0)
    profiles = table.T.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    row_sums = profiles.sum(axis=1)
    if (row_sums <= 0).any():
        missing = profiles.index[row_sums <= 0].tolist()
        raise ValueError(f"Zero-sum pathway profiles: {missing[:5]}")
    return profiles.div(row_sums, axis=0)


def finite_correlation(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 3:
        return {"n": int(keep.sum()), "pearson_r": np.nan, "spearman_rho": np.nan}
    return {
        "n": int(keep.sum()),
        "pearson_r": float(pearsonr(x[keep], y[keep]).statistic),
        "spearman_rho": float(spearmanr(x[keep], y[keep]).statistic),
    }


def compare_test_tables(primary: Path, screened: Path) -> tuple[pd.DataFrame, dict]:
    a = pd.read_csv(primary)
    b = pd.read_csv(screened)
    candidate_keys = [
        "Test_family",
        "Response",
        "Factor",
        "Dataset",
        "Unit",
        "Analysis",
        "Test",
        "Variable",
        "Module",
    ]
    keys = [column for column in candidate_keys if column in a.columns and column in b.columns]
    if not keys:
        a = a.reset_index().rename(columns={"index": "Row"})
        b = b.reset_index().rename(columns={"index": "Row"})
        keys = ["Row"]
    if a.duplicated(keys).any() or b.duplicated(keys).any():
        raise ValueError(
            "The selected test identifiers are not unique within both result "
            f"tables: {keys}"
        )
    merged = a.merge(b, on=keys, how="outer", suffixes=("_primary", "_screened"))

    for metric in ["p_value", "q_value", "P", "Q", "Partial_R2", "Marginal_R2", "Pseudo_F"]:
        pcol = f"{metric}_primary"
        scol = f"{metric}_screened"
        if pcol in merged and scol in merged:
            merged[f"{metric}_absolute_change"] = (
                pd.to_numeric(merged[scol], errors="coerce")
                - pd.to_numeric(merged[pcol], errors="coerce")
            ).abs()

    def effective_significance(suffix: str) -> pd.Series:
        q_column = (
            f"q_value_{suffix}"
            if f"q_value_{suffix}" in merged
            else f"Q_{suffix}"
        )
        p_column = (
            f"p_value_{suffix}"
            if f"p_value_{suffix}" in merged
            else f"P_{suffix}"
        )
        q_values = (
            pd.to_numeric(merged[q_column], errors="coerce")
            if q_column in merged
            else pd.Series(np.nan, index=merged.index)
        )
        p_values = (
            pd.to_numeric(merged[p_column], errors="coerce")
            if p_column in merged
            else pd.Series(np.nan, index=merged.index)
        )
        return q_values.where(q_values.notna(), p_values)

    pa = effective_significance("primary")
    pb = effective_significance("screened")
    keep = pa.notna() & pb.notna()
    sig = {
        "significance_metric": (
            "Benjamini–Hochberg q-value for adjusted test families; "
            "otherwise the prespecified p-value"
        ),
        "comparisons_with_both_values": int(keep.sum()),
        "significance_class_concordance_at_0.05": float(
            ((pa[keep] < 0.05) == (pb[keep] < 0.05)).mean()
        )
        if keep.any()
        else np.nan,
        "classification_changes": int(
            ((pa[keep] < 0.05) != (pb[keep] < 0.05)).sum()
        ),
    }
    return merged, sig


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--primary-picrust", required=True, type=Path)
    parser.add_argument("--screened-picrust", required=True, type=Path)
    parser.add_argument("--primary-results", required=True, type=Path)
    parser.add_argument("--screened-results", required=True, type=Path)
    args = parser.parse_args()

    qc_dir = args.project / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)

    primary = read_profiles(
        args.primary_picrust / "pathways_out" / "path_abun_unstrat.tsv.gz"
    )
    screened = read_profiles(
        args.screened_picrust / "pathways_out" / "path_abun_unstrat.tsv.gz"
    )
    common_samples = primary.index.intersection(screened.index)
    common_pathways = primary.columns.union(screened.columns)
    primary = primary.reindex(
        index=common_samples, columns=common_pathways, fill_value=0.0
    )
    screened = screened.reindex(
        index=common_samples, columns=common_pathways, fill_value=0.0
    )

    primary_dist = pdist(primary.to_numpy(), metric="braycurtis")
    screened_dist = pdist(screened.to_numpy(), metric="braycurtis")
    distance_correlation = finite_correlation(primary_dist, screened_dist)

    numerator = np.abs(primary.to_numpy() - screened.to_numpy()).sum(axis=1)
    denominator = (primary.to_numpy() + screened.to_numpy()).sum(axis=1)
    per_sample_bc = np.divide(
        numerator,
        denominator,
        out=np.zeros_like(numerator),
        where=denominator > 0,
    )
    sample_comparison = pd.DataFrame(
        {
            "SampleID": common_samples,
            "Bray_Curtis_primary_vs_blank_screened": per_sample_bc,
        }
    )
    sample_comparison.to_csv(
        qc_dir / "blank_screened_functional_profile_sensitivity_by_sample.csv",
        index=False,
    )

    test_comparison, significance_summary = compare_test_tables(
        args.primary_results / "tables" / "Table_S10_dependence_aware_functional_tests.csv",
        args.screened_results
        / "tables"
        / "Table_S10_dependence_aware_functional_tests.csv",
    )
    test_comparison.to_csv(
        qc_dir / "blank_screened_functional_test_sensitivity.csv", index=False
    )

    primary_modules = pd.read_csv(
        args.primary_results / "tables" / "functional_module_scores_sample.csv",
        index_col=0,
    )
    screened_modules = pd.read_csv(
        args.screened_results / "tables" / "functional_module_scores_sample.csv",
        index_col=0,
    )
    common_module_samples = primary_modules.index.intersection(screened_modules.index)
    common_modules = primary_modules.columns.intersection(screened_modules.columns)
    module_rows = []
    for module in common_modules:
        x = pd.to_numeric(
            primary_modules.loc[common_module_samples, module], errors="coerce"
        ).to_numpy()
        y = pd.to_numeric(
            screened_modules.loc[common_module_samples, module], errors="coerce"
        ).to_numpy()
        corr = finite_correlation(x, y)
        module_rows.append(
            {
                "Module": module,
                **corr,
                "median_absolute_score_change": float(
                    np.nanmedian(np.abs(x - y))
                ),
                "maximum_absolute_score_change": float(
                    np.nanmax(np.abs(x - y))
                ),
            }
        )
    module_comparison = pd.DataFrame(module_rows)
    module_comparison.to_csv(
        qc_dir / "blank_screened_module_score_sensitivity.csv", index=False
    )

    summary = {
        "interpretation": (
            "Sensitivity analysis only: the single Water blank was excluded "
            "from biological analyses, and ASVs conservatively flagged as "
            "blank-enriched were removed without claiming formal contaminant status."
        ),
        "common_samples": int(len(common_samples)),
        "union_pathways": int(len(common_pathways)),
        "pairwise_pathway_bray_curtis_correlation": distance_correlation,
        "within_sample_primary_vs_screened_bray_curtis": {
            "median": float(np.median(per_sample_bc)),
            "q1": float(np.quantile(per_sample_bc, 0.25)),
            "q3": float(np.quantile(per_sample_bc, 0.75)),
            "maximum": float(np.max(per_sample_bc)),
        },
        "functional_test_significance": significance_summary,
        "module_score_minimum_pearson_r": float(
            module_comparison["pearson_r"].min()
        ),
        "module_score_minimum_spearman_rho": float(
            module_comparison["spearman_rho"].min()
        ),
    }
    (qc_dir / "blank_screened_functional_sensitivity_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
