#!/usr/bin/env python3
"""Reproduce the library-depth and repeated-rarefaction sensitivity analysis.

The analysis follows the manuscript exactly:

* start from the 190 formal biological sequencing libraries;
* remove ASVs assigned to mitochondria, chloroplasts, or explicit eukaryotes;
* rarefy every eligible source library independently at seven depths;
* retain a canonical pond-month-layer profile only if all of its source
  libraries meet the requested depth;
* average source-library relative-abundance profiles before calculating
  Shannon diversity; and
* compare paired bottom and surface profiles with a two-sided Wilcoxon
  signed-rank test.

The reference summaries submitted with the manuscript are retained under
``data/supporting`` and ``reference_results/sensitivity``.  The fixed seed
makes future reruns deterministic.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon


DEFAULT_DEPTHS = (1500, 3000, 5000, 10000, 15000, 20000, 25000)
ORGANELLE_EUKARYOTE_PATTERN = r"mitochond|chloroplast|d__eukary"


def shannon(relative: np.ndarray) -> float:
    positive = relative[relative > 0]
    return float(-(positive * np.log(positive)).sum())


def rarefy_counts(
    counts: np.ndarray, depth: int, rng: np.random.Generator
) -> np.ndarray:
    """Draw without replacement from a sparse integer count vector."""

    nonzero = np.flatnonzero(counts)
    draw = rng.multivariate_hypergeometric(counts[nonzero], depth)
    relative = np.zeros(counts.size, dtype=np.float64)
    relative[nonzero] = draw / float(depth)
    return relative


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--library-asv", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--repeats", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260728)
    parser.add_argument(
        "--depths", type=int, nargs="+", default=list(DEFAULT_DEPTHS)
    )
    args = parser.parse_args()

    outdir = args.output_dir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # The workflow feature table contains a QIIME 2 semantic-type row.
    feature = pd.read_csv(args.library_asv, skiprows=[1], low_memory=False)
    metadata = pd.read_csv(args.metadata)
    libraries = metadata["SampleID_raw"].astype(str).tolist()
    missing = [library for library in libraries if library not in feature.columns]
    if missing:
        raise ValueError(f"Feature table is missing formal libraries: {missing}")

    taxonomy = feature["Taxon"].fillna("").astype(str)
    positive = feature[libraries].fillna(0).sum(axis=1) > 0
    excluded = taxonomy.str.contains(
        ORGANELLE_EUKARYOTE_PATTERN, case=False, regex=True
    )
    retained = feature.loc[positive & ~excluded, libraries]
    counts = retained.T.to_numpy(dtype=np.int64)
    totals = counts.sum(axis=1)
    if counts.shape[0] != 190 or counts.shape[1] != 23896:
        raise ValueError(
            "Expected a 190-library by 23,896-ASV prokaryotic matrix; "
            f"found {counts.shape}"
        )

    meta = metadata.set_index("SampleID_raw").loc[libraries].reset_index()
    library_position = {library: i for i, library in enumerate(libraries)}
    profile_groups = {
        profile: group["SampleID_raw"].astype(str).tolist()
        for profile, group in meta.groupby("SampleID", sort=False)
    }
    profile_info = (
        meta.drop_duplicates("SampleID")
        .set_index("SampleID")[["Pond", "Date", "Layer"]]
    )

    sample_rows: list[dict[str, object]] = []
    aggregate_rows: list[dict[str, object]] = []
    iteration_rows: list[dict[str, object]] = []

    seed_sequence = np.random.SeedSequence(args.seed)
    depth_sequences = seed_sequence.spawn(len(args.depths))

    for depth, depth_seed in zip(args.depths, depth_sequences):
        all_depth_eligible_libraries = [
            library
            for library in libraries
            if totals[library_position[library]] >= depth
        ]
        eligible_profiles = {
            profile: source_libraries
            for profile, source_libraries in profile_groups.items()
            if all(totals[library_position[x]] >= depth for x in source_libraries)
        }
        eligible_libraries = sorted(
            {x for source in eligible_profiles.values() for x in source}
        )
        repeat_values: dict[str, list[float]] = {
            profile: [] for profile in eligible_profiles
        }
        per_repeat_statistics: list[dict[str, float | int]] = []

        repeat_sequences = depth_seed.spawn(args.repeats)
        for repeat_index, repeat_seed in enumerate(repeat_sequences, start=1):
            rng = np.random.default_rng(repeat_seed)
            rarefied: dict[str, np.ndarray] = {}
            for library in eligible_libraries:
                position = library_position[library]
                rarefied[library] = rarefy_counts(counts[position], depth, rng)

            profile_shannon: dict[str, float] = {}
            for profile, source_libraries in eligible_profiles.items():
                averaged = np.mean(
                    [rarefied[library] for library in source_libraries], axis=0
                )
                value = shannon(averaged)
                repeat_values[profile].append(value)
                profile_shannon[profile] = value

            paired = []
            for (pond, date), group in profile_info.loc[
                list(eligible_profiles)
            ].groupby(["Pond", "Date"]):
                by_layer = {
                    row.Layer: profile_shannon[profile]
                    for profile, row in group.iterrows()
                }
                if "D" in by_layer and "U" in by_layer:
                    paired.append(by_layer["D"] - by_layer["U"])
            differences = np.asarray(paired, dtype=float)
            test = wilcoxon(differences, alternative="two-sided")
            per_repeat_statistics.append(
                {
                    "Depth": depth,
                    "Repeat": repeat_index,
                    "Pairs": len(differences),
                    "Median_bottom_minus_surface": float(
                        np.median(differences)
                    ),
                    "Wilcoxon_W": float(test.statistic),
                    "p_value": float(test.pvalue),
                }
            )

        for profile, values in repeat_values.items():
            source_libraries = eligible_profiles[profile]
            row = profile_info.loc[profile]
            sample_rows.append(
                {
                    "Rarefaction depth per original library": depth,
                    "CanonicalSampleID": profile,
                    "Pond": row["Pond"],
                    "PondMonth": f"{row['Pond']}-{int(row['Date'])}",
                    "Layer": row["Layer"],
                    "Source libraries": ";".join(source_libraries),
                    "Number of source libraries": len(source_libraries),
                    "Mean rarefied Shannon": float(np.mean(values)),
                    "SD across rarefaction repeats": float(
                        np.std(values, ddof=1)
                    ),
                    "Rarefaction repeats": args.repeats,
                }
            )

        mean_scores = {
            profile: float(np.mean(values))
            for profile, values in repeat_values.items()
        }
        paired_means = []
        for (pond, date), group in profile_info.loc[
            list(eligible_profiles)
        ].groupby(["Pond", "Date"]):
            by_layer = {
                row.Layer: mean_scores[profile]
                for profile, row in group.iterrows()
            }
            if "D" in by_layer and "U" in by_layer:
                paired_means.append(by_layer["D"] - by_layer["U"])
        paired_means_array = np.asarray(paired_means, dtype=float)
        mean_test = wilcoxon(paired_means_array, alternative="two-sided")
        iteration = pd.DataFrame(per_repeat_statistics)
        aggregate_rows.append(
            {
                "Depth_per_library": depth,
                # This count is the number of source libraries meeting the
                # depth threshold, including libraries whose paired source
                # record caused their canonical profile to be excluded.
                "Libraries_retained": len(all_depth_eligible_libraries),
                "Canonical_samples": len(eligible_profiles),
                "Pairs": len(paired_means_array),
                "Median_bottom_minus_surface": float(
                    np.median(paired_means_array)
                ),
                "Wilcoxon_W": float(mean_test.statistic),
                "p_value": float(mean_test.pvalue),
                "Iterations_with_p_below_0.05": int(
                    (iteration["p_value"] < 0.05).sum()
                ),
                "Iterations": args.repeats,
                "Seed": args.seed,
            }
        )
        iteration_rows.extend(per_repeat_statistics)

    pd.DataFrame(sample_rows).to_csv(
        outdir / "rarefaction_sample_summary.csv", index=False
    )
    pd.DataFrame(aggregate_rows).to_csv(
        outdir / "rarefaction_aggregate_summary.csv", index=False
    )
    pd.DataFrame(iteration_rows).to_csv(
        outdir / "rarefaction_iteration_tests.csv", index=False
    )
    print(pd.DataFrame(aggregate_rows).to_string(index=False))


if __name__ == "__main__":
    main()
