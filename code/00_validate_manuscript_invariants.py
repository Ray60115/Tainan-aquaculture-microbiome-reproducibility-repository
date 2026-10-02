#!/usr/bin/env python3
"""Fast audit of numerical invariants reported in the submitted manuscript.

This does not replace the 9,999-permutation reruns.  It verifies cohort
construction, feature filtering, reference-result values, and the exact
sample/column boundaries that are easiest to misstate in a public archive.
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform


ANNOTATION = ["id", "Sequence", "Taxon"]


def close(actual: float, expected: float, tolerance: float = 1e-6) -> None:
    if not np.isclose(float(actual), float(expected), rtol=0, atol=tolerance):
        raise AssertionError(f"Expected {expected}, found {actual}")


def pcoa_explained(relative: np.ndarray) -> np.ndarray:
    distance = squareform(pdist(relative, metric="braycurtis"))
    n = len(distance)
    center = np.eye(n) - np.ones((n, n)) / n
    gower = -0.5 * center @ (distance**2) @ center
    values = np.linalg.eigvalsh(gower)[::-1]
    positive = values[values > max(1e-12, values[0] * 1e-10)]
    return 100 * positive / positive.sum()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output or root / "validation" / "manuscript_invariant_audit.json"

    primary = root / "data" / "primary"
    manifests = root / "data" / "manifests"
    taxonomy_results = root / "reference_results" / "taxonomy"
    functional_results = root / "reference_results" / "functional"
    sensitivity_results = root / "reference_results" / "sensitivity"

    metadata = pd.read_csv(primary / "sample_metadata_168.csv")
    library_map = pd.read_csv(primary / "library_to_profile_metadata_190.csv")
    workflow_manifest = pd.read_csv(root / "workflow" / "sample_input_204.tsv", sep="\t")
    qc_dir = root / "workflow" / "qc"
    dada2_qc = pd.read_csv(qc_dir / "dada2_qc.tsv", sep="\t")
    dada2_qc = dada2_qc.loc[
        ~dada2_qc["sample-id"].astype(str).str.startswith("#")
    ].copy()
    read_qc = pd.read_csv(qc_dir / "read_qc_summary.tsv", sep="\t")
    trim_qc = pd.read_csv(qc_dir / "trim_qc_summary.tsv", sep="\t")
    workflow_ids = set(workflow_manifest["sample-id"].astype(str))
    assert len(dada2_qc) == 204 and set(dada2_qc["sample-id"].astype(str)) == workflow_ids
    assert len(read_qc) == 204 and set(read_qc["sample"].astype(str)) == workflow_ids
    assert len(trim_qc) == 204 and set(trim_qc["sample"].astype(str)) == workflow_ids
    strict = pd.read_csv(primary / "asv_table_168_profiles.csv.gz", low_memory=False)
    controls = pd.read_csv(primary / "asv_table_9_nonformal_columns.csv.gz", low_memory=False)
    raw = pd.read_csv(
        primary / "asv_table_workflow_feature_table_199_columns.csv.gz",
        skiprows=[1],
        low_memory=False,
    )

    samples = metadata["SampleID"].astype(str).tolist()
    formal_libraries = library_map["SampleID_raw"].astype(str).tolist()
    strict_sample_columns = [c for c in strict.columns if c not in ANNOTATION]
    control_columns = [c for c in controls.columns if c not in ANNOTATION]
    raw_sample_columns = [c for c in raw.columns if c not in ANNOTATION]
    workflow_ids = workflow_manifest["sample-id"].astype(str).tolist()

    assert len(workflow_ids) == len(set(workflow_ids)) == 204
    assert len(formal_libraries) == len(set(formal_libraries)) == 190
    assert len(samples) == len(set(samples)) == 168
    assert strict_sample_columns == samples
    assert len(control_columns) == 9 and set(control_columns) == {
        "T2-202302-1-G", "T2-202412-1-G", "T3-202302-1-G", "T3-202411-G",
        "T3-202412-1-G", "T4-202302-1-G", "T4-202411-G", "T4-202412-1-G", "Water",
    }
    assert len(raw_sample_columns) == 199
    assert set(formal_libraries).issubset(raw_sample_columns)
    assert set(workflow_ids) - set(raw_sample_columns) == {
        "Probiotic1", "Probiotic2", "Probiotic3", "SeaWater1", "SeaWater2"
    }
    assert len(strict) == len(controls) == len(raw) == 29461
    assert set(strict["id"]) == set(controls["id"]) == set(raw["id"])

    profile_counts = library_map.groupby("SampleID").size()
    assert len(profile_counts) == 168
    assert int((profile_counts == 1).sum()) == 146
    assert int((profile_counts == 2).sum()) == 22
    assert library_map["Layer"].value_counts().to_dict() == {"U": 102, "D": 88}
    assert metadata["Pond"].value_counts().to_dict() == {"T2": 47, "T4": 42, "T1": 40, "T3": 39}
    assert metadata["PondMonth"].nunique() == 93
    layer_sets = metadata.groupby("PondMonth")["Layer"].agg(set)
    assert int(layer_sets.eq({"D", "U"}).sum()) == 75
    assert int(layer_sets.eq({"U"}).sum()) == 16
    assert int(layer_sets.eq({"D"}).sum()) == 2
    assert int(metadata["Season_verified"].notna().sum()) == 160
    pm_meta = metadata.drop_duplicates("PondMonth")
    assert int(pm_meta["Season_verified"].notna().sum()) == 89
    assert int((metadata["Pond"].ne("T1") & metadata["Season_verified"].notna()).sum()) == 122
    assert int((pm_meta["Pond"].ne("T1") & pm_meta["Season_verified"].notna()).sum()) == 70

    taxonomy = strict["Taxon"].fillna("").astype(str)
    counts = strict[samples].fillna(0).astype(float)
    positive = counts.sum(axis=1) > 0
    excluded = taxonomy.str.contains(r"mitochond|chloroplast|d__eukary", case=False, regex=True)
    retained = positive & ~excluded
    assert int((~positive).sum()) == 4242
    assert int((positive & taxonomy.str.contains("mitochond", case=False)).sum()) == 728
    assert int((positive & taxonomy.str.contains("chloroplast", case=False)).sum()) == 590
    assert int((positive & taxonomy.str.contains("d__eukary", case=False)).sum()) == 5
    assert int(retained.sum()) == 23896
    domain = taxonomy.str.extract(r"(?:^|;\s*)(d__[^;]+)", expand=False).fillna("")
    assert int((retained & domain.str.casefold().eq("d__archaea")).sum()) == 15
    retained_counts = counts.loc[retained]
    archaeal_fraction = (
        counts.loc[retained & domain.str.casefold().eq("d__archaea")].to_numpy().sum()
        / retained_counts.to_numpy().sum()
    )
    close(100 * archaeal_fraction, 0.0073231219721, 5e-10)

    relative = retained_counts.T.to_numpy(float)
    relative /= relative.sum(axis=1, keepdims=True)
    explained = pcoa_explained(relative)
    close(explained[0], 6.58859603, 5e-7)
    close(explained[1], 6.36100800, 5e-7)

    shannon = pd.read_csv(taxonomy_results / "shannon_dependence_aware_tests.csv")
    layer = shannon.loc[shannon["Analysis"].eq("Sampling layer")].iloc[0]
    season = shannon.loc[shannon["Analysis"].eq("Season within grow-out ponds")].iloc[0]
    close(layer["Statistic"], 604)
    close(layer["Effect"], 0.2566888)
    close(layer["p_value"], 1.4553689069e-5, 5e-12)
    close(season["Effect"], 0.02793096)
    close(season["p_value"], 0.5759)

    permanova = pd.read_csv(taxonomy_results / "permanova_dependence_aware_results.csv")
    season_comp = permanova.loc[
        permanova["Dataset"].eq("Grow-out ponds") & permanova["Factor"].eq("Season")
    ].iloc[0]
    close(season_comp["Marginal_R2"], 0.0993539)
    close(season_comp["Partial_R2"], 0.1055753)
    close(season_comp["p_value"], 0.0002)
    layer_comp = permanova.loc[
        permanova["Dataset"].eq("All ponds") & permanova["Factor"].eq("Layer")
    ].iloc[0]
    close(layer_comp["Marginal_R2"], 0.0038076)
    close(layer_comp["Partial_R2"], 0.046264)
    close(layer_comp["p_value"], 0.0001)

    turnover = pd.read_csv(taxonomy_results / "turnover_dependence_aware_test.csv").iloc[0]
    assert int(turnover["n"]) == 67
    close(turnover["partial_R2"], 0.00498791)
    close(turnover["p_value"], 0.9629)

    vectors = pd.read_csv(taxonomy_results / "environmental_vector_fits_dependence_aware.csv")
    ph = vectors.loc[vectors["Environmental_variable"].eq("pH")].iloc[0]
    water_temp = vectors.loc[vectors["Environmental_variable"].eq("Temperature")].iloc[0]
    close(ph["Vector_fit_R2"], 0.546583)
    close(ph["q_value"], 0.001)
    close(water_temp["Vector_fit_R2"], 0.324514)
    close(water_temp["q_value"], 0.025)

    blocks = pd.read_csv(taxonomy_results / "blockwise_unique_R2_pond_month.csv")
    close(blocks.loc[blocks["Dataset"].eq("All ponds"), "Full_model_R2"].iloc[0], 0.307283)
    close(blocks.loc[blocks["Dataset"].eq("Grow-out ponds"), "Full_model_R2"].iloc[0], 0.333867)
    close(blocks.loc[(blocks["Dataset"].eq("Grow-out ponds")) & (blocks["Predictor_block"].eq("Water quality")), "Unique_R2"].iloc[0], 0.105258)

    rf = pd.read_csv(taxonomy_results / "rf_cross_validated_performance.csv")
    full_rf = rf.loc[rf["Model"].eq("+ Disturbances")].iloc[0]
    close(full_rf["Mean_R2"], 0.587448)
    close(full_rf["Incremental_mean_R2"], 0.003261)
    close(full_rf["PCoA1_MAE"], 0.072345)
    close(full_rf["PCoA2_MAE"], 0.086557)

    all_models = pd.read_csv(sensitivity_results / "multialgorithm_all_ponds_submitted.csv")
    grow_models = pd.read_csv(sensitivity_results / "multialgorithm_grow_out_only_submitted.csv")
    assert set(all_models["Samples"]) == {160} and set(all_models["Pond_month_groups"]) == {89}
    assert set(grow_models["Samples"]) == {122} and set(grow_models["Pond_month_groups"]) == {70}

    functional = json.loads((functional_results / "functional_analysis_summary.json").read_text())
    assert functional["formal_samples"] == 168 and functional["pond_months"] == 93
    assert functional["predicted_pathways"] == 556
    assert functional["picrust2_qc"]["placed_ASVs_with_domain_assignment"] == 23895
    assert functional["picrust2_qc"]["ASVs_retained_in_metagenome_step"] == 23394
    close(functional["profile_pcoa_positive_axis_variance_percent"]["PCoA1"], 39.0809776013)
    close(functional["profile_pcoa_positive_axis_variance_percent"]["PCoA2"], 12.7871398304)
    close(functional["turnover_association"]["pond_adjusted_rank_correlation"], 0.807762)
    close(functional["turnover_association"]["restricted_two_sided_p"], 0.0001)

    lineage = json.loads((taxonomy_results / "complete_genus_lineage_summary.json").read_text())
    assert lineage["resolved_genus_level_lineages"] == 2673
    assert lineage["ubiquitous_lineages_93_of_93"] == 7
    assert lineage["persistent_lineages_at_least_84_of_93"] == 82
    close(lineage["persistent_combined_mean_relative_abundance_percent"], 53.771634)
    close(lineage["abundance_prevalence_spearman_rho"], 0.859270)

    complete_lineages = pd.read_csv(
        root / "data" / "supporting" / "Additional_file_1_complete_genus_lineages.csv"
    )
    selected = complete_lineages.loc[
        complete_lineages["Selected_top30_seasonal_prevalence_range"]
    ]
    assert len(selected) == 30
    close(
        selected["Grow_out_70_mean_relative_abundance_percent"].sum(),
        2.2854285920,
    )
    abundance_bin = complete_lineages[
        complete_lineages["All_pond_month_mean_relative_abundance_percent"].ge(0.01)
        & complete_lineages["All_pond_month_mean_relative_abundance_percent"].lt(0.03)
    ]
    close(
        abundance_bin["Seasonal_prevalence_range_percentage_points"].median(),
        30.3140096618,
    )

    old_names = [
        primary / "asv_table_168_averaged.csv.gz",
        primary / "asv_table_workflow_204_inputs.csv.gz",
        functional_results / "Table_S13_pathway_module_definitions.csv",
        functional_results / "Table_S14_dependence_aware_functional_tests.csv",
    ]
    assert not any(path.exists() for path in old_names)
    assert (functional_results / "Table_S10_dependence_aware_functional_tests.csv").exists()
    assert (functional_results / "Table_S11_pathway_module_definitions.csv").exists()
    assert (root / "data" / "supporting" / "Additional_file_1_complete_genus_lineages.csv").exists()

    warnings = []
    blank_summary = (
        root
        / "reference_results"
        / "functional_sensitivity"
        / "blank_screened_functional_sensitivity_summary.json"
    )
    assert blank_summary.exists()
    blank_sensitivity = json.loads(blank_summary.read_text())
    assert blank_sensitivity["common_samples"] == 168
    assert blank_sensitivity["union_pathways"] == 556
    assert (
        blank_sensitivity["pairwise_pathway_bray_curtis_correlation"]["n"]
        == 14028
    )
    close(
        blank_sensitivity["pairwise_pathway_bray_curtis_correlation"][
            "spearman_rho"
        ],
        0.9998773041652059,
    )
    assert (
        blank_sensitivity["functional_test_significance"][
            "classification_changes"
        ]
        == 0
    )
    close(
        blank_sensitivity["functional_test_significance"][
            "significance_class_concordance_at_0.05"
        ],
        1.0,
    )
    combined_ec = root / "picrust2_primary" / "combined_EC_predicted.tsv.gz"
    assert combined_ec.exists()
    with gzip.open(combined_ec, "rt", encoding="utf-8") as handle:
        combined_ec_rows = sum(1 for _ in handle) - 1
    assert combined_ec_rows == 23895
    assert (
        root
        / "picrust2_blank_screened"
        / "pathways_out"
        / "path_abun_unstrat.tsv.gz"
    ).exists()
    screened_functional = json.loads(
        (
            root
            / "reference_results"
            / "functional_blank_screened"
            / "results"
            / "functional_analysis_summary.json"
        ).read_text()
    )
    assert screened_functional["picrust2_qc"]["input_ASVs"] == 23787
    assert screened_functional["predicted_pathways"] == 555
    crosswalk_path = root / "data" / "manifests" / "sample_SAMN_SRR_crosswalk.tsv"
    if not crosswalk_path.exists():
        warnings.append("The archived 190-library sample–SAMN–SRR crosswalk is missing.")
    else:
        crosswalk = pd.read_csv(crosswalk_path, sep="\t", dtype=str)
        assert len(crosswalk) == 190
        assert crosswalk["workflow_sample_id"].nunique() == 190
        assert crosswalk["biosample_accession"].nunique() == 190
        assert crosswalk["sra_run_accession"].nunique() == 190
        assert set(crosswalk["workflow_sample_id"]) == set(formal_libraries)
        public_check = root / "validation" / "sra_public_access_verified.txt"
        if not public_check.exists():
            warnings.append(
                "The complete 190-library sample–SAMN–SRR crosswalk is archived, but public "
                "accessibility without login has not yet been verified."
            )

    result = {
        "status": "PASS_WITH_DOCUMENTED_LIMITATIONS" if warnings else "PASS",
        "checked_invariants": {
            "workflow_inputs": 204,
            "workflow_QC_rows": 204,
            "workflow_feature_columns": 199,
            "formal_libraries": 190,
            "formal_libraries_with_SAMN_and_SRR": 190 if crosswalk_path.exists() else 0,
            "canonical_profiles": 168,
            "pond_months": 93,
            "complete_layer_pairs": 75,
            "season_model_samples": 160,
            "season_model_pond_months": 89,
            "grow_out_season_samples": 122,
            "grow_out_season_pond_months": 70,
            "starting_ASVs": 29461,
            "retained_prokaryotic_ASVs": 23896,
            "resolved_genus_lineages": 2673,
            "selected_top30_combined_grow_out_abundance_percent": 2.2854285920,
            "abundance_bin_0.01_to_0.03_median_prevalence_range": 30.3140096618,
            "predicted_MetaCyc_pathways": 556,
            "blank_screened_PICRUSt2_ASVs": 23787,
            "blank_screened_MetaCyc_pathways": 555,
            "blank_screened_pairwise_Bray_Curtis_spearman_rho": 0.9998773041652059,
            "blank_screened_significance_classification_changes": 0,
        },
        "warnings": warnings,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
