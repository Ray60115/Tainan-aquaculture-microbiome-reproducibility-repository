#!/usr/bin/env python3
"""Compare PICRUSt2 best-domain assignments with supplied ASV taxonomy."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def normalize_picrust_domain(value: object) -> str:
    text = str(value).strip().casefold()
    if text in {"bac", "bacteria", "bacterial"}:
        return "Bacteria"
    if text in {"arc", "archaea", "archaeal"}:
        return "Archaea"
    return str(value)


def supplied_domain(taxon: object) -> str:
    text = str(taxon).strip()
    if text.startswith("d__Bacteria"):
        return "Bacteria"
    if text.startswith("d__Archaea"):
        return "Archaea"
    return "Other/unspecified"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--marker", required=True, type=Path)
    parser.add_argument("--taxonomy-audit", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    marker = pd.read_csv(args.marker, sep="\t", compression="infer", index_col=0)
    taxonomy = pd.read_csv(args.taxonomy_audit)
    taxonomy = taxonomy[
        taxonomy["Included_in_primary_PICRUSt2_input"].astype(bool)
    ].copy()
    taxonomy["Supplied_domain"] = taxonomy["Taxon"].map(supplied_domain)

    if "best_domain" not in marker.columns:
        raise ValueError("PICRUSt2 marker table lacks best_domain.")
    marker_domain = pd.DataFrame(
        {
            "PICRUSt2_best_domain": marker["best_domain"].map(
                normalize_picrust_domain
            ),
            "PICRUSt2_NSTI": pd.to_numeric(
                marker["metadata_NSTI"], errors="coerce"
            ),
        }
    )
    marker_domain.index = marker_domain.index.astype(str)
    merged = taxonomy.merge(
        marker_domain,
        left_on="ASV_ID",
        right_index=True,
        how="left",
        validate="one_to_one",
    )
    merged["Domain_assignment_available"] = merged["PICRUSt2_best_domain"].notna()
    merged["Domain_concordant"] = (
        merged["Supplied_domain"].eq(merged["PICRUSt2_best_domain"])
        & merged["Domain_assignment_available"]
    )
    merged["Retained_under_max_NSTI_2"] = (
        merged["Domain_assignment_available"]
        & merged["PICRUSt2_NSTI"].le(2)
    )

    assigned = merged[merged["Domain_assignment_available"]]
    confusion = pd.crosstab(
        assigned["Supplied_domain"],
        assigned["PICRUSt2_best_domain"],
        margins=True,
    )
    assigned_reads = float(assigned["Biological_total_abundance"].sum())
    concordant_reads = float(
        assigned.loc[assigned["Domain_concordant"], "Biological_total_abundance"].sum()
    )
    nsti_retained = assigned[assigned["Retained_under_max_NSTI_2"]]
    nsti_retained_reads = float(nsti_retained["Biological_total_abundance"].sum())
    nsti_retained_concordant_reads = float(
        nsti_retained.loc[
            nsti_retained["Domain_concordant"], "Biological_total_abundance"
        ].sum()
    )
    summary = {
        "input_ASVs": int(len(merged)),
        "ASVs_with_PICRUSt2_domain_assignment": int(len(assigned)),
        "ASVs_without_domain_assignment": int((~merged["Domain_assignment_available"]).sum()),
        "ASV_level_concordance_fraction": float(assigned["Domain_concordant"].mean()),
        "read_weighted_concordance_fraction": (
            concordant_reads / assigned_reads if assigned_reads > 0 else None
        ),
        "ASVs_retained_under_max_NSTI_2": int(len(nsti_retained)),
        "retained_ASV_level_concordance_fraction": float(
            nsti_retained["Domain_concordant"].mean()
        ),
        "retained_read_weighted_concordance_fraction": (
            nsti_retained_concordant_reads / nsti_retained_reads
            if nsti_retained_reads > 0
            else None
        ),
        "retained_domain_mismatch_reads": float(
            nsti_retained_reads - nsti_retained_concordant_reads
        ),
        "supplied_domain_counts": merged["Supplied_domain"].value_counts().to_dict(),
        "picrust2_best_domain_counts": assigned[
            "PICRUSt2_best_domain"
        ].value_counts().to_dict(),
        "interpretation": (
            "This is a QC comparison between two reference-based assignments, "
            "not an independent truth standard. Domain mismatches with high "
            "NSTI are excluded from metagenome prediction by the prespecified "
            "maximum NSTI of 2."
        ),
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    merged.to_csv(
        args.output_dir / "picrust2_domain_assignment_audit.csv",
        index=False,
    )
    confusion.to_csv(args.output_dir / "picrust2_domain_assignment_confusion.csv")
    (args.output_dir / "picrust2_domain_assignment_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
