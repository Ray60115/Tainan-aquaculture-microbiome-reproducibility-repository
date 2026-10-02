#!/usr/bin/env python3
"""Prepare and audit the Tainan 168-sample PICRUSt2 input files.

The biological ASV table contains exactly 168 metadata-matched, averaged U/D
profiles.  Excluded ``-G`` records and the blank-water control are supplied in
a separate non-formal-column table so that a reader never mistakes them for
members of the 168-profile analytical cohort.

Primary PICRUSt2 input:
  * includes only the 168 metadata-matched biological samples;
  * excludes mitochondria, chloroplast, and explicit eukaryotic features;
  * retains bacterial and archaeal features with positive biological abundance.

The blank-water control is not treated as a biological sample. Because there is
only one negative control, blank-enriched ASVs are used for a sensitivity screen
rather than definitive statistical contaminant classification.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from scipy.stats import pearsonr, spearmanr


ANNOTATION_COLUMNS = ("id", "Sequence", "Taxon")
ORGANELLE_PATTERN = re.compile(
    r"mitochond|chloroplast|d__eukary", flags=re.IGNORECASE
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def shannon_rows(counts: np.ndarray) -> np.ndarray:
    totals = counts.sum(axis=1, keepdims=True)
    proportions = np.divide(
        counts,
        totals,
        out=np.zeros_like(counts, dtype=float),
        where=totals > 0,
    )
    logp = np.zeros_like(proportions)
    positive = proportions > 0
    logp[positive] = np.log(proportions[positive])
    return -(proportions * logp).sum(axis=1)


def write_fasta(path: Path, identifiers: pd.Series, sequences: pd.Series) -> None:
    tmp_path = path.with_name(f".{path.name}.tmp")
    with tmp_path.open("w", encoding="utf-8", newline="\n") as handle:
        for identifier, sequence in zip(identifiers, sequences, strict=True):
            handle.write(f">{identifier}\n{sequence}\n")
    tmp_path.replace(path)


def write_abundance(path: Path, identifiers: pd.Series, matrix: pd.DataFrame) -> None:
    output = matrix.copy()
    output.insert(0, "#OTU ID", identifiers.to_numpy())
    tmp_path = path.with_name(f".{path.name}.tmp")
    output.to_csv(tmp_path, sep="\t", index=False, float_format="%.10g")
    tmp_path.replace(path)


def write_dataframe(path: Path, frame: pd.DataFrame, **kwargs) -> None:
    """Write a table atomically so interrupted runs cannot leave partial files."""
    tmp_path = path.with_name(f".{path.name}.tmp")
    frame.to_csv(tmp_path, **kwargs)
    tmp_path.replace(path)


def distance_sensitivity(
    primary_counts: pd.DataFrame,
    screened_counts: pd.DataFrame,
    sample_ids: list[str],
) -> dict[str, float]:
    primary = primary_counts.T.to_numpy(dtype=float)
    screened = screened_counts.T.to_numpy(dtype=float)
    primary /= primary.sum(axis=1, keepdims=True)
    screened /= screened.sum(axis=1, keepdims=True)
    d_primary = pdist(primary, metric="braycurtis")
    d_screened = pdist(screened, metric="braycurtis")
    pearson = pearsonr(d_primary, d_screened)
    spearman = spearmanr(d_primary, d_screened)

    h_primary = shannon_rows(primary_counts.T.to_numpy(dtype=float))
    h_screened = shannon_rows(screened_counts.T.to_numpy(dtype=float))
    delta = h_screened - h_primary
    shannon_table = pd.DataFrame(
        {
            "SampleID": sample_ids,
            "Shannon_primary": h_primary,
            "Shannon_blank_screened": h_screened,
            "Difference_screened_minus_primary": delta,
        }
    )
    return {
        "bray_curtis_pearson_r": float(pearson.statistic),
        "bray_curtis_pearson_p": float(pearson.pvalue),
        "bray_curtis_spearman_rho": float(spearman.statistic),
        "bray_curtis_spearman_p": float(spearman.pvalue),
        "shannon_absolute_difference_median": float(np.median(np.abs(delta))),
        "shannon_absolute_difference_maximum": float(np.max(np.abs(delta))),
    }, shannon_table


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--asv", required=True, type=Path)
    parser.add_argument("--controls", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    outdir = args.outdir.resolve()
    input_dir = outdir / "input"
    qc_dir = outdir / "qc"
    input_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)

    asv = pd.read_csv(args.asv)
    controls = pd.read_csv(args.controls)
    metadata = pd.read_csv(args.metadata)

    missing_annotation = [c for c in ANNOTATION_COLUMNS if c not in asv.columns]
    if missing_annotation:
        raise ValueError(f"Missing ASV annotation columns: {missing_annotation}")
    if "SampleID" not in metadata.columns:
        raise ValueError("Metadata lacks SampleID.")
    if metadata["SampleID"].duplicated().any():
        duplicates = metadata.loc[
            metadata["SampleID"].duplicated(keep=False), "SampleID"
        ].tolist()
        raise ValueError(f"Duplicated metadata SampleID values: {duplicates}")

    sample_ids = metadata["SampleID"].astype(str).tolist()
    if len(sample_ids) != 168:
        raise ValueError(f"Expected 168 formal samples, found {len(sample_ids)}.")
    missing_samples = [sample for sample in sample_ids if sample not in asv.columns]
    if missing_samples:
        raise ValueError(f"ASV table is missing formal samples: {missing_samples}")

    if "Layer" in metadata.columns:
        layers = set(metadata["Layer"].dropna().astype(str))
        if not layers.issubset({"U", "D"}):
            raise ValueError(f"Unexpected formal sample layers: {sorted(layers)}")

    nonformal_columns = [
        c for c in controls.columns if c not in ANNOTATION_COLUMNS
    ]
    g_columns = [c for c in nonformal_columns if c.endswith("-G")]
    if "Water" not in controls.columns:
        raise ValueError("Expected blank-water control column 'Water' was not found.")
    unexpected_biological_columns = [
        c for c in asv.columns if c not in ANNOTATION_COLUMNS and c not in sample_ids
    ]
    if unexpected_biological_columns:
        raise ValueError(
            "The 168-profile ASV table contains non-formal columns: "
            f"{unexpected_biological_columns}"
        )
    if set(asv["id"].astype(str)) != set(controls["id"].astype(str)):
        raise ValueError("Biological and non-formal ASV identifier sets differ.")

    counts_all = asv[sample_ids].fillna(0).astype(float)
    blank = (
        controls.set_index("id")["Water"]
        .reindex(asv["id"].astype(str))
        .fillna(0)
        .astype(float)
        .reset_index(drop=True)
    )
    if (counts_all < 0).any().any() or (blank < 0).any():
        raise ValueError("Negative abundance values were detected.")
    if (counts_all.sum(axis=0) <= 0).any():
        empty = counts_all.columns[counts_all.sum(axis=0) <= 0].tolist()
        raise ValueError(f"Formal samples with zero library size: {empty}")

    sequences = asv["Sequence"].astype(str).str.upper()
    invalid_sequence = ~sequences.str.fullmatch(r"[ACGTN]+")
    duplicate_ids = asv["id"].astype(str).duplicated(keep=False)
    duplicate_sequences = sequences.duplicated(keep=False)

    taxonomy = asv["Taxon"].fillna("").astype(str)
    organelle = taxonomy.str.contains(ORGANELLE_PATTERN)
    positive_biological = counts_all.sum(axis=1) > 0
    primary_keep = positive_biological & ~organelle & ~invalid_sequence

    biological_library = counts_all.sum(axis=0)
    biological_relative = counts_all.div(biological_library, axis=1)
    blank_relative = blank / blank.sum()
    biological_max_relative = biological_relative.max(axis=1)
    biological_mean_relative = biological_relative.mean(axis=1)
    biological_prevalence = (counts_all > 0).mean(axis=1)

    # A deliberately conservative sensitivity flag. With only one blank, this
    # is not treated as definitive contaminant classification.
    blank_enriched = (
        (blank > 0)
        & positive_biological
        & (blank_relative > biological_max_relative)
    )
    blank_screened_keep = primary_keep & ~blank_enriched

    primary = asv.loc[primary_keep].copy()
    primary_counts = primary[sample_ids].astype(float)
    screened = asv.loc[blank_screened_keep].copy()
    screened_counts = screened[sample_ids].astype(float)

    primary_fasta = input_dir / "study_seqs_168_prokaryotic.fna"
    primary_tsv = input_dir / "study_abundance_168_prokaryotic.tsv"
    screened_fasta = input_dir / "study_seqs_168_prokaryotic_blank_screened.fna"
    screened_tsv = input_dir / "study_abundance_168_prokaryotic_blank_screened.tsv"
    write_fasta(primary_fasta, primary["id"].astype(str), primary["Sequence"].astype(str))
    write_abundance(primary_tsv, primary["id"].astype(str), primary_counts)
    write_fasta(
        screened_fasta,
        screened["id"].astype(str),
        screened["Sequence"].astype(str),
    )
    write_abundance(
        screened_tsv,
        screened["id"].astype(str),
        screened_counts,
    )

    organelle_sample_fraction = counts_all.loc[organelle].sum(axis=0) / biological_library
    shared_blank_fraction = biological_relative.loc[blank > 0].sum(axis=0)
    blank_enriched_fraction = biological_relative.loc[blank_enriched].sum(axis=0)

    audit = pd.DataFrame(
        {
            "ASV_ID": asv["id"].astype(str),
            "Taxon": taxonomy,
            "Sequence_length": sequences.str.len(),
            "Biological_total_abundance": counts_all.sum(axis=1),
            "Biological_prevalence": biological_prevalence,
            "Biological_mean_relative_abundance": biological_mean_relative,
            "Biological_max_relative_abundance": biological_max_relative,
            "Blank_abundance": blank,
            "Blank_relative_abundance": blank_relative,
            "Detected_in_blank": blank > 0,
            "Blank_enriched_sensitivity_flag": blank_enriched,
            "Organelle_or_eukaryotic": organelle,
            "Invalid_sequence_characters": invalid_sequence,
            "Included_in_primary_PICRUSt2_input": primary_keep,
            "Included_in_blank_screened_input": blank_screened_keep,
        }
    )
    write_dataframe(qc_dir / "asv_and_blank_audit.csv", audit, index=False)

    sample_qc = metadata.copy()
    sample_qc["Original_library_size"] = biological_library.to_numpy()
    sample_qc["Prokaryotic_library_size"] = primary_counts.sum(axis=0).to_numpy()
    sample_qc["Prokaryotic_retained_fraction"] = (
        sample_qc["Prokaryotic_library_size"] / sample_qc["Original_library_size"]
    )
    sample_qc["Fraction_from_ASVs_detected_in_blank"] = (
        shared_blank_fraction.to_numpy()
    )
    sample_qc["Fraction_from_blank_enriched_ASVs"] = (
        blank_enriched_fraction.to_numpy()
    )
    write_dataframe(qc_dir / "sample_input_qc.csv", sample_qc, index=False)

    sensitivity_summary, shannon_table = distance_sensitivity(
        primary_counts,
        screened_counts,
        sample_ids,
    )
    write_dataframe(
        qc_dir / "blank_screening_shannon_sensitivity.csv",
        shannon_table,
        index=False,
    )

    summary = {
        "formal_sample_count": len(sample_ids),
        "formal_layers": sorted(set(metadata["Layer"].astype(str)))
        if "Layer" in metadata.columns
        else None,
        "nonformal_columns_excluded": nonformal_columns,
        "g_columns_excluded": g_columns,
        "blank_column_excluded": "Water",
        "blank_library_size": float(blank.sum()),
        "total_ASV_records": int(len(asv)),
        "positive_biological_ASVs": int(positive_biological.sum()),
        "blank_positive_ASVs": int((blank > 0).sum()),
        "blank_only_ASVs": int(((blank > 0) & ~positive_biological).sum()),
        "shared_blank_biological_ASVs": int(((blank > 0) & positive_biological).sum()),
        "blank_enriched_sensitivity_ASVs": int(blank_enriched.sum()),
        "organelle_or_eukaryotic_ASVs_positive_in_biological_samples": int(
            (organelle & positive_biological).sum()
        ),
        "organelle_or_eukaryotic_read_fraction_overall": float(
            counts_all.loc[organelle].to_numpy().sum()
            / counts_all.to_numpy().sum()
        ),
        "primary_PICRUSt2_ASV_count": int(primary_keep.sum()),
        "blank_screened_PICRUSt2_ASV_count": int(blank_screened_keep.sum()),
        "biological_library_size": {
            "minimum": float(biological_library.min()),
            "median": float(biological_library.median()),
            "mean": float(biological_library.mean()),
            "maximum": float(biological_library.max()),
        },
        "prokaryotic_retained_fraction": {
            "minimum": float(sample_qc["Prokaryotic_retained_fraction"].min()),
            "median": float(sample_qc["Prokaryotic_retained_fraction"].median()),
            "mean": float(sample_qc["Prokaryotic_retained_fraction"].mean()),
            "maximum": float(sample_qc["Prokaryotic_retained_fraction"].max()),
        },
        "blank_enriched_ASV_read_fraction_in_biological_samples": {
            "minimum": float(blank_enriched_fraction.min()),
            "median": float(blank_enriched_fraction.median()),
            "mean": float(blank_enriched_fraction.mean()),
            "maximum": float(blank_enriched_fraction.max()),
        },
        "duplicate_ASV_ID_records": int(duplicate_ids.sum()),
        "duplicate_sequence_records": int(duplicate_sequences.sum()),
        "invalid_sequence_records": int(invalid_sequence.sum()),
        "sequence_length": {
            "minimum": int(sequences.str.len().min()),
            "median": float(sequences.str.len().median()),
            "maximum": int(sequences.str.len().max()),
        },
        "blank_screening_sensitivity": sensitivity_summary,
        "input_sha256": {
            args.asv.name: sha256_file(args.asv),
            args.controls.name: sha256_file(args.controls),
            primary_fasta.name: sha256_file(primary_fasta),
            primary_tsv.name: sha256_file(primary_tsv),
            screened_fasta.name: sha256_file(screened_fasta),
            screened_tsv.name: sha256_file(screened_tsv),
        },
    }
    summary_path = qc_dir / "input_and_blank_qc_summary.json"
    summary_tmp = summary_path.with_name(f".{summary_path.name}.tmp")
    with summary_tmp.open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)
    summary_tmp.replace(summary_path)

    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
