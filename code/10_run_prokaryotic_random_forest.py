#!/usr/bin/env python3
"""Rerun grouped random-forest reconstruction on organelle-filtered profiles."""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform


ORGANELLE_PATTERN = re.compile(
    r"mitochond|chloroplast|d__eukary", flags=re.IGNORECASE
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--helper-dir", required=True, type=Path)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--asv", required=True, type=Path)
    args = parser.parse_args()

    os.environ["TAINAN_RESULTS_DIR"] = str(args.results.resolve())
    sys.path.insert(0, str(args.helper_dir.resolve()))
    import run_random_forest_reconstruction as rf

    def load_filtered_data() -> pd.DataFrame:
        meta = pd.read_csv(
            args.results / "tables" / "validated_sample_metadata_with_shannon.csv"
        )
        asv = pd.read_csv(args.asv)
        samples = meta["SampleID"].tolist()
        taxonomy = asv["Taxon"].fillna("").astype(str)
        positive = asv[samples].fillna(0).sum(axis=1) > 0
        organelle = taxonomy.str.contains(ORGANELLE_PATTERN)
        kept = asv.loc[positive & ~organelle]
        counts = kept[samples].T.to_numpy(dtype=float)
        relative = rf.relative_abundance_rows(counts)
        distance = squareform(pdist(relative, metric="braycurtis"))
        scores, _, _ = rf.pcoa(distance)
        meta["PCoA1"] = scores[:, 0]
        meta["PCoA2"] = scores[:, 1]
        work = meta.loc[
            meta["DisturbanceComplete"].astype(bool)
            & meta["Season_verified"].notna()
        ].copy()
        work = work.reset_index(drop=True)
        month = work["Date"].astype(str).str[4:].astype(int)
        work["Month_sin"] = np.sin(2 * np.pi * month / 12)
        work["Month_cos"] = np.cos(2 * np.pi * month / 12)
        work["Year"] = work["Date"].astype(str).str[:4].astype(int)
        return work

    rf.load_analysis_data = load_filtered_data
    # The helper's two legacy figure files predate the submitted four-panel
    # layout.  Results are written here; 11_make_merged_rf_figure.py builds the
    # current manuscript Fig. 11 from those exact tables.
    rf.make_figures = lambda performance, predictions, importance: None
    rf.main()


if __name__ == "__main__":
    main()
