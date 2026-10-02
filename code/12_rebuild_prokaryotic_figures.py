#!/usr/bin/env python3
"""Rebuild manuscript Figures 1b–f, 2, 3, 5, 7, 8, and 9."""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

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
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-tainan-prok")
    sys.path.insert(0, str(args.helper_dir.resolve()))
    import rebuild_dependence_aware_figures as rebuild

    # The archived helper resolves its original input directory relative to
    # itself. Point the climate input to the recovered package explicitly so
    # figure regeneration remains portable after extraction.
    rebuild.CLIMATE_FILE = args.asv.resolve().parent / "climate_daily_complete_2024.csv"

    def load_filtered_community():
        meta = pd.read_csv(
            args.results / "tables" / "validated_sample_metadata_with_shannon.csv"
        )
        asv = pd.read_csv(args.asv)
        samples = meta["SampleID"].tolist()
        taxonomy = asv["Taxon"].fillna("").astype(str)
        positive = asv[samples].fillna(0).sum(axis=1) > 0
        organelle = taxonomy.str.contains(ORGANELLE_PATTERN)
        kept = asv.loc[positive & ~organelle]
        counts = kept[samples].T
        counts.index = samples
        relative = rebuild.relative_abundance_rows(counts.to_numpy(dtype=float))
        distance = squareform(pdist(relative, metric="braycurtis"))
        return meta, relative, distance

    rebuild.load_community = load_filtered_community
    rebuild.rebuild_figure2()
    rebuild.rebuild_figure3()
    rebuild.rebuild_figure4()
    # The current Fig. 4 prevalence panel is rebuilt by
    # 14_plot_current_seasonal_prevalence.py.  The archived abundance-based
    # seasonal panel is intentionally not emitted by the integrated workflow.
    rebuild.rebuild_figure6()
    rebuild.rebuild_figure9()
    rebuild.rebuild_figure10()
    rebuild.rebuild_environmental_vectors()


if __name__ == "__main__":
    main()
