#!/usr/bin/env python3
"""Rebuild Figure 3 from the finalized prokaryotic ASV cohort."""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path

import pandas as pd
from scipy.spatial.distance import pdist, squareform


CODE = Path(__file__).resolve().parent
ROOT = CODE.parent
HELPER_DIR = CODE
RESULTS = ROOT / "reproduced_results" / "taxonomy"
ASV_FILE = ROOT / "data" / "primary" / "asv_table_168_profiles.csv.gz"
ORGANELLE_PATTERN = re.compile(
    r"mitochond|chloroplast|d__eukary", flags=re.IGNORECASE
)


def main() -> None:
    os.environ["TAINAN_RESULTS_DIR"] = str(RESULTS)
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-tainan-fig3-repair")
    sys.path.insert(0, str(HELPER_DIR))

    import rebuild_dependence_aware_figures as rebuild

    def load_filtered_community():
        meta = pd.read_csv(
            RESULTS / "tables" / "validated_sample_metadata_with_shannon.csv"
        )
        asv = pd.read_csv(ASV_FILE)
        samples = meta["SampleID"].tolist()
        taxonomy = asv["Taxon"].fillna("").astype(str)
        positive = asv[samples].fillna(0).sum(axis=1) > 0
        organelle = taxonomy.str.contains(ORGANELLE_PATTERN)
        kept = asv.loc[positive & ~organelle]
        counts = kept[samples].T
        relative = rebuild.relative_abundance_rows(counts.to_numpy(dtype=float))
        distance = squareform(pdist(relative, metric="braycurtis"))
        return meta, relative, distance

    rebuild.load_community = load_filtered_community
    rebuild.rebuild_figure3()


if __name__ == "__main__":
    main()
