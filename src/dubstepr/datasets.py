"""Bundled example data."""

from importlib import resources

import pandas as pd


def pbmc_norm_small_data() -> pd.DataFrame:
    """Small PBMC dataset (230 genes x 80 cells), normalised and log-transformed data from
    Seurat's ``pbmc_small`` object, as shipped with R DUBStepR."""
    path = resources.files("dubstepr") / "data" / "pbmc_norm_small_data.tsv.gz"
    with resources.as_file(path) as p:
        return pd.read_csv(p, sep="\t", index_col=0, float_precision="round_trip")
