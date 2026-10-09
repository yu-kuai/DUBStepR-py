"""Input handling and small helpers shared across the package."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
import scipy.sparse as sp

logger = logging.getLogger("dubstepr")


def as_gene_matrix(data, gene_names=None, cell_names=None):
    """Coerce input to a (genes x cells) CSR float64 matrix plus gene/cell names.

    Accepts a pandas DataFrame (genes as index, cells as columns), a numpy
    array or a scipy sparse matrix. For array/sparse input, ``gene_names``
    must be supplied; ``cell_names`` is optional.
    """
    if isinstance(data, pd.DataFrame):
        if gene_names is None:
            gene_names = data.index
        if cell_names is None:
            cell_names = data.columns
        mat = sp.csr_matrix(data.to_numpy(dtype=np.float64))
    elif sp.issparse(data):
        mat = sp.csr_matrix(data, dtype=np.float64)
    else:
        mat = sp.csr_matrix(np.asarray(data, dtype=np.float64))

    if gene_names is None:
        raise ValueError("gene_names must be provided when data is not a pandas DataFrame")
    gene_names = np.asarray(gene_names, dtype=object).astype(str)
    if gene_names.shape[0] != mat.shape[0]:
        raise ValueError(
            f"len(gene_names)={gene_names.shape[0]} does not match number of rows {mat.shape[0]} "
            "(input must be genes x cells)"
        )
    if pd.Index(gene_names).has_duplicates:
        raise ValueError("gene names must be unique")
    if cell_names is None:
        cell_names = np.array([f"cell{i + 1}" for i in range(mat.shape[1])], dtype=object)
    cell_names = np.asarray(cell_names, dtype=object).astype(str)
    mat.sum_duplicates()
    mat.sort_indices()
    return mat, gene_names, cell_names


def stable_order_desc(x):
    """Indices sorting ``x`` in decreasing order, ties kept in input order, NaNs dropped.

    Mirrors R's ``sort(x, decreasing = TRUE)`` (radix order, stable, NA removed).
    """
    x = np.asarray(x, dtype=np.float64)
    idx = np.flatnonzero(~np.isnan(x))
    return idx[np.argsort(-x[idx], kind="stable")]
