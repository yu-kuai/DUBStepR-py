"""Gene filtering (port of ``getFilteredData.R``)."""

from __future__ import annotations

import re
from functools import lru_cache
from importlib import resources

import numpy as np
import pandas as pd

from ._utils import as_gene_matrix, logger

SPECIES = ("human", "mouse", "rat")


@lru_cache(maxsize=None)
def load_gene_annotation(species: str = "human") -> pd.DataFrame:
    """Return the Ensembl gene table shipped with DUBStepR (``allgene_list[[species]]``).

    Columns: ``Gene.stable.ID``, ``Gene.name``, ``Gene.Synonym``, ``Pseudogene.Status``.
    """
    if species not in SPECIES:
        raise ValueError(f"species must be one of {SPECIES}, got {species!r}")
    path = resources.files("dubsteppy") / "data" / f"allgene_{species}.tsv.gz"
    with resources.as_file(path) as p:
        return pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False, na_filter=False)


def _identifiers(table: pd.DataFrame) -> set:
    ids = set(table["Gene.stable.ID"]) | set(table["Gene.name"]) | set(table["Gene.Synonym"])
    ids.discard("")  # missing values in the R table; never valid gene names
    return ids


@lru_cache(maxsize=None)
def _excluded_identifiers(species: str):
    table = load_gene_annotation(species)
    mito = table[table["Gene.name"].str.contains(r"^MT-", flags=re.IGNORECASE, regex=True)]
    ribo = table[table["Gene.name"].str.contains(r"^RPS|^RPL", flags=re.IGNORECASE, regex=True)]
    pseudo = table[table["Pseudogene.Status"] == "Pseudogene"]
    return _identifiers(mito) | _identifiers(ribo), _identifiers(pseudo)


def filter_genes_mask(data, gene_names, min_cells=None, species="human"):
    """Boolean mask over genes (rows) kept by DUBStepR's filtering step."""
    n_cells = data.shape[1]
    if min_cells is None:
        min_cells = 0.05 * n_cells
    logger.info("Dimensions of input data: %d x %d", data.shape[0], n_cells)

    # Expression filter: expressed (> 0) in strictly more than min_cells cells
    n_expr = np.asarray((data > 0).sum(axis=1)).ravel()
    keep = n_expr > min_cells
    logger.info("Expression Filtering Done.")

    names = pd.Series(gene_names)
    # Spike-ins
    keep &= ~names.str.contains(r"^ERCC-", flags=re.IGNORECASE, regex=True).to_numpy()
    # Mitochondrial, ribosomal and pseudogenes (matched on Ensembl ID, symbol or synonym)
    mito_ribo, pseudo = _excluded_identifiers(species)
    keep &= ~names.isin(mito_ribo | pseudo).to_numpy()
    logger.info("Mitochondrial, Ribosomal and Pseudo Genes Filtering Done.")
    logger.info("Dimensions of filtered data: %d x %d", int(keep.sum()), n_cells)
    return keep


def get_filtered_data(data, gene_names=None, min_cells=None, species="human"):
    """Filter lowly expressed, spike-in, mitochondrial, ribosomal and pseudo genes.

    Parameters
    ----------
    data
        Genes x cells expression matrix (DataFrame, ndarray or scipy sparse).
    gene_names
        Gene names (required unless ``data`` is a DataFrame).
    min_cells
        Genes must be expressed in more than this many cells. Default ``0.05 * n_cells``.
    species
        ``"human"``, ``"mouse"`` or ``"rat"``.

    Returns
    -------
    (filtered CSR matrix, kept gene names)
    """
    mat, gene_names, _ = as_gene_matrix(data, gene_names)
    keep = filter_genes_mask(mat, gene_names, min_cells, species)
    return mat[keep], gene_names[keep]
