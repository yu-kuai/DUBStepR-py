"""Main DUBStepR entry points (port of ``DUBStepR.R`` and ``logNormalize.R``)."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ._utils import as_gene_matrix, logger
from .filtering import filter_genes_mask
from .ggc import get_ggc
from .optimal import get_optimal_feature_set
from .stepwise import run_stepwise_reg


@dataclass
class DUBStepRResult:
    """Output of :func:`dubstepr`.

    Attributes
    ----------
    corr_info
        DataFrame with columns ``feature.genes`` (all GGC genes in DUBStepR order) and
        ``corr.range`` (their z-transformed correlation range).
    elbow_pt
        Number of genes selected by stepwise regression (elbow of the scree plot).
    optimal_feature_genes
        Optimal feature set (``None`` if ``optimise_features=False``).
    density_index
        Density index per candidate feature-set size (``None`` if not optimised).
    corr_range
        z-transformed correlation range of every filtered gene.
    scree
        Scree values from the stepwise regression.
    """

    corr_info: pd.DataFrame
    elbow_pt: int
    optimal_feature_genes: Optional[np.ndarray] = None
    density_index: Optional[pd.Series] = None
    corr_range: Optional[pd.Series] = field(default=None, repr=False)
    scree: Optional[np.ndarray] = field(default=None, repr=False)

    @property
    def feature_genes(self) -> np.ndarray:
        """All correlated genes in DUBStepR order."""
        return self.corr_info["feature.genes"].to_numpy()


def dubstepr(data, gene_names=None, min_cells=None, species="human", optimise_features=True,
             k=10, num_pcs=20, error=0.0, align_bins=False):
    """Select feature genes for clustering with DUBStepR.

    Parameters
    ----------
    data
        Log-normalised expression matrix, **genes x cells** (pandas DataFrame with genes as
        index, numpy array, or scipy sparse matrix). For AnnData use :func:`dubstepr_anndata`.
    gene_names
        Gene names; required unless ``data`` is a DataFrame.
    min_cells
        Keep genes expressed in more than this many cells. Default ``0.05 * n_cells``.
    species
        ``"human"`` (default), ``"mouse"`` or ``"rat"``; used to remove mitochondrial,
        ribosomal and pseudo genes.
    optimise_features
        Determine the optimal feature set using the density index.
    k
        Number of nearest neighbours for the density index.
    num_pcs
        Number of principal components for the density index.
    error
        Approximation tolerance for the kNN search (set to 1 for > 10000 cells, as in R).
    align_bins
        See :func:`dubstepr.ggc.get_ggc`. ``False`` reproduces R DUBStepR exactly.

    Returns
    -------
    DUBStepRResult
    """
    logger.info("Running DUBStepR...")
    mat, gene_names, _ = as_gene_matrix(data, gene_names)
    keep = filter_genes_mask(mat, gene_names, min_cells=min_cells, species=species)
    filt, filt_genes = mat[keep], gene_names[keep]

    corr_range, ggc, ggc_genes = get_ggc(filt, filt_genes, align_bins=align_bins)
    feature_genes, elbow_pt, scree = run_stepwise_reg(ggc, ggc_genes)
    corr_info = pd.DataFrame({
        "feature.genes": feature_genes,
        "corr.range": corr_range.loc[feature_genes].to_numpy(),
    })

    result = DUBStepRResult(corr_info=corr_info, elbow_pt=elbow_pt, corr_range=corr_range, scree=scree)
    if optimise_features:
        if filt.shape[1] > 10000:
            error = 1.0
        result.optimal_feature_genes, result.density_index = get_optimal_feature_set(
            filt, filt_genes, feature_genes, elbow_pt=elbow_pt, k=k, num_pcs=num_pcs, error=error)
    return result


def dubstepr_anndata(adata, layer=None, key_added="dubstepr", **kwargs):
    """Run :func:`dubstepr` on an AnnData object (cells x genes, log-normalised).

    Writes ``adata.var[f"{key_added}_feature"]`` (optimal feature set, or all correlated
    genes if ``optimise_features=False``), ``adata.var[f"{key_added}_rank"]`` (DUBStepR
    order, NaN for genes not ranked) and ``adata.uns[key_added]``. Returns the result.
    """
    x = adata.layers[layer] if layer is not None else adata.X
    x = x.T.tocsr() if sp.issparse(x) else np.asarray(x).T
    res = dubstepr(x, gene_names=adata.var_names.to_numpy(), **kwargs)
    chosen = res.optimal_feature_genes if res.optimal_feature_genes is not None else res.feature_genes
    adata.var[f"{key_added}_feature"] = adata.var_names.isin(chosen)
    rank = pd.Series(np.arange(1, len(res.feature_genes) + 1, dtype=float), index=res.feature_genes)
    adata.var[f"{key_added}_rank"] = rank.reindex(adata.var_names).to_numpy()
    adata.uns[key_added] = {
        "elbow_pt": res.elbow_pt,
        "feature_genes": res.feature_genes,
        "optimal_feature_genes": res.optimal_feature_genes,
        "density_index": None if res.density_index is None else res.density_index.to_dict(),
    }
    return res


def log_normalize(raw_data, scale_factor=10000):
    """``log(1 + counts / (library_size / scale_factor))`` per cell; genes x cells input."""
    if isinstance(raw_data, pd.DataFrame):
        depth = raw_data.sum(axis=0) / scale_factor
        return np.log1p(raw_data / depth)
    if sp.issparse(raw_data):
        x = sp.csc_matrix(raw_data, dtype=np.float64, copy=True)
        depth = np.asarray(x.sum(axis=0)).ravel() / scale_factor
        x = x @ sp.diags(1 / depth)
        x.data = np.log1p(x.data)
        return x.tocsr()
    x = np.asarray(raw_data, dtype=np.float64)
    return np.log1p(x / (x.sum(axis=0) / scale_factor))
