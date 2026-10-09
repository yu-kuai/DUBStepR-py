"""Optimal feature set via the density index (port of ``getOptimalFeatureSet.R``)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.linalg
import scipy.sparse as sp
from scipy.spatial import cKDTree

from ._utils import logger


def seurat_scale(x, scale_max=10.0):
    """Seurat ``ScaleData`` (FastSparseRowScale): per-gene z-score with n-1 sd, values
    clipped above at ``scale_max``, NaN set to 0. ``x`` is genes x cells (dense)."""
    n = x.shape[1]
    mean = x.mean(axis=1, keepdims=True)
    sd = np.sqrt(((x - mean) ** 2).sum(axis=1, keepdims=True) / (n - 1))
    with np.errstate(divide="ignore", invalid="ignore"):
        z = (x - mean) / sd
    z = np.minimum(z, scale_max)
    z[np.isnan(z)] = 0.0
    return z


def seurat_pca(scaled, npcs=50):
    """Seurat ``RunPCA`` on scaled data (genes x cells), solved exactly.

    Seurat drops zero-variance features, then takes the top ``min(npcs, n_features - 1)``
    singular triplets of the uncentred cells x features matrix (via irlba).
    Returns cell embeddings (U * d) and stdev (d / sqrt(n_cells - 1)).
    """
    var = scaled.var(axis=1, ddof=1)
    scaled = scaled[var > 0]
    a = scaled.T  # cells x features
    n_cells, n_feat = a.shape
    npcs = min(npcs, n_feat - 1)
    if n_feat <= n_cells:
        evals, evecs = scipy.linalg.eigh(a.T @ a, subset_by_index=[n_feat - npcs, n_feat - 1])
        evals, evecs = evals[::-1], evecs[:, ::-1]
        d = np.sqrt(np.clip(evals, 0, None))
        emb = a @ evecs
    else:
        u, d, _ = np.linalg.svd(a, full_matrices=False)
        d = d[:npcs]
        emb = u[:, :npcs] * d
    stdev = d / np.sqrt(max(1, n_cells - 1))
    return emb, stdev


def density_index(log_feature_data, k=10, num_pcs=20, error=0.0):
    """Mean distance to the k nearest neighbours in PCA space, scaled by the total
    standard deviation of the first ``num_pcs`` PCs."""
    emb, stdev = seurat_pca(seurat_scale(log_feature_data))
    pca = emb[:, : min(num_pcs, emb.shape[1])]
    dists, _ = cKDTree(pca).query(pca, k=k + 1, eps=error)
    mean_nn = dists[:, 1:].mean()  # drop the zero self-distance
    length_scale = np.sqrt(np.sum(stdev[:num_pcs] ** 2))
    return mean_nn / length_scale


def get_optimal_feature_set(filt_data, gene_names, ordered_genes, elbow_pt=25, k=10, num_pcs=20,
                            error=0.0, step=25):
    """Evaluate the density index for the top ``elbow_pt, elbow_pt + 25, ...`` ordered genes
    and return the gene set with the minimum density index.

    Returns
    -------
    optimal_feature_genes : numpy.ndarray
    density_index : pandas.Series (indexed by number of genes)
    """
    logger.info("Determining optimal feature set...")
    row = pd.Index(gene_names).get_indexer(ordered_genes)
    sizes = list(range(elbow_pt, len(ordered_genes) + 1, step))
    di = []
    for num_genes in sizes:
        x = filt_data[row[:num_genes]]
        x = x.toarray() if sp.issparse(x) else np.asarray(x, dtype=np.float64)
        di.append(density_index(x, k=k, num_pcs=num_pcs, error=error))
        logger.debug("num_genes=%d density_index=%.6g", num_genes, di[-1])
    di = pd.Series(di, index=sizes, name="density.index")
    best = sizes[int(np.argmin(di.to_numpy()))]
    logger.info("Done.")
    return np.asarray(ordered_genes)[:best], di
