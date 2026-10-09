"""Gene-gene correlation and correlation-range scoring (port of ``getGGC.R`` and
``getCorrelationRange.R``)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ._utils import logger, stable_order_desc


def cor_sparse(data):
    """Pearson correlation between rows of a (genes x cells) matrix.

    Same formula as ``qlcMatrix::corSparse(t(data))``:
    ``cov = (X'X - n * mu mu') / (n - 1)``, ``cor = cov / (sd sd')``.
    """
    X = data.T.tocsc() if sp.issparse(data) else np.asarray(data, dtype=np.float64).T
    n = X.shape[0]
    mu = np.asarray(X.mean(axis=0)).ravel()
    xtx = X.T @ X
    xtx = xtx.toarray() if sp.issparse(xtx) else np.asarray(xtx)
    cov = (xtx - n * np.outer(mu, mu)) / (n - 1)
    sd = np.sqrt(np.diag(cov))
    with np.errstate(divide="ignore", invalid="ignore"):
        return cov / np.outer(sd, sd)


def get_correlation_range(corr):
    """Correlation range per gene: 3rd-largest correlation - 0.75 * smallest correlation.

    (The largest is the gene's self-correlation of 1.) Returns a Series sorted in
    decreasing order, ties kept in input order (as R's ``sort(decreasing = TRUE)``).
    """
    n = corr.shape[0]
    c = np.where(np.isnan(corr), -np.inf, corr)
    max3 = np.partition(c, n - 3, axis=0)[n - 3]
    min_corr = np.where(np.isnan(corr), np.inf, corr).min(axis=0)
    return max3 - 0.75 * min_corr


def r_cut_codes(x, breaks: int):
    """Interval codes from R's ``cut(x, breaks = <int>)`` (right-closed, 0-based)."""
    x = np.asarray(x, dtype=np.float64)
    nb = int(breaks) + 1
    lo, hi = float(np.min(x)), float(np.max(x))
    dx = hi - lo
    if dx == 0:
        dx = abs(lo) if lo != 0 else 1.0
        edges = np.linspace(lo - dx / 1000, hi + dx / 1000, nb)
    else:
        edges = np.linspace(lo, hi, nb)
        edges[0], edges[-1] = lo - dx / 1000, hi + dx / 1000
    # (edges[i], edges[i+1]] -> i
    return np.searchsorted(edges, x, side="left") - 1


def get_ggc(data, gene_names, z_threshold=0.7, align_bins=False):
    """Compute the GGC and z-scored correlation range of each gene.

    Parameters
    ----------
    data
        Filtered log-normalised (genes x cells) matrix.
    gene_names
        Row names of ``data``.
    z_threshold
        Genes with z-transformed log correlation range above this are kept in the GGC.
    align_bins
        The R implementation passes the correlation-range vector (sorted by range) to
        ``tapply`` with the expression bins in the original gene order, so each gene is
        z-scored within the bin of whichever gene originally sat at its sorted position.
        ``False`` (default) reproduces that behaviour so results match R DUBStepR;
        ``True`` assigns each gene to its own mean-expression bin.

    Returns
    -------
    corr_range : pandas.Series
        z-transformed log correlation range, in R's order (by bin, then decreasing range).
    ggc : numpy.ndarray
        Correlation matrix restricted to genes with ``corr_range > z_threshold``.
    ggc_genes : numpy.ndarray
        Names of the GGC rows/columns.
    """
    logger.info("Computing GGC...")
    gene_names = np.asarray(gene_names)
    n_genes = data.shape[0]
    num_bins = min(20, n_genes - 1)
    gene_mean = np.asarray(data.mean(axis=1)).ravel()
    gene_bins = r_cut_codes(gene_mean, num_bins)

    corr = cor_sparse(data)

    crange = get_correlation_range(corr)
    order = stable_order_desc(crange)
    log_range = np.log(1 + crange[order])
    ranked_genes = order  # indices into gene_names, sorted by range

    bins = gene_bins[order] if align_bins else gene_bins[: len(order)]

    z = np.empty_like(log_range)
    out_idx = []
    all_na = True
    for b in np.unique(bins):  # factor levels in increasing order; empty bins drop out
        pos = np.flatnonzero(bins == b)
        vals = log_range[pos]
        if len(vals) > 1:
            zv = (vals - vals.mean()) / vals.std(ddof=1)
            all_na = False
        else:
            zv = np.array([np.nan])
        zv[np.isnan(zv)] = 0.0
        z[pos] = zv
        out_idx.append(pos)
    if all_na:
        raise RuntimeError("Feature correlation range could not be obtained.")
    out_idx = np.concatenate(out_idx)
    corr_range = pd.Series(z[out_idx], index=gene_names[ranked_genes[out_idx]], name="corr.range")

    top = ranked_genes[out_idx][corr_range.to_numpy() > z_threshold]
    ggc = corr[np.ix_(top, top)]
    logger.info("Done.")
    return corr_range, ggc, gene_names[top]
