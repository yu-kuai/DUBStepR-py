"""Stepwise regression on the GGC (port of ``runStepwiseReg.R`` and ``findElbow.R``)."""

from __future__ import annotations

import warnings

import numpy as np

from ._utils import logger, stable_order_desc


def find_elbow(y, return_index=True):
    """Index (1-based, as in R) of the point furthest from the line joining the first and
    last points of ``y``.

    Port of ``findElbow`` (B. A. Hanson), which measures each point's distance to a
    20-unit segment of that line centred on the point's x.
    """
    y = np.asarray(y, dtype=np.float64)
    n = len(y)
    x = np.arange(1, n + 1, dtype=np.float64)
    m = (y[-1] - y[0]) / (x[-1] - x[0])
    b = y[0] - m * x[0]

    use = slice(1, n - 1)
    refpts = m * x[use] + b
    if not (np.all(refpts > y[use]) or np.all(refpts < y[use])):
        warnings.warn("Your curve doesn't appear to be concave", stacklevel=2)

    px, py = x[use], y[use]
    x1, x2 = px - 10, px + 10
    y1, y2 = x1 * m + b, x2 * m + b
    line_mag = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
    if np.any(line_mag < 1e-8):
        warnings.warn("At least one line segment given by x1, y1, x2, y2 is very short.", stacklevel=2)
    u = ((px - x1) * (x2 - x1) + (py - y1) * (y2 - y1)) / (line_mag * line_mag)
    if np.any(u < 0.00001) or np.any(u > 1):
        # In R this branch compares whole vectors in an `if` (error in R >= 4.2). Use the
        # evidently intended per-point distance to the nearest segment endpoint.
        warnings.warn("Projection falls outside the segment; using endpoint distances.", stacklevel=2)
        dist = np.minimum(np.hypot(px - x1, py - y1), np.hypot(px - x2, py - y2))
    else:
        ix = x1 + u * (x2 - x1)
        iy = y1 + u * (y2 - y1)
        dist = np.hypot(px - ix, py - iy)

    full = np.concatenate([[np.nan], dist, [np.nan]])
    if return_index:
        return int(np.nanargmax(full)) + 1
    return full


def _add_correlated_features(ggc, seed_idx):
    """Greedily grow the feature set from ``seed_idx`` by repeatedly adding the gene with
    the highest correlation to any selected gene (ties: earlier selected gene, then
    earlier GGC row), exactly as the R neighbour-adding loop does."""
    n = ggc.shape[0]
    # Column j of the GGC holds gene j's correlations (R iterates over columns).
    cols = np.where(np.isnan(ggc), -np.inf, ggc)
    available = np.ones(n, dtype=bool)
    available[seed_idx] = False

    selected = list(seed_idx)
    best_val, best_idx = [], []

    def best_for(j):
        col = np.where(available, cols[:, j], -np.inf)
        k = int(np.argmax(col))
        return col[k], k

    for j in selected:
        v, k = best_for(j)
        best_val.append(v)
        best_idx.append(k)

    for _ in range(n - len(seed_idx)):
        bv = np.asarray(best_val)
        pick = int(np.argmax(bv))
        if bv[pick] == -np.inf:
            break  # no remaining candidate (only NaN correlations left)
        new = best_idx[pick]
        available[new] = False
        selected.append(new)
        for i, j in enumerate(selected[:-1]):
            if best_idx[i] == new:
                best_val[i], best_idx[i] = best_for(j)
        v, k = best_for(new)
        best_val.append(v)
        best_idx.append(k)
    return selected


def run_stepwise_reg(ggc, ggc_genes, num_steps=100, num_regressions=30):
    """Order genes by stepwise regression on the GGC.

    Returns
    -------
    feature_genes : numpy.ndarray
        All GGC genes in DUBStepR order (stepwise-selected seed genes up to the elbow,
        followed by correlated genes added greedily).
    elbow_pt : int
        Number of seed genes (elbow of the scree curve).
    scree : numpy.ndarray
        Scree values (index 0 is the Frobenius norm of the centred GGC).
    """
    logger.info("Running Stepwise Regression...")
    ggc = np.asarray(ggc, dtype=np.float64)
    ggc_genes = np.asarray(ggc_genes)
    g_cent = ggc - ggc.mean(axis=0, keepdims=True)

    scree = np.empty(num_steps + 1)
    scree[0] = np.linalg.norm(g_cent)
    feature_idx = []
    for i in range(1, num_steps + 1):
        if i < num_regressions:
            gg = g_cent.T @ g_cent
            gc_norm = np.sqrt((gg * gg).sum(axis=1))
            g_norm = np.sqrt(np.abs(np.diag(gg)))
            with np.errstate(divide="ignore", invalid="ignore"):
                var_exp = gc_norm / g_norm
            j = int(stable_order_desc(var_exp)[0])
            if j not in feature_idx:
                feature_idx.append(j)
            g = g_cent[:, j]
            gtg = g @ g
            g_cent = g_cent - np.outer(g, gg[j, :]) / gtg
            scree[i] = var_exp[j]
        else:
            # R: scree_values[num_regressions] is the value at position 30, i.e. step 29
            scree[i] = scree[num_regressions - 1]
    logger.info("Done.")

    scree_nz = scree[scree != 0]
    y = np.full(num_steps, np.nan)
    y[: min(num_steps, len(scree_nz))] = np.log(scree_nz[:num_steps])
    if np.isnan(y).any():
        raise RuntimeError("Too few non-zero scree values to locate the elbow.")
    elbow_id = find_elbow(y)
    elbow_id = min(elbow_id, len(feature_idx))

    logger.info("Adding correlated features...")
    ordered = _add_correlated_features(ggc, feature_idx[:elbow_id])
    logger.info("Done.")
    return ggc_genes[ordered], elbow_id, scree
