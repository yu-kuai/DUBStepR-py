"""Compare dubsteppy against outputs of the R package DUBStepR 1.2.0.

Reference files in tests/data were produced by reference/run_reference_*.R
(R 4.4.3, Seurat 5.5.1 with v3 assays, qlcMatrix 0.9.9, RANN 2.6.2, irlba 2.4.1).
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import dubsteppy

DATA = Path(__file__).parent / "data"


def _read_ref(prefix):
    ci = pd.read_csv(DATA / f"{prefix}_corr_info.tsv", sep="\t", keep_default_na=False)
    opt = (DATA / f"{prefix}_optimal_genes.txt").read_text().split()
    di = pd.read_csv(DATA / f"{prefix}_density_index.tsv", sep="\t")
    return ci, opt, di


def _check(res, prefix, elbow):
    ci, opt, di = _read_ref(prefix)
    assert res.elbow_pt == elbow
    assert list(res.corr_info["feature.genes"]) == list(ci["feature.genes"])
    np.testing.assert_allclose(res.corr_info["corr.range"], ci["corr.range"], rtol=0, atol=1e-10)
    assert list(res.optimal_feature_genes) == opt
    assert list(res.density_index.index) == list(di["n"])
    # R uses irlba (approximate truncated SVD); dubsteppy solves the PCA exactly
    np.testing.assert_allclose(res.density_index.to_numpy(), di["di"].to_numpy(), rtol=1e-6)


@pytest.fixture(scope="module")
def small_result():
    return dubsteppy.dubstepr(dubsteppy.pbmc_norm_small_data())


def test_small_dataset_matches_r(small_result):
    _check(small_result, "small", elbow=17)


def test_small_dataset_corr_range_matches_r(small_result):
    z = pd.read_csv(DATA / "small_zrange.tsv", sep="\t", keep_default_na=False)
    assert list(small_result.corr_range.index) == list(z["g"])
    np.testing.assert_allclose(small_result.corr_range.to_numpy(), z["z"].to_numpy(), atol=1e-10)


def test_vignette_dataset_matches_r():
    mat = sp.load_npz(DATA / "vignette_pbmc1k_lognorm.npz")
    genes = (DATA / "vignette_pbmc1k_genes.txt").read_text().split("\n")
    res = dubsteppy.dubstepr(mat, gene_names=genes)
    _check(res, "vignette", elbow=15)


def test_input_types_agree(small_result):
    df = dubsteppy.pbmc_norm_small_data()
    res = dubsteppy.dubstepr(sp.csc_matrix(df.to_numpy()), gene_names=df.index, optimise_features=False)
    assert list(res.feature_genes) == list(small_result.feature_genes)
    assert res.optimal_feature_genes is None


def test_anndata_wrapper(small_result):
    anndata = pytest.importorskip("anndata")
    df = dubsteppy.pbmc_norm_small_data()
    adata = anndata.AnnData(X=sp.csr_matrix(df.T.to_numpy()), obs=pd.DataFrame(index=df.columns),
                            var=pd.DataFrame(index=df.index))
    res = dubsteppy.dubstepr_anndata(adata)
    assert list(res.optimal_feature_genes) == list(small_result.optimal_feature_genes)
    assert adata.var["dubstepr_feature"].sum() == len(res.optimal_feature_genes)
    assert adata.var.loc[res.feature_genes[0], "dubstepr_rank"] == 1


def test_filtering_removes_mito_ribo_ercc():
    genes = ["MT-CO1", "RPL13", "rps3", "ERCC-00002", "CD3E", "ENSG00000198804"]
    mat = np.ones((len(genes), 20))
    filt, kept = dubsteppy.get_filtered_data(mat, gene_names=genes)
    # As in R, annotation lookup is exact-case: "rps3" is not a listed symbol, so it is kept
    assert list(kept) == ["rps3", "CD3E"]


def test_find_elbow():
    y = np.log(np.r_[100.0, 50, 20, 10, 9, 8.5, 8.2, 8.0, 7.9, 7.8])
    assert dubsteppy.find_elbow(y) == 4


def test_log_normalize():
    counts = np.array([[1.0, 0.0], [3.0, 5.0]])
    out = dubsteppy.log_normalize(counts)
    np.testing.assert_allclose(out, np.log1p(counts / (counts.sum(0) / 1e4)))
    np.testing.assert_allclose(dubsteppy.log_normalize(sp.csr_matrix(counts)).toarray(), out)
