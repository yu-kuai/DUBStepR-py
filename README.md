# dubsteppy

Python port of the R package [DUBStepR](https://github.com/prabhakarlab/DUBStepR) 1.2.0
(Determining the Underlying Basis using Stepwise Regression), a correlation-based feature
selection method for clustering single-cell RNA-seq data.

> Ranjan B. et al. *DUBStepR is a scalable correlation-based feature selection method for
> accurately clustering single-cell data.* Nat Commun 12, 5849 (2021).

The port reproduces the R package's results: the same feature genes in the same order, the
same elbow point and the same optimal feature set on both the dataset bundled with DUBStepR
and the PBMC 1k dataset from its vignette (see [Validation](#validation)).

## Installation

**1. Clone the repository**

```bash
git clone https://github.com/yu-kuai/DUBStepPy.git
cd DUBStepPy
```

**2. Create and activate a conda/mamba environment**

```bash
mamba env create -f environment.yml    # or: conda env create -f environment.yml
mamba activate dubsteppy               # or: conda activate dubsteppy
```

This installs Python with numpy, scipy, pandas, anndata and pytest from conda-forge. To use an
existing environment instead, skip this step; pip installs the dependencies in step 3.

**3. Install dubsteppy**

```bash
pip install .                # add -e for an editable (development) install
```

**4. (Optional) Check the installation**

```bash
pytest tests                 # 8 tests, about 10 s; compares against R DUBStepR outputs
python -c "import dubsteppy; print(dubsteppy.dubstepr(dubsteppy.pbmc_norm_small_data()).elbow_pt)"  # 17
```

Without conda: `pip install ".[anndata]"` into any Python >= 3.9 environment.

## Usage

Input is **log-normalised** expression, as for the R package.

```python
import dubsteppy

# genes x cells (DataFrame, or ndarray / scipy sparse + gene_names)
data = dubsteppy.pbmc_norm_small_data()
res = dubsteppy.dubstepr(data)            # same defaults as R DUBStepR()

res.optimal_feature_genes                 # optimal feature set
res.corr_info                             # feature.genes (DUBStepR order) + corr.range
res.elbow_pt
res.density_index                         # indexed by number of genes
```

With AnnData (cells x genes, e.g. after `sc.pp.normalize_total(target_sum=1e4); sc.pp.log1p()`):

```python
res = dubsteppy.dubstepr_anndata(adata, layer=None)
adata.var["dubstepr_feature"]   # bool, optimal feature set -> use as highly_variable
adata.var["dubstepr_rank"]      # DUBStepR order of correlated genes
```

### Parameters (R name → Python name)

| R | Python | default |
|---|---|---|
| `input.data` | `data` (+ `gene_names`) | genes x cells |
| `min.cells` | `min_cells` | `0.05 * n_cells` |
| `species` | `species` | `"human"` (`"mouse"`, `"rat"`) |
| `optimise.features` | `optimise_features` | `True` |
| `k` | `k` | 10 |
| `num.pcs` | `num_pcs` | 20 |
| `error` | `error` | 0 (forced to 1 for > 10,000 cells, as in R) |
| — | `align_bins` | `False` (see below) |

The individual steps are exposed too: `get_filtered_data`, `get_ggc`, `get_correlation_range`,
`run_stepwise_reg`, `find_elbow`, `get_optimal_feature_set`, `log_normalize`,
`load_gene_annotation`.

## Validation

`tests/test_against_r.py` compares against R DUBStepR 1.2.0 (R 4.4.3, Seurat 5.5.1,
qlcMatrix 0.9.9, RANN 2.6.2, irlba 2.4.1). The R scripts used are in `reference/`.

| dataset | genes in GGC | elbow | optimal set | gene order | corr.range | density index |
|---|---|---|---|---|---|---|
| `pbmc_norm_small_data` (230 × 80) | 62 / 62 | 17 / 17 | 42 / 42 | identical | < 1e-14 | rel. diff < 4e-9 |
| vignette PBMC 1k (33,538 × 996) | 870 / 870 | 15 / 15 | 65 / 65 | identical | < 1e-12 | rel. diff < 5e-9 |

The small density-index differences come from R's approximate truncated SVD (irlba), which
dubsteppy replaces with an exact eigendecomposition. Runtime on the vignette dataset: about 4 s
in Python and about 38 s in R.

```bash
pytest tests
```

A broader stress test (16 subsampled, synthetic and parameter-varied datasets, all matching R)
lives in `tests/stress/`. See [TESTING.md](TESTING.md) for the full report, including why the
density index at the first step can differ from R by up to about 1%: R's irlba is inaccurate
there, and dubsteppy matches R's own exact SVD.

## Notes on fidelity to the R implementation

* **Expression bins in `getGGC`.** In R, the correlation-range vector is sorted before
  `tapply` is called with the expression bins in the original gene order, so each gene is
  z-scored inside the bin of whichever gene originally held its rank position. dubsteppy
  reproduces this by default so the results match. Pass `align_bins=True` to bin each gene
  by its own mean expression, which is probably what the authors intended. Results will then
  differ from R.
* **Seurat v5.** R DUBStepR 1.2.0 fails under Seurat v5 (`No layer matching pattern 'data'`).
  The R reference outputs were generated with `options(Seurat.object.assay.version = "v3")`,
  which restores the Seurat v3/v4 behaviour the package was written for. dubsteppy has no
  Seurat dependency. It reimplements `ScaleData` (z-score, clipped above at 10) and
  `RunPCA` (uncentred, with zero-variance features dropped).
* **`findElbow`.** If any projection falls outside the reference segment, R ≥ 4.2 throws an
  error (`if` on a vector). dubsteppy uses the per-point endpoint distance and issues a
  warning instead.
* **Approximate kNN** (`error > 0`, used automatically above 10,000 cells) uses
  `scipy.spatial.cKDTree(eps=...)`. It has the same (1+eps) guarantee as RANN/ANN, but the
  approximate neighbours can differ, so above 10,000 cells the density index will not match R
  bit for bit.
* Gene names containing `_` are renamed by Seurat's `CreateSeuratObject` and then silently
  dropped from R's PCA. dubsteppy keeps them.

## Container

`container/dubsteppy.def` builds an Apptainer image with R DUBStepR (plus Seurat, qlcMatrix,
RANN, irlba and matrixcalc) and Python (numpy, scipy, pandas, anndata, pytest). It is used to
regenerate the references and run the tests:

```bash
apptainer build dubsteppy.sif container/dubsteppy.def
apptainer exec dubsteppy.sif bash -c 'PYTHONPATH=src python -m pytest tests'
apptainer exec dubsteppy.sif Rscript reference/run_reference_small.R
```

## License

dubsteppy is derived from R DUBStepR and is distributed under **the original DUBStepR user
license** (see `LICENSE`, (c) Shyam Prabhakar, Shyam's Lab, Genome Institute of Singapore).
That license is *not* an open-source license: academic, non-commercial use is permitted, but
redistribution is limited to the user's research team and immediate collaborators, and
commercial use needs a written licence from the author. (The R package's DESCRIPTION says
"MIT + file LICENSE", but the LICENSE file itself contains these restrictive terms.) The gene
annotation tables in `src/dubsteppy/data` and the `pbmc_norm_small_data` example are taken from
R DUBStepR.

If you use dubsteppy, cite: Ranjan B, Sun W, Park J, et al. DUBStepR is a scalable
correlation-based feature selection method for accurately clustering single-cell data.
*Nat Commun* 12, 5849 (2021). https://doi.org/10.1038/s41467-021-26085-2
