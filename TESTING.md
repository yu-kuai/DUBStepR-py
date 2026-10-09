# Testing and validation

dubsteppy is a pure-Python rewrite of the R package
[DUBStepR](https://github.com/prabhakarlab/DUBStepR) 1.2.0. It does not call R. To confirm that
the rewrite is faithful, its output was compared with the original R package on 2 real reference
datasets (pytest suite) and 16 additional datasets and parameter settings (stress test).

**Result:** across all 18 comparisons, dubsteppy reproduces R DUBStepR's ranked feature genes,
elbow point and optimal feature set exactly. Correlation-range scores agree to within 2e-12. The
density index agrees to within 1e-8 except at a single step in 5 runs, where R's approximate SVD
(irlba) is inaccurate and dubsteppy matches the exact SVD (see
[Density-index differences](#density-index-differences)).

## Environment

All comparisons were run inside the Apptainer image built from `container/dubsteppy.def`
(`container/dubsteppy.sif`, not tracked in git):

| | version |
|---|---|
| R | 4.4.3 |
| DUBStepR | 1.2.0 (GitHub `prabhakarlab/DUBStepR`) |
| Seurat / SeuratObject | 5.5.1 / 5.4.0 |
| qlcMatrix / RANN / irlba | 0.9.9 / 2.6.2 / 2.4.1 |
| Python / numpy / scipy / pandas / anndata | 3.12.15 / 2.5.3 / 1.18.1 / 3.0.6 / 0.13.4 |

DUBStepR 1.2.0 fails under Seurat v5 (`No layer matching pattern 'data' found`), so R was run with
`options(Seurat.object.assay.version = "v3")`. This restores the Seurat v3/v4 behaviour the
package was written for.

## What is compared

Every comparison checks each output of `DUBStepR()`:

| output | R | dubsteppy | criterion |
|---|---|---|---|
| ranked feature genes | `corr.info$feature.genes` | `res.corr_info["feature.genes"]` | identical, same order |
| correlation range | `corr.info$corr.range` | `res.corr_info["corr.range"]` | abs. diff |
| elbow point | `elbow.pt` | `res.elbow_pt` | identical |
| optimal feature set | `optimal.feature.genes` | `res.optimal_feature_genes` | identical, same order |
| density index | `density.index` | `res.density_index` | same steps; rel. diff |

## 1. pytest suite (`tests/test_against_r.py`)

Run it with `pytest tests` (8 tests, about 8 s). It needs only the files in the repo.

| test | what it checks | result |
|---|---|---|
| `test_small_dataset_matches_r` | R's bundled `pbmc_norm_small_data` (230 genes x 80 cells), default parameters: all outputs above | identical; DI rel. diff 3.4e-9 |
| `test_small_dataset_corr_range_matches_r` | z-scored correlation range of all 216 filtered genes, in R's order (an intermediate of `getGGC`) | identical order; max diff 8e-15 |
| `test_vignette_dataset_matches_r` | Vignette dataset: 10x PBMC 1k (33,538 genes x 996 cells), log-normalised with Seurat `NormalizeData`, as in the vignette | 870/870 genes, elbow 15/15, optimal 65/65 identical; DI rel. diff 4.3e-9 over 35 steps |
| `test_input_types_agree` | DataFrame and scipy sparse input give the same ranking; `optimise_features=False` works | pass |
| `test_anndata_wrapper` | `dubstepr_anndata` gives the same genes and writes `adata.var` correctly | pass |
| `test_filtering_removes_mito_ribo_ercc` | removal of mitochondrial, ribosomal, ERCC spike-in and Ensembl-ID-annotated genes | pass |
| `test_find_elbow` | elbow detection on a known curve | pass |
| `test_log_normalize` | `log1p(count / library_size * 1e4)`, dense and sparse | pass |

The R outputs in `tests/data/` were produced by `reference/run_reference_small.R` and
`reference/run_reference_vignette.R` (originally in a conda R 4.4.3 env with the same package
versions). Rerunning them in the container reproduces the committed files byte for byte, except
for one `corr.range` value in the vignette output (NAMPT), which differs in the 15th significant
digit from BLAS rounding. The vignette test uses `tests/data/vignette_pbmc1k_lognorm.npz`,
the exact matrix R received. Genes with zero counts were removed to keep the file small; they
can never pass DUBStepR's expression filter, so this does not change the result.

## 2. Stress test (`tests/stress/`)

Run it with `tests/stress/run_stress.sh [workdir]`. It takes about 12 min, almost all of it in
R. The script:

1. `make_jobs.py` generates 16 inputs from fixed seeds.
2. `run_r.R` runs R DUBStepR on each.
3. `compare.py` runs dubsteppy and compares the outputs.
4. `irlba_first_step_check.R` investigates the density-index differences.

Generated files go to `workdir` and are not committed. The table below is saved as
`tests/stress/results_2026-10-09.txt`.

**The 16 inputs:**

| job | data | what it exercises |
|---|---|---|
| `sub150`, `sub300`, `sub500`, `sub750` | random subsets of 150 to 750 cells from the vignette data | different cell numbers, so different filtering, GGC and kNN structure |
| `genes30`, `genes60` | random 30% / 60% of vignette genes | different GGC size and expression-bin layout |
| `k5_pcs10`, `k20_pcs30` | full vignette, `k`=5/`num.pcs`=10 and `k`=20/`num.pcs`=30 | density-index parameters |
| `mincells20`, `mincells150` | full vignette, `min.cells`=20 / 150 | expression filter threshold (larger and smaller GGC) |
| `mouse`, `rat` | full vignette, `species`="mouse" / "rat" | species gene-annotation filtering |
| `noopt` | full vignette, `optimise.features=FALSE` | ranking-only mode |
| `synth400x2000`, `synth1200x3000`, `synth2500x4000` | simulated negative-binomial counts (4 / 8 / 12 clusters with 40 marker genes each, cell-size variation), log-normalised | data unlike PBMC; larger elbow (up to 29) and optimal sets (up to 379 genes) |

**Results** (`ggc`, `elbow` and `opt` are shown as dubsteppy/R; `zdiff` is the max abs. diff of
corr.range; `di_reldiff` is the max rel. diff of the density index):

| job | ggc genes | order identical | zdiff | elbow | optimal set | optimal identical | DI steps | di_reldiff |
|---|---|---|---|---|---|---|---|---|
| sub150 | 745/745 | ✅ | 8.8e-13 | 14/14 | 64/64 | ✅ | 30 | 1.9e-10 |
| sub300 | 711/711 | ✅ | 6.2e-13 | 15/15 | 40/40 | ✅ | 28 | 5.0e-09 |
| sub500 | 765/765 | ✅ | 6.6e-13 | 14/14 | 39/39 | ✅ | 31 | 8.5e-11 |
| sub750 | 819/819 | ✅ | 3.1e-13 | 14/14 | 39/39 | ✅ | 33 | 7.3e-10 |
| genes30 | 246/246 | ✅ | 4.2e-13 | 13/13 | 13/13 | ✅ | 10 | 1.1e-08 |
| genes60 | 536/536 | ✅ | 9.7e-13 | 14/14 | 14/14 | ✅ | 21 | 1.2e-03 † |
| k5_pcs10 | 870/870 | ✅ | 7.3e-13 | 15/15 | 115/115 | ✅ | 35 | 1.6e-15 |
| k20_pcs30 | 870/870 | ✅ | 7.3e-13 | 15/15 | 15/15 | ✅ | 35 | 4.4e-09 |
| mincells20 | 1131/1131 | ✅ | 7.7e-13 | 17/17 | 42/42 | ✅ | 45 | 2.4e-09 |
| mincells150 | 354/354 | ✅ | 7.7e-13 | 16/16 | 91/91 | ✅ | 14 | 4.3e-10 |
| mouse | 889/889 | ✅ | 1.6e-12 | 18/18 | 93/93 | ✅ | 35 | 4.6e-09 |
| rat | 889/889 | ✅ | 1.6e-12 | 18/18 | 93/93 | ✅ | 35 | 7.7e-03 † |
| noopt | 870/870 | ✅ | 7.3e-13 | 15/15 | — | — | — | — |
| synth400x2000 | 164/164 | ✅ | 7.3e-14 | 11/11 | 136/136 | ✅ | 7 | 3.6e-03 † |
| synth1200x3000 | 288/288 | ✅ | 9.4e-14 | 18/18 | 243/243 | ✅ | 11 | 4.0e-03 † |
| synth2500x4000 | 429/429 | ✅ | 7.2e-13 | 29/29 | 379/379 | ✅ | 17 | 3.1e-15 |

† See below. In each of these runs only the first density-index step differs; every other step
agrees to within about 1e-8.

Runtime: R took 6 to 74 s per job; dubsteppy took 0.2 to 8 s.

### Density-index differences

In 5 runs the density index differs from R by 0.1 to 0.8%, always only at the **first step**
(the elbow-sized gene set, e.g. 11 to 18 genes). At that step Seurat's `RunPCA` asks irlba for
`n_genes − 1` singular vectors of an `n_cells × n_genes` matrix. irlba is a truncated method and
is unreliable when asked for almost all singular values (it warns *"You're computing too large a
percentage of total singular values"*). dubsteppy computes the PCA exactly.

`tests/stress/irlba_first_step_check.R` recomputes that step in R both ways:

| job | step (genes) | R with irlba (= DUBStepR) | R with exact `svd()` | dubsteppy | worst irlba singular-value error |
|---|---|---|---|---|---|
| synth400x2000 | 11 | 0.649847 | 0.652203 | 0.652203 | 0.6% |
| rat | 18 | 0.401993 | 0.405094 | 0.405094 | 26.6% |
| genes60 | 14 | 0.365625 | 0.366080 | 0.366080 | 7.3% |
| synth1200x3000 | 18 | 0.722952 | 0.725863 | 0.725863 | 4.2% |

dubsteppy matches R's exact SVD, so these differences are approximation error in the original
implementation, not a porting error. In none of the 16 runs did this change the optimal feature
set. In principle it could, if the first step's density index is within about 1% of the minimum.
The pytest suite uses `rtol=1e-6` for the density index, and both reference datasets are within it.

## Not covered

* **More than 10,000 cells.** DUBStepR then switches to approximate kNN (`error = 1`). dubsteppy
  uses `scipy.spatial.cKDTree(eps=1)`. It has the same (1+eps) guarantee as RANN/ANN, but
  approximate neighbours can differ, so results are expected to be close but not identical. Not
  tested against R.
* **Gene names containing `_`.** Seurat renames them and R DUBStepR then silently drops them from
  the PCA; dubsteppy keeps them, so results differ by design.
* `align_bins=True` (the corrected expression binning) deliberately differs from R and is
  therefore not compared.
