# For the first density-index step (elbow-sized gene set), compare R DUBStepR's irlba-based PCA with an
# exact SVD. Usage: Rscript irlba_first_step_check.R <stress workdir>  (after run_stress.sh)
OUT <- paste0(sub("/$", "", commandArgs(trailingOnly = TRUE)[1]), "/")
suppressMessages({library(DUBStepR); library(Matrix); library(Seurat)})
options(Seurat.object.assay.version = "v3")
check <- function(name, species) {
  m <- as(readMM(paste0(OUT, name, ".mtx")), "CsparseMatrix"); rownames(m) <- readLines(paste0(OUT, name, ".genes")); colnames(m) <- paste0("c", seq_len(ncol(m)))
  ci <- read.delim(paste0(OUT, name, ".R.corrinfo")); elbow <- as.integer(readLines(paste0(OUT, name, ".R.elbow")))
  genes <- ci$g[1:elbow]
  so <- suppressWarnings(CreateSeuratObject(counts = m[genes, ]))
  so <- ScaleData(so, features = genes, verbose = FALSE)
  A <- t(GetAssayData(so, layer = "scale.data")); npcs <- min(50, ncol(A) - 1)
  di <- function(emb, d) { emb <- emb[, 1:min(20, ncol(emb))]; nn <- RANN::nn2(emb, k = 11)$nn.dists[, -1]
    mean(nn) / sqrt(sum((d / sqrt(nrow(A) - 1))[1:min(20, length(d))]^2)) }
  set.seed(42); ir <- suppressWarnings(irlba::irlba(A, nv = npcs))
  ex <- svd(A, nu = npcs, nv = npcs)
  cat(sprintf("%-14s step=%d  irlba DI=%.12f  exact-SVD DI=%.12f  max |d_irlba-d_exact|/d = %.2e\n", name, elbow,
      di(ir$u %*% diag(ir$d), ir$d), di(ex$u %*% diag(ex$d[1:npcs]), ex$d[1:npcs]), max(abs(ir$d - ex$d[1:npcs]) / ex$d[1:npcs])))
}
check("synth400x2000", "human"); check("rat", "rat"); check("genes60", "human"); check("synth1200x3000", "human")
