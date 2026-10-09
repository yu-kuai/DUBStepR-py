# Reference outputs of R DUBStepR on the vignette dataset (10x PBMC 1k, seuratObj.rda from the
# DUBStepR repository's vignettes/ folder), following the vignette.
# Usage: Rscript run_reference_vignette.R path/to/seuratObj.rda [out_dir]
library(DUBStepR); library(Matrix); library(Seurat)
options(Seurat.object.assay.version = "v3")
args <- commandArgs(trailingOnly = TRUE)
out <- file.path(if (length(args) > 1) args[2] else "tests/data", "vignette_")

load(args[1])
seuratObj <- UpdateSeuratObject(seuratObj)
seuratObj <- NormalizeData(object = seuratObj, normalization.method = "LogNormalize")
d <- GetAssayData(seuratObj, assay = "RNA", layer = "data")
writeMM(as(d, "CsparseMatrix"), paste0(out, "data.mtx"))   # input for the Python side
writeLines(rownames(d), paste0(out, "genes.txt"))

res <- DUBStepR(input.data = d, min.cells = 0.05 * ncol(d), optimise.features = TRUE, k = 10, num.pcs = 20, error = 0)
write.table(res$corr.info, paste0(out, "corr_info.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(res$optimal.feature.genes, paste0(out, "optimal_genes.txt"))
write.table(data.frame(n = names(res$density.index), di = sprintf("%.17g", res$density.index)),
            paste0(out, "density_index.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("elbow", res$elbow.pt, "n_corr", nrow(res$corr.info), "n_opt", length(res$optimal.feature.genes), "\n")
