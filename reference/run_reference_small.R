# Reference outputs of R DUBStepR on its bundled pbmc_norm_small_data.
# Usage: Rscript run_reference_small.R [out_dir]   (default: tests/data)
library(DUBStepR); library(Matrix)
options(Seurat.object.assay.version = "v3")  # DUBStepR 1.2.0 needs Seurat v3-style assays
args <- commandArgs(trailingOnly = TRUE)
out <- file.path(if (length(args)) args[1] else "tests/data", "small_")

d <- DUBStepR::pbmc_norm_small_data
res <- DUBStepR(d)
write.table(res$corr.info, paste0(out, "corr_info.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(res$optimal.feature.genes, paste0(out, "optimal_genes.txt"))
write.table(data.frame(n = names(res$density.index), di = sprintf("%.17g", res$density.index)),
            paste0(out, "density_index.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

filt <- DUBStepR:::getFilteredData(d, 0.05 * ncol(d), "human")
g <- DUBStepR:::getGGC(filt)
write.table(data.frame(g = names(g$corr.range), z = sprintf("%.17g", g$corr.range)),
            paste0(out, "zrange.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("elbow", res$elbow.pt, "n_corr", nrow(res$corr.info), "n_opt", length(res$optimal.feature.genes), "\n")
