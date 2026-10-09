suppressMessages({library(DUBStepR); library(Matrix); library(jsonlite)})
options(Seurat.object.assay.version = "v3")
OUT <- paste0(sub("/$", "", commandArgs(trailingOnly = TRUE)[1]), "/")
jobs <- fromJSON(paste0(OUT, "jobs.json"), simplifyDataFrame = FALSE)
for (j in jobs) {
  m <- as(readMM(paste0(OUT, j$name, ".mtx")), "CsparseMatrix")
  rownames(m) <- readLines(paste0(OUT, j$name, ".genes")); colnames(m) <- paste0("c", seq_len(ncol(m)))
  mc <- if (identical(j$min_cells, "NA")) 0.05 * ncol(m) else j$min_cells
  t0 <- Sys.time()
  res <- tryCatch(suppressMessages(suppressWarnings(DUBStepR(m, min.cells = mc, species = j$species,
            optimise.features = j$optimise, k = j$k, num.pcs = j$num_pcs))), error = function(e) e)
  if (inherits(res, "error")) { writeLines(conditionMessage(res), paste0(OUT, j$name, ".R.error")); cat(j$name, "ERROR\n"); next }
  write.table(data.frame(g = res$corr.info$feature.genes, z = sprintf("%.17g", res$corr.info$corr.range)),
              paste0(OUT, j$name, ".R.corrinfo"), sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(as.character(res$elbow.pt), paste0(OUT, j$name, ".R.elbow"))
  if (j$optimise) {
    writeLines(res$optimal.feature.genes, paste0(OUT, j$name, ".R.opt"))
    write.table(data.frame(n = names(res$density.index), di = sprintf("%.17g", res$density.index)),
                paste0(OUT, j$name, ".R.di"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  cat(sprintf("%-18s R done in %.1fs\n", j$name, as.numeric(Sys.time() - t0, units = "secs")))
}
