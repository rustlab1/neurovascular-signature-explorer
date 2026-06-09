#!/usr/bin/env Rscript
# build_umap.R -- extract a downsampled harmonized UMAP embedding (coordinates +
# cell-type + condition) from the per-disease Seurat objects, for the atlas view.
# Writes umap.js (window.NV_UMAP). Each disease has its own harmonized embedding,
# so the atlas shows one disease at a time (selectable in the page).
suppressMessages(library(Seurat))

# Point NV_DATA_DIR at the source analysis directory (scRNAseq_human/Data).
DATA <- Sys.getenv("NV_DATA_DIR", "source_data")
OUT  <- "umap.js"
N    <- 5000  # cells sampled per disease

files <- list(
  AD  = file.path(DATA, "Analysis/1_Obj_Clusterin_QC/Objects/AD_cortex_small_snRNAseq_harmony_reclustered.rds"),
  HD  = file.path(DATA, "Analysis/1_Obj_Clusterin_QC/Objects/HD_small_pericytes_snRNAseq_harmony_reclustered.rds"),
  FTD = file.path(DATA, "Analysis/1_Obj_Clusterin_QC/Objects/DM_small_snRNAseq_harmony_reclustered.rds")
)

set.seed(42)
all_rows <- character(0)

for (dcode in names(files)) {
  f <- files[[dcode]]
  cat("Reading", dcode, "...\n")
  o <- readRDS(f)
  red <- if ("umap_harmony" %in% Reductions(o)) "umap_harmony" else "umap"
  emb <- Embeddings(o, reduction = red)
  md  <- o@meta.data
  lin <- as.character(md$lineage)
  cond <- ifelse(grepl("ctrl|control", as.character(md$condition), ignore.case = TRUE), 0L, 1L)

  keep <- which(!is.na(lin) & !is.na(emb[, 1]) & !is.na(emb[, 2]))
  if (length(keep) > N) keep <- sample(keep, N)

  x <- round(emb[keep, 1], 2)
  y <- round(emb[keep, 2], 2)
  rows <- sprintf('[%s,%s,"%s","%s",%d]', x, y, lin[keep], dcode, cond[keep])
  all_rows <- c(all_rows, rows)
  cat("  ", dcode, ":", length(keep), "cells (", red, ")\n")
  rm(o, emb, md); gc(verbose = FALSE)
}

con <- file(OUT, "w")
cat("window.NV_UMAP=[", paste(all_rows, collapse = ","), "];\n", sep = "", file = con)
close(con)
cat("Wrote", OUT, "(", length(all_rows), "points )\n")
