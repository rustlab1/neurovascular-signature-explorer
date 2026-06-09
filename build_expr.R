#!/usr/bin/env Rscript
# build_expr.R -- average expression per gene x cell-type x condition (control /
# disease) for ALL expressed genes, from the per-disease Seurat objects. Powers
# the "absolute expression" view of gene search (confirm markers, see levels for
# non-DEG genes). Writes expr.js (window.NV_EXPR).
#
# Values are mean log-normalized expression, int-encoded x100 to keep the file small.
suppressMessages({ library(Seurat); library(Matrix) })

DATA <- Sys.getenv("NV_DATA_DIR", "source_data")
OUT  <- "expr.js"
FLOOR <- 0.1   # keep genes whose max group-mean expression >= this (drops the unexpressed tail)

LIN_ORDER <- c("Endothelial","Pericyte","vSMC","M_Fibro","P_Fibro","P_Mac",
               "Microglia","Astro","Oligo","OPC","Neuron_exc","Neuron_inh")

files <- list(
  AD  = file.path(DATA, "Analysis/1_Obj_Clusterin_QC/Objects/AD_cortex_small_snRNAseq_harmony_reclustered.rds"),
  HD  = file.path(DATA, "Analysis/1_Obj_Clusterin_QC/Objects/HD_small_pericytes_snRNAseq_harmony_reclustered.rds"),
  FTD = file.path(DATA, "Analysis/1_Obj_Clusterin_QC/Objects/DM_small_snRNAseq_harmony_reclustered.rds")
)

jsarr <- function(v) paste(v, collapse = ",")
qstr  <- function(v) paste0('"', v, '"', collapse = ",")

blocks <- character(0)
all_genes <- character(0)

for (dcode in names(files)) {
  cat("Reading", dcode, "...\n")
  o <- readRDS(files[[dcode]])
  DefaultAssay(o) <- "RNA"
  o <- tryCatch(JoinLayers(o), error = function(e) o)
  m <- tryCatch(GetAssayData(o, assay = "RNA", slot = "data"), error = function(e) NULL)
  if (is.null(m) || length(m@x) == 0 || max(m@x) == 0) {
    o <- NormalizeData(o, verbose = FALSE)
    m <- GetAssayData(o, assay = "RNA", slot = "data")
  }
  lin  <- as.character(o$lineage)
  cond <- ifelse(grepl("ctrl|control", as.character(o$condition), ignore.case = TRUE), "C", "D")
  lins <- LIN_ORDER[LIN_ORDER %in% unique(lin)]
  genes <- rownames(m); ng <- length(genes); nl <- length(lins)
  avgC <- matrix(0, ng, nl); avgD <- matrix(0, ng, nl)
  for (j in seq_along(lins)) {
    cC <- which(lin == lins[j] & cond == "C")
    cD <- which(lin == lins[j] & cond == "D")
    if (length(cC) > 0) avgC[, j] <- Matrix::rowMeans(m[, cC, drop = FALSE])
    if (length(cD) > 0) avgD[, j] <- Matrix::rowMeans(m[, cD, drop = FALSE])
  }
  keep <- which(pmax(apply(avgC, 1, max), apply(avgD, 1, max)) >= FLOOR)
  genes <- genes[keep]; avgC <- round(avgC[keep, , drop = FALSE] * 100); avgD <- round(avgD[keep, , drop = FALSE] * 100)
  # flatten gene-major: index = (g-1)*nl + l
  flatC <- as.integer(t(avgC)); flatD <- as.integer(t(avgD))
  blocks <- c(blocks, sprintf('"%s":{"lin":[%s],"genes":[%s],"avgC":[%s],"avgD":[%s]}',
                              dcode, qstr(lins), qstr(genes), jsarr(flatC), jsarr(flatD)))
  all_genes <- union(all_genes, genes)
  cat("  ", dcode, ":", length(genes), "genes,", nl, "lineages\n")
  rm(o, m, avgC, avgD); gc(verbose = FALSE)
}

all_genes <- sort(all_genes)
con <- file(OUT, "w")
cat("window.NV_EXPR={", paste(blocks, collapse = ","),
    ',"allGenes":[', qstr(all_genes), "]};\n", sep = "", file = con)
close(con)
cat("Wrote", OUT, "(", length(all_genes), "genes total )\n")
