#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(pathview))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: run_pathview.R <ko_table_dir> <kegg_dir> <pathway_id...>")
}

ko_dir <- args[1]
kegg_dir <- args[2]
pathways <- args[3:length(args)]

for (pw in pathways) {
  tsv <- file.path(ko_dir, paste0(pw, "_ko.tsv"))
  if (!file.exists(tsv)) {
    warning(paste("Missing KO table:", tsv))
    next
  }
  df <- read.delim(tsv, stringsAsFactors = FALSE)
  rownames(df) <- df$ko
  frac_cols <- grep("_frac$", names(df), value = TRUE)
  if (length(frac_cols) < 2) {
    warning(paste("Missing fraction columns in", tsv))
    next
  }
  mat <- as.matrix(df[, frac_cols])
  colnames(mat) <- sub("_frac$", "", colnames(mat))

  # Difference map: group2 - group1
  if (ncol(mat) >= 2) {
    diff <- mat[, 2] - mat[, 1]
    mat <- as.matrix(diff)
    colnames(mat) <- paste0(colnames(df)[which(names(df) == frac_cols[2])])
  }

  pathview(
    gene.data = mat,
    pathway.id = gsub("^ko", "", pw),
    species = "ko",
    kegg.native = TRUE,
    kegg.dir = kegg_dir,
    out.suffix = "ISS_minus_Earth",
    key.pos = "bottomright"
  )
}
