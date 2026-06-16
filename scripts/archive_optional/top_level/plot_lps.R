#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[1] else "."
metadata_path <- if (length(args) >= 2 && nzchar(args[2])) args[2] else NA
out_suffix <- if (length(args) >= 3) args[3] else ""
remove_oantigen_title <- length(args) >= 4 && tolower(args[4]) %in% c("1", "true", "yes")

out_dir <- file.path(base_dir, "plots")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

need_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Missing package: ", pkg, ". Install with install.packages(\"", pkg, "\")"))
  }
}

need_pkg("ggplot2")
need_pkg("ggrepel")
need_pkg("pheatmap")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

rename_isolates <- c(
  "Pantoea_piersonii_F8_6S_D2_EB4_spades" = "Pantoea_piersonii_F9_6S_D2_EB4_spades"
)

read_metadata <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    return(NULL)
  }
  meta <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(meta) <- trimws(names(meta))
  if (!"Sample ID" %in% names(meta)) {
    stop("Metadata file must include a 'Sample ID' column.")
  }
  rownames(meta) <- trimws(meta[["Sample ID"]])
  meta
}

metadata_df <- read_metadata(metadata_path)

apply_rename_vector <- function(x) {
  out <- x
  for (old in names(rename_isolates)) {
    out[out == old] <- rename_isolates[[old]]
  }
  out
}

apply_rename_columns <- function(mat) {
  colnames(mat) <- apply_rename_vector(colnames(mat))
  mat
}

flight_rank <- function(isolate) {
  if (grepl("^Pantoea_piersonii_IIIF", isolate)) {
    return(3L)
  }
  hit <- regmatches(isolate, regexpr("_F[0-9]+", isolate))
  if (length(hit) == 1 && nchar(hit) > 0) {
    return(as.integer(sub("_F", "", hit)))
  }
  999L
}

metadata_group <- function(isolate, metadata) {
  if (is.null(metadata) || !isolate %in% rownames(metadata)) {
    return("Unknown")
  }
  meta <- metadata[isolate, , drop = FALSE]
  category <- meta[["Category"]]
  location <- meta[["Location"]]
  is_piersonii <- grepl("piersonii", isolate)
  if (!is_piersonii) {
    return("Non-piersonii")
  }
  if (identical(location, "ISS")) {
    return("ISS")
  }
  if (identical(category, "Clinical")) {
    return("Earth clinical piersonii")
  }
  "Earth non-clinical piersonii"
}

metadata_order <- function(isolates, metadata) {
  if (is.null(metadata)) {
    return(isolates)
  }
  order_df <- data.frame(
    isolate = isolates,
    block = 4,
    category = "",
    source = "",
    flight = 999,
    stringsAsFactors = FALSE
  )
  matched <- isolates %in% rownames(metadata)
  order_df$category[matched] <- metadata[isolates[matched], "Category"]
  location <- rep("", length(isolates))
  location[matched] <- metadata[isolates[matched], "Location"]
  order_df$source[matched] <- metadata[isolates[matched], "Source"]
  is_piersonii <- grepl("piersonii", isolates)
  order_df$block[!is_piersonii] <- 0
  order_df$block[is_piersonii & location != "ISS" & order_df$category != "Clinical"] <- 1
  order_df$block[is_piersonii & location != "ISS" & order_df$category == "Clinical"] <- 2
  order_df$block[is_piersonii & location == "ISS"] <- 3
  order_df$sort_source <- ifelse(order_df$block == 3, order_df$source, "")
  order_df$flight <- vapply(isolates, flight_rank, integer(1))
  order_df$isolate[order(order_df$block, order_df$category, order_df$sort_source, order_df$isolate)]
}

output_name <- function(tag, plot_name) {
  paste0(tag, "_", plot_name, out_suffix, ".png")
}

plot_set <- function(tag, label) {
  summary_path <- file.path(base_dir, paste0(tag, "_summary.tsv"))
  jaccard_path <- file.path(base_dir, paste0(tag, "_jaccard.tsv"))
  prevalence_path <- file.path(base_dir, paste0(tag, "_prevalence_diff.tsv"))
  presence_path <- file.path(base_dir, paste0(tag, "_presence.tsv"))
  stats_path <- file.path(base_dir, paste0(tag, "_group_stats.tsv"))
  fisher_path <- file.path(base_dir, paste0(tag, "_gene_fisher.tsv"))

  if (!file.exists(summary_path) || !file.exists(jaccard_path) ||
      !file.exists(prevalence_path) || !file.exists(presence_path)) {
    message("Skipping ", tag, ": missing inputs")
    return()
  }

  summary_df <- read.delim(summary_path, stringsAsFactors = FALSE)
  summary_df$isolate <- apply_rename_vector(summary_df$isolate)
  summary_df$group <- factor(summary_df$group, levels = c("Earth", "ISS"))

  stats_df <- NULL
  if (file.exists(stats_path)) {
    stats_df <- read.delim(stats_path, stringsAsFactors = FALSE)
  }
  mean_diff <- NA
  mean_p <- NA
  perm_F <- NA
  perm_p <- NA
  if (!is.null(stats_df)) {
    md <- stats_df[stats_df$metric == "mean_count_diff", ]
    if (nrow(md) == 1) {
      mean_diff <- md$value
      mean_p <- md$p_value
    }
    pm <- stats_df[stats_df$metric == "permanova_F", ]
    if (nrow(pm) == 1) {
      perm_F <- pm$value
      perm_p <- pm$p_value
    }
  }

  # 1) Counts plot
  subtitle_counts <- NULL
  if (!is.na(mean_p)) {
    subtitle_counts <- paste0("mean diff=", mean_diff, ", perm p=", mean_p)
  }

  p_counts <- ggplot(summary_df, aes(x = group, y = lps_gene_count, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
    geom_boxplot(width = 0.18, outlier.shape = NA, color = "#333333") +
    geom_jitter(width = 0.12, size = 1.6, alpha = 0.6, color = "#222222") +
    scale_fill_manual(values = c("Earth" = "#d9a441", "ISS" = "#4f7cac")) +
    labs(title = paste0(label, ": LPS gene diversity"),
         subtitle = subtitle_counts,
         x = NULL, y = "LPS gene families per genome") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))

  ggsave(filename = file.path(out_dir, output_name(tag, "counts_by_group")),
         plot = p_counts, width = 8.5, height = 5.5, dpi = 200)

  # 2) MDS plot
  jac_df <- read.delim(jaccard_path, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(jac_df) <- jac_df[[1]]
  jac_mat <- as.matrix(jac_df[, -1, drop = FALSE])

  mds <- cmdscale(as.dist(jac_mat), k = 2)
  mds_df <- data.frame(isolate = rownames(mds),
                       dim1 = mds[, 1],
                       dim2 = mds[, 2],
                       stringsAsFactors = FALSE)
  mds_df <- merge(mds_df, summary_df[, c("isolate", "group")], by = "isolate", all.x = TRUE)

  # label top outliers to keep plot readable
  center <- colMeans(mds_df[, c("dim1", "dim2")])
  mds_df$dist <- sqrt((mds_df$dim1 - center[1])^2 + (mds_df$dim2 - center[2])^2)
  label_df <- mds_df[order(-mds_df$dist), ]
  label_df <- head(label_df, 12)

  subtitle_mds <- NULL
  if (!is.na(perm_p)) {
    subtitle_mds <- paste0("PERMANOVA F=", perm_F, ", perm p=", perm_p)
  }

  set.seed(7)
  jitter_x <- diff(range(mds_df$dim1)) * 0.01
  jitter_y <- diff(range(mds_df$dim2)) * 0.01
  mds_df$jdim1 <- mds_df$dim1 + rnorm(nrow(mds_df), 0, jitter_x)
  mds_df$jdim2 <- mds_df$dim2 + rnorm(nrow(mds_df), 0, jitter_y)

  p_mds <- ggplot(mds_df, aes(x = jdim1, y = jdim2, color = group)) +
    geom_point(size = 1.8, alpha = 0.7) +
    geom_text_repel(data = label_df, aes(x = dim1, y = dim2, label = isolate),
                    size = 2.5, max.overlaps = Inf, box.padding = 0.2,
                    min.segment.length = 0, segment.color = "#666666") +
    scale_color_manual(values = c("Earth" = "#d9a441", "ISS" = "#4f7cac")) +
    labs(title = paste0(label, ": LPS gene content (MDS, Jaccard)"),
         subtitle = subtitle_mds,
         x = "MDS1", y = "MDS2") +
    theme_classic(base_size = 12) +
    theme(legend.position = "top",
          plot.title = element_text(face = "bold"))

  ggsave(filename = file.path(out_dir, output_name(tag, "mds_jaccard")),
         plot = p_mds, width = 8.5, height = 6, dpi = 200)

  # 3) Heatmap (top genes by q-value, fallback to abs diff)
  prev_df <- read.delim(prevalence_path, stringsAsFactors = FALSE)
  fisher_df <- NULL
  if (file.exists(fisher_path)) {
    fisher_df <- read.delim(fisher_path, stringsAsFactors = FALSE)
  }

  top_n <- 30
  top_genes <- NULL
  title_suffix <- ""
  if (!is.null(fisher_df)) {
    fisher_df <- fisher_df[order(fisher_df$q_value, fisher_df$p_value), ]
    sig <- fisher_df[fisher_df$q_value <= 0.05, ]
    if (nrow(sig) > 0) {
      top_genes <- sig$gene_id[seq_len(min(top_n, nrow(sig)))]
      title_suffix <- " (q<=0.05)"
    }
  }
  if (is.null(top_genes)) {
    prev_df <- prev_df[order(-prev_df$abs_diff), ]
    top_n <- min(top_n, nrow(prev_df))
    top_genes <- prev_df$gene_id[seq_len(top_n)]
    title_suffix <- " (top abs diff)"
  }

  pres_df <- read.delim(presence_path, check.names = FALSE, stringsAsFactors = FALSE)
  names(pres_df) <- apply_rename_vector(names(pres_df))
  rownames(pres_df) <- pres_df$gene_id
  pres_mat <- as.matrix(pres_df[top_genes, -(1:3), drop = FALSE])

  iso_order <- summary_df$isolate[order(summary_df$group, summary_df$isolate)]
  if (!is.null(metadata_df)) {
    iso_order <- metadata_order(colnames(pres_mat), metadata_df)
  }
  pres_mat <- pres_mat[, iso_order, drop = FALSE]

  annotation_col <- data.frame(group = summary_df$group[match(iso_order, summary_df$isolate)])
  rownames(annotation_col) <- iso_order

  ann_colors <- list(group = c("Earth" = "#d9a441", "ISS" = "#4f7cac"))

  pheatmap(pres_mat,
           color = c("#f0f0f0", "#1b7837"),
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           fontsize_row = 9,
           fontsize_col = 6,
           main = paste0(label, ": Top LPS gene differences", title_suffix),
           filename = file.path(out_dir, output_name(tag, "heatmap_top_diff")),
           width = 12, height = 7)

  # 4) O-antigen focused heatmap (clean, presentation-ready)
  oantigen_re <- "^(rfb|wbb|wzm|wzt|waaL|rfaZ)"
  o_rows <- rownames(pres_df)[grepl(oantigen_re, rownames(pres_df), ignore.case = TRUE)]
  if (length(o_rows) > 0) {
    o_mat <- as.matrix(pres_df[o_rows, -(1:3), drop = FALSE])
    o_mat <- o_mat[, iso_order, drop = FALSE]
    pheatmap(o_mat,
             color = c("#f0f0f0", "#1b7837"),
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             fontsize_row = 9,
             fontsize_col = 6,
             main = if (remove_oantigen_title) "" else paste0(label, ": O-antigen genes (presence/absence)"),
             filename = file.path(out_dir, output_name(tag, "heatmap_oantigen")),
             width = 12, height = 7)
  }

  # 5) Volcano-style plot of gene association (q-values)
  if (!is.null(fisher_df)) {
    fisher_df$odds_ratio <- as.numeric(fisher_df$odds_ratio)
    fisher_df$q_value <- as.numeric(fisher_df$q_value)
    fisher_df$log2_or <- log2(fisher_df$odds_ratio)
    fisher_df$log2_or[is.infinite(fisher_df$log2_or)] <- 6
    fisher_df$log2_or[fisher_df$odds_ratio == 0] <- -6
    fisher_df$neglog10q <- -log10(fisher_df$q_value + 1e-12)
    fisher_df$class <- ifelse(fisher_df$log2_or > 0, "ISS-enriched", "Earth-enriched")
    fisher_df$label <- ifelse(fisher_df$q_value <= 0.01, fisher_df$gene_id, NA)
    fisher_df$is_oantigen <- grepl(oantigen_re, fisher_df$gene_id, ignore.case = TRUE)

    p_volcano <- ggplot(fisher_df, aes(x = log2_or, y = neglog10q)) +
      geom_point(aes(color = class, alpha = is_oantigen), size = 2) +
      scale_color_manual(values = c("ISS-enriched" = "#4f7cac", "Earth-enriched" = "#d95f02")) +
      scale_alpha_manual(values = c("TRUE" = 1.0, "FALSE" = 0.4), guide = "none") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#888888") +
      geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "#888888") +
      geom_text_repel(aes(label = label), size = 3, max.overlaps = 30) +
      labs(title = paste0(label, ": Gene association (ISS vs Earth)"),
           subtitle = "log2(odds ratio) vs -log10(q)",
           x = "log2(odds ratio)", y = "-log10(q)") +
      theme_classic(base_size = 12) +
      theme(legend.position = "top",
            plot.title = element_text(face = "bold"))

    ggsave(filename = file.path(out_dir, output_name(tag, "volcano")),
           plot = p_volcano, width = 8.5, height = 6, dpi = 200)
  }
}

plot_set("lps_strict", "All isolates (strict)")
plot_set("lps_piersonii_strict", "P. piersonii only (strict)")

message("Plots written to: ", out_dir)
