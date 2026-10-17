# QC filtering of the post-CellBender merged Seurat object (24h, 48h and 72h blastoids)
#
# Input:  3.Results/2.CellBender/merged_filtered_seurat.rds (from h5_to_seurat.R)
# Steps:  1. Remove the "_cellbender" suffix from sample names, cell barcodes and layers (one counts layer per sample, counts.<Sample>)
#         2. Compute QC metrics per sample
#         3. Per-sample filtering: hard bounds + MAD-based bounds
#         4. Per-sample QC histograms (kept vs removed) and cross-sample QC plots and figures before and after filtering
# Output: 3.Results/3.QC_filtering/
#           merged_postfilt_seurat.rds
#           QC_figure_prefilt.pdf, QC_figure_postfilt.pdf
#           QC_plots_combined_prefilt.pdf, QC_plots_combined_postfilt.pdf
#           Thresholding_combined.txt, qc_filtering_summary_mqc.yaml
#           per_sample/<Sample>/<Sample>_{QC_histograms.pdf, thresholds.txt, cells_to_keep.rds, stats.rds}
#
# Usage:  Rscript 2.Scripts/QC.R   (from any directory)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(ggridges)
  library(patchwork)
})

# Repository root: the folder above 2.Scripts (from the script's own path with Rscript or source(),
# otherwise the working directory)
Root <- local({
  f <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), fixed = TRUE)
  if (length(f) == 0) f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (length(f) == 1) normalizePath(file.path(dirname(f), "..")) else getwd()
})
if (!dir.exists(file.path(Root, "2.Scripts"))) stop("Repository root not found (no 2.Scripts folder in ", Root, ")")

# Results folder: 3.Results, or the folder given in SC_RESULTS (relative to the repository root unless absolute)
Results <- Sys.getenv("SC_RESULTS", "3.Results")
if (!grepl("^(/|[A-Za-z]:)", Results)) Results <- file.path(Root, Results)
source(file.path(Root, "2.Scripts", "functions.R"))

Input_rds <- file.path(Results, "2.CellBender", "merged_filtered_seurat.rds")
Out_dir <- file.path(Results, "3.QC_filtering")
dir.create(file.path(Out_dir, "per_sample"), recursive = TRUE, showWarnings = FALSE)

# Thresholds (NULL switches a filter off)
nFeatures_lower <- 2000
nFeatures_upper <- 10000
nCount_upper <- 50000
mito_percent_hard <- 2.5
pct_counts_mt_nmads <- NULL
log1p_total_counts_nmads <- 5
log1p_n_genes_by_counts_nmads <- 5
top_n_genes <- 20
pct_counts_top20_nmads <- 5
mad_constant <- 1.4826

t0 <- Sys.time()


# Load data ----
log_msg("Reading ", Input_rds)
merged_df <- readRDS(Input_rds)

# Strip the "_cellbender" suffix from sample names, cell barcodes and layer names
merged_df$Sample <- sub("_cellbender$", "", merged_df$Sample)
merged_df$orig.ident <- factor(sub("_cellbender$", "", merged_df$orig.ident))
merged_df <- RenameCells(merged_df, new.names = sub("_cellbender_", "_", colnames(merged_df)))
merged_df[["RNA"]] <- JoinLayers(merged_df[["RNA"]])
merged_df[["RNA"]] <- split(merged_df[["RNA"]], f = merged_df$Sample)
Idents(merged_df) <- "Sample"
samples_list <- unique(merged_df$Sample)
n_genes_input <- nrow(merged_df)
log_msg("Input: ", ncol(merged_df), " cells, ", n_genes_input, " genes; samples: ", paste(samples_list, collapse = ", "))


# Compute QC metrics ----
top_col <- paste0("pct_counts_in_top_", top_n_genes, "_genes")
metrics <- do.call(rbind, lapply(samples_list, function(s) {
  layer <- paste0("counts.", s)
  log_msg("Computing QC metrics: ", s)
  m <- qc_metrics(LayerData(merged_df, assay = "RNA", layer = layer), top_n_genes)
  m[intersect(colnames(merged_df)[merged_df$Sample == s], rownames(m)), ]
}))
names(metrics)[names(metrics) == "pct_top"] <- top_col
metrics$log1p_total_counts <- log1p(metrics$nCount_RNA)
metrics$log1p_n_genes_by_counts <- log1p(metrics$nFeature_RNA)
stopifnot(setequal(rownames(metrics), colnames(merged_df)))
merged_df <- AddMetaData(merged_df, metadata = metrics[colnames(merged_df), ])
for (v in c("pct_counts_ribo", "pct_counts_hb")) {
  if (!has_variation(merged_df[[v]][, 1])) log_msg(v, " is constant (", unique(merged_df[[v]][, 1]), ") in every cell")
}


# Per-sample filtering ----
sample_stats <- list()
all_cells_to_keep <- character(0)

for (s in samples_list) {
  s_dir <- file.path(Out_dir, "per_sample", s)
  dir.create(s_dir, recursive = TRUE, showWarnings = FALSE)
  md <- merged_df[[]][merged_df$Sample == s, ]

  bounds_nFeature <- hard_bounds(nFeatures_lower, nFeatures_upper)
  bounds_mito_mad <- mad_bounds(md$pct_counts_mt, pct_counts_mt_nmads, mad_constant)
  bounds_log1p_counts <- mad_bounds(md$log1p_total_counts, log1p_total_counts_nmads, mad_constant)
  bounds_log1p_genes <- mad_bounds(md$log1p_n_genes_by_counts, log1p_n_genes_by_counts_nmads, mad_constant)
  bounds_top20 <- mad_bounds(md[[top_col]], pct_counts_top20_nmads, mad_constant)

  pass <- list(
    nFeature_RNA = in_bounds(md$nFeature_RNA, bounds_nFeature),
    nCount_RNA_cap = below_cap(md$nCount_RNA, nCount_upper),
    pct_counts_mt_cap = below_cap(md$pct_counts_mt, mito_percent_hard),
    pct_counts_mt_mad = in_bounds(md$pct_counts_mt, bounds_mito_mad),
    log1p_total_counts_mad = in_bounds(md$log1p_total_counts, bounds_log1p_counts),
    log1p_n_genes_by_counts_mad = in_bounds(md$log1p_n_genes_by_counts, bounds_log1p_genes),
    pct_counts_top_mad = in_bounds(md[[top_col]], bounds_top20)
  )
  sample_keep <- Reduce(`&`, pass)
  md$qc_status <- ifelse(sample_keep, "Kept", "Removed")
  log_msg(s, ": ", sum(sample_keep), " of ", nrow(md), " cells kept")

  # Thresholds report
  tf <- file.path(s_dir, paste0(s, "_thresholds.txt"))
  sink(tf)
  cat("==================================================\n")
  cat("Per-Sample QC Thresholding\n")
  cat("Sample:", s, "\n")
  cat("Generated on:", as.character(Sys.time()), "\n")
  cat("==================================================\n")
  cat("Cells in sample before filtering:", nrow(md), "\n")
  cat("\nThresholds applied:\n")
  cat("  nFeature_RNA hard bounds  :", nFeatures_lower, "<  x  <", nFeatures_upper, "\n")
  cat("  nCount_RNA hard cap       :", if (is.null(nCount_upper)) "off" else paste("<=", nCount_upper), "\n")
  cat("  pct_counts_mt hard cap    :", if (is.null(mito_percent_hard)) "off" else paste("<=", mito_percent_hard), "\n")
  cat("  pct_counts_mt MAD nmads   :", or_default(pct_counts_mt_nmads, "off"), "\n")
  cat("  log1p_total_counts nmads  :", or_default(log1p_total_counts_nmads, "off"),
      sprintf("(%.3f < x < %.3f; nCount_RNA %.0f-%.0f)", bounds_log1p_counts["lower"], bounds_log1p_counts["upper"],
              expm1(bounds_log1p_counts["lower"]), expm1(bounds_log1p_counts["upper"])), "\n")
  cat("  log1p_n_genes_by_counts nmads:", or_default(log1p_n_genes_by_counts_nmads, "off"),
      sprintf("(%.3f < x < %.3f; nFeature_RNA %.0f-%.0f)", bounds_log1p_genes["lower"], bounds_log1p_genes["upper"],
              expm1(bounds_log1p_genes["lower"]), expm1(bounds_log1p_genes["upper"])), "\n")
  cat("  pct_counts_top", top_n_genes, "_genes nmads:", or_default(pct_counts_top20_nmads, "off"), "\n")
  cat("\nCells removed by each filter (a cell can fail several):\n")
  for (nm in names(pass)) cat(sprintf("  %-28s: %d\n", nm, sum(!pass[[nm]])))
  cat("\nPer-filter pass counts:\n")
  cat("  Combined (final) pass     :", sum(sample_keep), "/", nrow(md), "(", round(100 * sum(sample_keep) / nrow(md), 1), "%)\n")
  sink()

  # Histograms: kept (green) vs removed (red)
  hist_plots <- list(
    nFeature = qc_histogram(md, "nFeature_RNA",
      vlines = list(lower = list(x = bounds_nFeature["lower"], color = "steelblue"),
                    upper = list(x = bounds_nFeature["upper"], color = "darkorange")),
      title = paste(s, "- nFeature_RNA"), xlab = "nFeature_RNA"),
    nCount = qc_histogram(md, "nCount_RNA",
      vlines = list(hard_upper = list(x = or_default(nCount_upper, Inf), color = "firebrick")),
      title = paste(s, "- nCount_RNA"), xlab = "nCount_RNA"),
    mito = qc_histogram(md, "pct_counts_mt",
      vlines = list(hard_upper = list(x = or_default(mito_percent_hard, Inf), color = "firebrick"),
                    mad_lower = list(x = bounds_mito_mad["lower"], color = "steelblue"),
                    mad_upper = list(x = bounds_mito_mad["upper"], color = "darkorange")),
      title = paste(s, "- pct_counts_mt"), xlab = "pct_counts_mt", xlim = c(0, 10)),
    hb = qc_histogram(md, "pct_counts_hb", title = paste(s, "- pct_counts_hb"), xlab = "pct_counts_hb", xlim = c(0, 0.5)),
    ribo = qc_histogram(md, "pct_counts_ribo", title = paste(s, "- pct_counts_ribo"), xlab = "pct_counts_ribo"),
    log1p_counts = qc_histogram(md, "log1p_total_counts",
      vlines = list(lower = list(x = bounds_log1p_counts["lower"], color = "steelblue"),
                    upper = list(x = bounds_log1p_counts["upper"], color = "darkorange")),
      title = paste(s, "- log1p_total_counts"), xlab = "log1p(nCount_RNA)"),
    log1p_genes = qc_histogram(md, "log1p_n_genes_by_counts",
      vlines = list(lower = list(x = bounds_log1p_genes["lower"], color = "steelblue"),
                    upper = list(x = bounds_log1p_genes["upper"], color = "darkorange")),
      title = paste(s, "- log1p_n_genes_by_counts"), xlab = "log1p(nFeature_RNA)"),
    top20 = qc_histogram(md, top_col,
      vlines = list(lower = list(x = bounds_top20["lower"], color = "steelblue"),
                    upper = list(x = bounds_top20["upper"], color = "darkorange")),
      title = paste0(s, " - pct_counts_in_top_", top_n_genes, "_genes"), xlab = top_col)
  )
  pdf(file.path(s_dir, paste0(s, "_QC_histograms.pdf")), width = 17, height = 11)
  print(wrap_plots(hist_plots, ncol = 3) +
    plot_annotation(title = paste("QC Histograms -", s),
                    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))))
  invisible(dev.off())

  cells_to_keep <- rownames(md)[sample_keep]
  saveRDS(cells_to_keep, file.path(s_dir, paste0(s, "_cells_to_keep.rds")))
  all_cells_to_keep <- c(all_cells_to_keep, cells_to_keep)

  sample_stats[[s]] <- list(
    sample_name = s,
    n_cells_before = nrow(md),
    n_cells_after = sum(sample_keep),
    pct_retained = round(100 * sum(sample_keep) / nrow(md), 2),
    n_removed_by_filter = vapply(pass, function(x) sum(!x), integer(1)),
    thresholds = list(
      nFeatures_lower = nFeatures_lower,
      nFeatures_upper = nFeatures_upper,
      nCount_upper = nCount_upper,
      mito_percent_hard = mito_percent_hard,
      pct_counts_mt_nmads = pct_counts_mt_nmads,
      log1p_total_counts_nmads = log1p_total_counts_nmads,
      log1p_n_genes_by_counts_nmads = log1p_n_genes_by_counts_nmads,
      pct_counts_top20_nmads = pct_counts_top20_nmads
    )
  )
  saveRDS(sample_stats[[s]], file.path(s_dir, paste0(s, "_stats.rds")))
}


# Merge and cross-sample plots ----
n_before <- ncol(merged_df)
n_after <- length(all_cells_to_keep)

log_msg("Writing pre-filter cross-sample QC plots")
qc_pre <- qc_combined_plots(merged_df)
write_combined_qc_pdf(qc_pre, file.path(Out_dir, "QC_plots_combined_prefilt.pdf"), "QC before filtering", n_before, n_before)
write_qc_figure(qc_pre, file.path(Out_dir, "QC_figure_prefilt.pdf"))
fits_pre <- qc_pre$fits
rm(qc_pre)

merged_df_filt <- subset(merged_df, cells = all_cells_to_keep)
rm(merged_df); invisible(gc())

log_msg("Writing post-filter cross-sample QC plots")
qc_post <- qc_combined_plots(merged_df_filt)
write_combined_qc_pdf(qc_post, file.path(Out_dir, "QC_plots_combined_postfilt.pdf"), "QC after filtering", n_before, n_after)
write_qc_figure(qc_post, file.path(Out_dir, "QC_figure_postfilt.pdf"))
fits_post <- qc_post$fits
rm(qc_post)


# Combined thresholding report ----
sink(file.path(Out_dir, "Thresholding_combined.txt"))
cat("==================================================\n")
cat("COMBINED QC THRESHOLDING REPORT\n")
cat("Generated on:", as.character(Sys.time()), "\n")
cat("==================================================\n")
cat("\nInput:", Input_rds, "\n")
cat("Before filtering:", n_before, "cells,", n_genes_input, "genes\n")
cat("Samples:", paste(samples_list, collapse = ", "), "\n")
cat("\n==================================================\n")
cat("PER-SAMPLE RESULTS\n")
cat("==================================================\n")
for (st in sample_stats) {
  cat("\n", st$sample_name, ":\n")
  cat("  Cells before: ", st$n_cells_before, "\n")
  cat("  Cells after:  ", st$n_cells_after, "\n")
  cat("  Retained:     ", st$pct_retained, "%\n")
  cat("  Removed by filter:", paste(names(st$n_removed_by_filter), st$n_removed_by_filter, sep = " = ", collapse = ", "), "\n")
}
cat("\n==================================================\n")
cat("OVERALL SUMMARY\n")
cat("==================================================\n")
cat("Total cells before filtering:", n_before, "\n")
cat("Total cells after filtering :", n_after, "\n")
cat("Retained:", round(100 * n_after / n_before, 1), "% of cells overall\n")
cat("\nSaturation fit, log10(nFeature_RNA) = Asym + (R0 - Asym) * exp(-exp(lrc) * log10(nCount_RNA)):\n")
for (stage in c("prefilt", "postfilt")) {
  fits <- if (stage == "prefilt") fits_pre else fits_post
  for (s in names(fits)) {
    cf <- fits[[s]]
    cat(sprintf("  %-8s %-26s Asym = %.3f, R0 = %.3f, lrc = %.3f\n", stage, s, cf["Asym"], cf["R0"], cf["lrc"]))
  }
}
cat("\nSaved filtered object to: merged_postfilt_seurat.rds\n")
cat("Combined PDFs written to: QC_plots_combined_prefilt.pdf, QC_plots_combined_postfilt.pdf\n")
cat("QC figures written to: QC_figure_prefilt.pdf, QC_figure_postfilt.pdf\n")
sink()


# MultiQC table ----
mqc <- c(
  "id: 'qc_filtering_summary'",
  "section_name: 'QC Filtering Summary'",
  "description: 'Per-sample QC filtering results with MAD-based thresholds'",
  "plot_type: 'table'",
  "pconfig:",
  "  id: 'qc_filtering_table'",
  "  title: 'QC Filtering Results'",
  "",
  "headers:",
  "  Cells_Before:",
  "    title: 'Cells Before'",
  "    description: 'Number of cells before QC filtering'",
  "    format: '{:,.0f}'",
  "  Cells_After:",
  "    title: 'Cells After'",
  "    description: 'Number of cells after QC filtering'",
  "    format: '{:,.0f}'",
  "  Pct_Retained:",
  "    title: '% Retained'",
  "    description: 'Percentage of cells retained after filtering'",
  "    suffix: '%'",
  "    format: '{:,.1f}'",
  "    scale: 'RdYlGn'",
  "    min: 50",
  "    max: 100",
  "",
  "data:",
  unlist(lapply(sample_stats, function(st) c(
    paste0("  ", st$sample_name, ":"),
    paste0("    Cells_Before: ", st$n_cells_before),
    paste0("    Cells_After: ", st$n_cells_after),
    paste0("    Pct_Retained: ", round(st$pct_retained, 1))
  )))
)
writeLines(mqc, file.path(Out_dir, "qc_filtering_summary_mqc.yaml"))


# Save filtered object ----
log_msg("Saving ", file.path(Out_dir, "merged_postfilt_seurat.rds"))
saveRDS(merged_df_filt, file.path(Out_dir, "merged_postfilt_seurat.rds"))
log_msg("Done: ", n_after, " of ", n_before, " cells kept (", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min)")
