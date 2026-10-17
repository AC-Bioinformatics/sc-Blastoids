#!/usr/bin/env Rscript
# ============================================================================
# Multi-Resolution Clustering Analysis
# ============================================================================
# Purpose: Cluster cells at multiple resolutions (0.1 to 1.5 in 0.1 increments, Louvain)
#          and create UMAP visualizations for comparison
# Input: 3.Results/5.Normalisation/merged_normalised_seurat.rds
# Output: 3.Results/6.Clustering/Louvain/
#           seurat_multi_res_clustered.rds (RNA_snn_res.<resolution> columns)
#           umap_resolution_<resolution>.pdf, umap_all_resolutions_combined.pdf, resolution_vs_clusters.pdf,
#           clustree_stability.pdf, cluster_summary.csv, clustering_summary.txt
# Usage: Rscript 2.Scripts/Clustering.R   (from any directory)
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

# ============================================================================
# Configuration
# ============================================================================

# Repository root: the folder above 2.Scripts (from the script's own path with Rscript or source(),
# otherwise the working directory)
Root <- local({
  f <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), fixed = TRUE)
  if (length(f) == 0) f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (length(f) == 1) normalizePath(file.path(dirname(f), "..")) else getwd()
})
if (!dir.exists(file.path(Root, "2.Scripts"))) stop("Repository root not found (no 2.Scripts folder in ", Root, ")")

# Input file path
input_file <- file.path(Root, "3.Results", "5.Normalisation", "merged_normalised_seurat.rds")

# Output directory
outdir <- file.path(Root, "3.Results", "6.Clustering", "Louvain")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Resolution parameters
resolutions <- seq(0.1, 1.5, by = 0.1)  # 0.1, 0.2, ..., 1.5

# Number of dimensions to use (based on PCA from normalisation)
n_dims <- 20  # Adjust this based on your elbow plot

# ============================================================================
# Load Data
# ============================================================================

message("========================================")
message("MULTI-RESOLUTION CLUSTERING ANALYSIS")
message("========================================\n")

message("Loading Seurat object from:")
message("  ", input_file)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

seurat_obj <- readRDS(input_file)

message("Loaded successfully")
message("  Cells: ", ncol(seurat_obj))
message("  Genes: ", nrow(seurat_obj))
message("  Assays: ", paste(Assays(seurat_obj), collapse = ", "))

# ============================================================================
# Check for PCA and UMAP
# ============================================================================

message("\nChecking for dimensionality reductions...")

# Check if PCA exists
if (!"pca" %in% names(seurat_obj@reductions)) {
  message(" PCA not found - running PCA now...")
  seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)
  message("PCA complete")
} else {
  message("PCA found")
}

# Check if UMAP exists
if (!"umap" %in% names(seurat_obj@reductions)) {
  message(" UMAP not found - running UMAP now...")
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_dims, verbose = FALSE)
  message("UMAP complete")
} else {
  message("UMAP found")
}

# ============================================================================
# Find Neighbors (only need to do this once)
# ============================================================================

message("\nFinding neighbors using ", n_dims, " dimensions...")
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_dims, verbose = FALSE) 
message("Neighbors computed")

# ============================================================================
# Cluster at Multiple Resolutions
# ============================================================================

message("\n========================================")
message("CLUSTERING AT MULTIPLE RESOLUTIONS")
message("========================================\n")

# Store clustering results and create plots
umap_plots <- list()
cluster_summary <- data.frame(
  Resolution = numeric(),
  N_Clusters = integer(),
  stringsAsFactors = FALSE
)

for (res in resolutions) {
  message(sprintf("Processing resolution %.1f...", res))

  # Run clustering
  seurat_obj <- FindClusters(seurat_obj,
                             resolution = res,
                             verbose = FALSE)

  # Store the cluster identities with resolution-specific name
  cluster_col_name <- paste0("RNA_snn_res.", res)
  n_clusters <- length(unique(seurat_obj[[cluster_col_name]][,1]))

  message(sprintf("  Found %d clusters", n_clusters))

  # Add to summary
  cluster_summary <- rbind(cluster_summary,
                          data.frame(Resolution = res,
                                   N_Clusters = n_clusters))

  # Create UMAP plot for this resolution
  p <- DimPlot(seurat_obj,
               reduction = "umap",
               group.by = cluster_col_name,
               label = TRUE,
               label.size = 4,
               pt.size = 0.5) +
    ggtitle(sprintf("Resolution %.1f\n(%d clusters)", res, n_clusters)) +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          legend.position = "right")

  # Save individual plot
  ggsave(filename = file.path(outdir, sprintf("umap_resolution_%.1f.pdf", res)),
         plot = p,
         width = 10,
         height = 8)

  # Store plot for combined figure
  umap_plots[[as.character(res)]] <- p

  message(sprintf("  Saved umap_resolution_%.1f.pdf", res))
}

# ============================================================================
# Create Combined Figure
# ============================================================================

message("\n========================================")
message("CREATING COMBINED COMPARISON FIGURE")
message("========================================\n")

# Combine all plots into a single figure
combined_plot <- wrap_plots(umap_plots, ncol = 3) +
  plot_annotation(
    title = "Multi-Resolution Clustering Comparison",
    subtitle = sprintf("Resolutions from %.1f to %.1f (step: 0.1) | %d dimensions used",
                      min(resolutions), max(resolutions), n_dims),
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                 plot.subtitle = element_text(size = 12, hjust = 0.5))
  )

# Save combined plot
ggsave(filename = file.path(outdir, "umap_all_resolutions_combined.pdf"),
       plot = combined_plot,
       width = 18,
       height = 20)

message("Saved combined plot: umap_all_resolutions_combined.pdf")

# ============================================================================
# Create Resolution vs Cluster Count Plot
# ============================================================================

message("\nCreating resolution vs cluster count plot...")

p_resolution_curve <- ggplot(cluster_summary,
                            aes(x = Resolution, y = N_Clusters)) +
  geom_line(color = "steelblue", size = 1.2) +
  geom_point(color = "steelblue", size = 3) +
  geom_text(aes(label = N_Clusters), vjust = -1, size = 3.5) +
  scale_x_continuous(breaks = resolutions) +
  scale_y_continuous(breaks = seq(0, max(cluster_summary$N_Clusters) + 5, by = 2)) +
  labs(title = "Number of Clusters vs Resolution",
       x = "Resolution",
       y = "Number of Clusters") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = file.path(outdir, "resolution_vs_clusters.pdf"),
       plot = p_resolution_curve,
       width = 10,
       height = 6)

message("Saved resolution_vs_clusters.pdf")

# ============================================================================
# Save Summary Statistics
# ============================================================================

message("\nWriting summary statistics...")

# Save cluster summary CSV
write.csv(cluster_summary,
          file.path(outdir, "cluster_summary.csv"),
          row.names = FALSE)

# Save detailed summary text file
summary_file_path <- file.path(outdir, "clustering_summary.txt")
summary_file <- file(summary_file_path, open = "wt")

cat("=======================================================\n", file = summary_file)
cat("     MULTI-RESOLUTION CLUSTERING SUMMARY              \n", file = summary_file)
cat("=======================================================\n\n", file = summary_file)

cat("INPUT FILE:\n", file = summary_file)
cat("  ", input_file, "\n\n", file = summary_file)

cat("PARAMETERS:\n", file = summary_file)
cat("  Number of dimensions used:", n_dims, "\n", file = summary_file)
cat("  Resolution range: ", min(resolutions), " to ", max(resolutions), "\n", file = summary_file)
cat("  Resolution step: 0.1\n", file = summary_file)
cat("  Total resolutions tested:", length(resolutions), "\n\n", file = summary_file)

cat("DATA:\n", file = summary_file)
cat("  Total cells:", ncol(seurat_obj), "\n", file = summary_file)
cat("  Total genes:", nrow(seurat_obj), "\n\n", file = summary_file)

cat("CLUSTERING RESULTS:\n", file = summary_file)
cat("---------------------------------------------------\n", file = summary_file)
cat(sprintf("%-12s | %-15s\n", "Resolution", "Clusters Found"), file = summary_file)
cat("---------------------------------------------------\n", file = summary_file)
for (i in 1:nrow(cluster_summary)) {
  cat(sprintf("%-12.1f | %-15d\n",
              cluster_summary$Resolution[i],
              cluster_summary$N_Clusters[i]), file = summary_file)
}
cat("---------------------------------------------------\n\n", file = summary_file)

cat("CLUSTER COLUMN NAMES IN SEURAT OBJECT:\n", file = summary_file)
for (res in resolutions) {
  cat("  RNA_snn_res.", res, "\n", sep = "", file = summary_file)
}

cat("\n=======================================================\n", file = summary_file)
cat("OUTPUT FILES CREATED:\n", file = summary_file)
cat("  → cluster_summary.csv\n", file = summary_file)
cat("  → resolution_vs_clusters.pdf\n", file = summary_file)
cat("  → umap_all_resolutions_combined.pdf\n", file = summary_file)
for (res in resolutions) {
  cat(sprintf("  → umap_resolution_%.1f.pdf\n", res), file = summary_file)
}

cat("\n=======================================================\n", file = summary_file)
cat("NEXT STEPS:\n", file = summary_file)
cat("  1. Review the combined UMAP plot to choose optimal resolution\n", file = summary_file)
cat("  2. Set the default resolution in your Seurat object:\n", file = summary_file)
cat("     Idents(seurat_obj) <- seurat_obj$RNA_snn_res.X.X\n", file = summary_file)
cat("  3. Find marker genes for the chosen resolution:\n", file = summary_file)
cat("     markers <- FindAllMarkers(seurat_obj)\n", file = summary_file)
cat("  4. Save the updated object:\n", file = summary_file)
cat("     saveRDS(seurat_obj, 'seurat_multi_res_clustered.rds')\n", file = summary_file)
cat("=======================================================\n", file = summary_file)

# Close the file explicitly
close(summary_file)
message("Saved: clustering_summary.txt")

# ============================================================================
# Save Updated Seurat Object
# ============================================================================

message("\nSaving updated Seurat object with all resolution clustering...")
saveRDS(seurat_obj, file.path(outdir, "seurat_multi_res_clustered.rds"))
message("Saved: seurat_multi_res_clustered.rds")

# ============================================================================
# Final Summary
# ============================================================================

message("\n========================================")
message("MULTI-RESOLUTION CLUSTERING COMPLETE!")
message("========================================\n")

message("SUMMARY:")
message(sprintf("   • Tested %d resolutions (%.1f to %.1f)",
               length(resolutions), min(resolutions), max(resolutions)))
message(sprintf("   • Cluster count range: %d to %d",
               min(cluster_summary$N_Clusters), max(cluster_summary$N_Clusters)))
message(sprintf("   • Analyzed %d cells", ncol(seurat_obj)))

message("\nOUTPUT DIRECTORY: ", normalizePath(outdir))

message("\nFILES CREATED:")
message("   → umap_all_resolutions_combined.pdf (main comparison plot)")
message("   → resolution_vs_clusters.pdf (cluster count curve)")
message("   → cluster_summary.csv (summary table)")
message("   → clustering_summary.txt (detailed summary)")
message("   → seurat_multi_res_clustered.rds (updated Seurat object)")
message(sprintf("   → %d individual UMAP PDFs", length(resolutions)))

message("\nCLUSTER COUNTS BY RESOLUTION:")
for (i in 1:nrow(cluster_summary)) {
  message(sprintf("   Resolution %.1f: %d clusters",
                 cluster_summary$Resolution[i],
                 cluster_summary$N_Clusters[i]))
}

message("\nTIP: Review umap_all_resolutions_combined.pdf to choose")
message("   the optimal resolution for your analysis.")
message("\n========================================\n")

message("Analysis complete!")

# ============================================================================
# Clustree Analysis - Visualize Cluster Stability Across Resolutions
# ============================================================================

message("\n========================================")
message("CLUSTREE ANALYSIS - CLUSTER STABILITY")
message("========================================\n")

library(clustree)

message("Creating clustree visualization...")

# Basic clustree plot showing cluster relationships across resolutions
p_clustree <- clustree(seurat_obj,
                       prefix = "RNA_snn_res.",
                       node_colour = "sc3_stability",
                       node_size = 4,
                       edge_width = 1.5) +
  labs(title = "Cluster Stability Across Resolutions",
       subtitle = "Node size = cluster size | Color = SC3 stability score") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10))

ggsave(filename = file.path(outdir, "clustree_stability.pdf"),
       plot = p_clustree,
       width = 12,
       height = 10)

message("Saved: clustree_stability.pdf")
