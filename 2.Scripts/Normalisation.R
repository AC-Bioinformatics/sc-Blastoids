# Log normalisation, PCA and UMAP of the QC-filtered object (24h, 48h and 72h blastoids)
#
# Input:  3.Results/3.QC_filtering/merged_postfilt_seurat.rds (RNA counts split into one layer per sample)
# Steps:  1. LogNormalize per sample layer (scale factor 10,000)
#         2. 3,000 variable features (vst) per sample layer, combined into 3,000 consensus features
#         3. ScaleData on the consensus features (no regression)
#         4. PCA with 50 PCs; UMAP on PCs 1-20 (cosine)
#         5. Elbow plot, variance explained per PC (split into between- and within-sample), PCA and UMAP coloured by sample,
#            PCA with each sample highlighted in turn and UMAP of each sample on its own
# Output: 3.Results/5.Normalisation/
#           merged_normalised_seurat.rds
#           01_elbow_plot.{pdf,png}, 02_variance_explained_per_PC.{pdf,png}
#           03_pca_plot.{pdf,png}, 03b_pca_plot_by_sample.{pdf,png}, 04_umap_plot.{pdf,png}, 04b_umap_plot_by_sample.{pdf,png}
#           variable_features.txt, variable_features_by_sample.csv, variance_explained_per_PC.txt,
#           Normalisation_log.txt, normalisation_mqc.yaml
#
# Usage:  Rscript 2.Scripts/Normalisation.R   (from the repository root)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

Root <- getwd()
if (!file.exists(file.path(Root, "2.Scripts", "functions.R"))) stop("Run this script from the repository root")
source(file.path(Root, "2.Scripts", "functions.R"))

Input_rds <- file.path(Root, "3.Results", "3.QC_filtering", "merged_postfilt_seurat.rds")
Out_dir <- file.path(Root, "3.Results", "5.Normalisation")
dir.create(Out_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters
norm_scale_factor <- 10000
norm_nfeatures <- 3000
norm_npcs <- 50
umap_dims <- 20
seed <- 42

t0 <- Sys.time()


# Load ----
log_msg("Reading ", Input_rds)
merged_df <- readRDS(Input_rds)
count_layers <- Layers(merged_df[["RNA"]], search = "counts")
if (length(count_layers) == 1) {
  merged_df[["RNA"]] <- split(merged_df[["RNA"]], f = merged_df$Sample)
  count_layers <- Layers(merged_df[["RNA"]], search = "counts")
}
log_msg("Cells: ", ncol(merged_df), ", genes: ", nrow(merged_df), ", count layers: ", paste(count_layers, collapse = ", "))


# Normalisation and variable features ----
log_msg("LogNormalize per sample layer")
merged_df <- NormalizeData(merged_df, normalization.method = "LogNormalize", scale.factor = norm_scale_factor, verbose = FALSE)

log_msg("Variable features per sample layer and consensus")
merged_df <- FindVariableFeatures(merged_df, selection.method = "vst", nfeatures = norm_nfeatures, verbose = FALSE)
var_features <- VariableFeatures(merged_df)
writeLines(var_features, file.path(Out_dir, "variable_features.txt"))

data_layers <- Layers(merged_df[["RNA"]], search = "data")
layer_features <- lapply(data_layers, function(l) head(VariableFeatures(merged_df, layer = l), norm_nfeatures))
names(layer_features) <- sub("^data\\.", "", data_layers)
max_len <- max(lengths(layer_features))
write.csv(as.data.frame(lapply(layer_features, function(x) c(x, rep(NA, max_len - length(x)))), check.names = FALSE),
          file.path(Out_dir, "variable_features_by_sample.csv"), row.names = FALSE)

log_msg("Scaling the consensus variable features")
merged_df <- ScaleData(merged_df, features = var_features, verbose = FALSE)


# PCA and UMAP ----
log_msg("Running PCA with ", norm_npcs, " PCs")
merged_df <- RunPCA(merged_df, features = var_features, npcs = norm_npcs, seed.use = seed, verbose = FALSE)

log_msg("Running UMAP on PCs 1-", umap_dims)
merged_df <- RunUMAP(merged_df, dims = 1:umap_dims, metric = "cosine", seed.use = seed, verbose = FALSE)

stdev <- Stdev(merged_df, "pca")
pc_var <- stdev^2
scaled_var <- scaled_total_variance(merged_df)
n_pca_features <- scaled_var$n_features

# Share of each PC's variance explained by sample: R^2 of the PC scores on Sample
sample_r2 <- pc_r2(Embeddings(merged_df, "pca"), model.matrix(~ Sample, data = merged_df@meta.data)[, -1])

var_df <- data.frame(
  PC = seq_len(norm_npcs),
  stdev = stdev,
  pct_total_variance = 100 * pc_var / scaled_var$total,
  pct_of_PCs = 100 * pc_var / sum(pc_var),
  sample_r2 = unname(sample_r2)
)
var_df$cumulative_pct_total <- cumsum(var_df$pct_total_variance)
var_df$pct_total_sample <- var_df$pct_total_variance * var_df$sample_r2
write.table(var_df, file.path(Out_dir, "variance_explained_per_PC.txt"), sep = "\t", quote = FALSE, row.names = FALSE)


# Plots ----
title_theme <- theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
                     plot.subtitle = element_text(hjust = 0.5, size = 12))

# Elbow plot with the number of PCs used downstream marked
p_elbow <- ElbowPlot(merged_df, ndims = norm_npcs) +
  ggtitle("Elbow Plot: Standard Deviation by PC") +
  labs(subtitle = sprintf("PCs 1-%d: %.1f%% of total variance (%.1f%% of the variance in all %d PCs)",
                          umap_dims, var_df$cumulative_pct_total[umap_dims],
                          100 * sum(pc_var[1:umap_dims]) / sum(pc_var), norm_npcs)) +
  theme_classic() + title_theme +
  theme(plot.subtitle = element_text(colour = "red")) +
  geom_vline(xintercept = umap_dims, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = umap_dims + 1, y = max(stdev) * 0.6, hjust = 0, color = "red", size = 5,
           label = sprintf("PCs 1-%d: %.1f%% of total variance", umap_dims, var_df$cumulative_pct_total[umap_dims]))
ggsave(file.path(Out_dir, "01_elbow_plot.pdf"), plot = p_elbow, width = 10, height = 6)
ggsave(file.path(Out_dir, "01_elbow_plot.png"), plot = p_elbow, width = 10, height = 6, dpi = 150)

# Variance explained per PC (% of the total variance of the scaled variable features), split by sample
sample_df <- rbind(
  data.frame(PC = var_df$PC, part = "Between samples", pct = var_df$pct_total_sample),
  data.frame(PC = var_df$PC, part = "Within samples", pct = var_df$pct_total_variance - var_df$pct_total_sample)
)
sample_df$part <- factor(sample_df$part, levels = c("Within samples", "Between samples"))
p_var <- ggplot(sample_df, aes(x = PC, y = pct, fill = part)) +
  geom_col() +
  geom_text(data = var_df[1:10, ], aes(x = PC, y = pct_total_variance, label = sprintf("%.1f", pct_total_variance)),
            inherit.aes = FALSE, vjust = -0.4, size = 3) +
  scale_fill_manual(values = c("Between samples" = "firebrick", "Within samples" = "steelblue"), name = NULL) +
  labs(title = "Variance explained by each PC",
       subtitle = sprintf(paste0("%% of the total variance of the %s scaled variable features (LogNormalize)\n",
                                 "Red: part of each PC explained by sample (PC ~ Sample), R^2 x PC variance"),
                          format(n_pca_features, big.mark = ",")),
       x = "Principal Component", y = "Variance explained (%)") +
  theme_classic() + title_theme +
  theme(legend.position = "inside", legend.position.inside = c(0.85, 0.75))
p_cumvar <- ggplot(var_df, aes(x = PC, y = cumulative_pct_total)) +
  geom_line(linewidth = 1.2, color = "steelblue") +
  geom_point(size = 2, color = "steelblue") +
  geom_vline(xintercept = umap_dims, linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(title = "Cumulative variance explained",
       subtitle = sprintf("PCs 1-%d: %.1f%%   |   PCs 1-%d: %.1f%%", umap_dims, var_df$cumulative_pct_total[umap_dims],
                          norm_npcs, var_df$cumulative_pct_total[norm_npcs]),
       x = "Principal Component", y = "Cumulative variance explained (%)") +
  theme_classic() + title_theme
p_variance <- p_var / (p_cumvar | p_elbow)
ggsave(file.path(Out_dir, "02_variance_explained_per_PC.pdf"), plot = p_variance, width = 18, height = 12)
ggsave(file.path(Out_dir, "02_variance_explained_per_PC.png"), plot = p_variance, width = 18, height = 12, dpi = 150)

# PCA and UMAP coloured by sample
pca_plot <- DimPlot(merged_df, reduction = "pca", group.by = "Sample", shuffle = TRUE, seed = seed) +
  ggtitle("PCA Plot Post-Normalisation")
ggsave(file.path(Out_dir, "03_pca_plot.pdf"), plot = pca_plot, width = 7, height = 6)
ggsave(file.path(Out_dir, "03_pca_plot.png"), plot = pca_plot, width = 7, height = 6, dpi = 150)

umap_plot <- DimPlot(merged_df, reduction = "umap", group.by = "Sample", shuffle = TRUE, seed = seed) +
  ggtitle(sprintf("UMAP on PCs 1-%d", umap_dims)) + labs(x = "UMAP 1", y = "UMAP 2")
ggsave(file.path(Out_dir, "04_umap_plot.pdf"), plot = umap_plot, width = 7, height = 6)
ggsave(file.path(Out_dir, "04_umap_plot.png"), plot = umap_plot, width = 7, height = 6, dpi = 150)

# One panel per sample, in DimPlot's default sample colours
samples <- sort(unique(merged_df$Sample))
sample_cols <- setNames(scales::hue_pal()(length(samples)), samples)

pca_by_sample <- wrap_plots(dimplot_highlight_by_sample(merged_df, "pca", sample_cols), nrow = 1)
ggsave(file.path(Out_dir, "03b_pca_plot_by_sample.pdf"), plot = pca_by_sample, width = 18, height = 6)
ggsave(file.path(Out_dir, "03b_pca_plot_by_sample.png"), plot = pca_by_sample, width = 18, height = 6, dpi = 300)

umap_by_sample <- wrap_plots(dimplot_only_by_sample(merged_df, "umap", sample_cols), nrow = 1)
ggsave(file.path(Out_dir, "04b_umap_plot_by_sample.pdf"), plot = umap_by_sample, width = 18, height = 6)
ggsave(file.path(Out_dir, "04b_umap_plot_by_sample.png"), plot = umap_by_sample, width = 18, height = 6, dpi = 300)


# Log ----
sink(file.path(Out_dir, "Normalisation_log.txt"))
cat("==================================================\n")
cat("Normalisation Settings & Stats\n")
cat("Generated on:", as.character(Sys.time()), "\n")
cat("==================================================\n")
cat("Input:", Input_rds, "\n")
cat("Method: LogNormalize, scale factor", norm_scale_factor, "\n")
cat("Variable features: vst,", norm_nfeatures, "per sample layer, consensus of", length(data_layers), "layers\n")
cat("Variables regressed: none\n")
cat("npcs (PCA):", norm_npcs, "\n")
cat("UMAP: PCs 1-", umap_dims, ", cosine\n", sep = "")
cat("seed:", seed, "\n")
cat("--------------------------------------------------\n")
cat("Total cells:", ncol(merged_df), "\n")
cat("Total genes:", nrow(merged_df), "\n")
print(table(merged_df$Sample))
cat("Variable features (consensus):", length(var_features), "\n")
cat("Top 10 features:", paste(head(var_features, 10), collapse = ", "), "\n")
cat("--------------------------------------------------\n")
cat(sprintf("Variance explained, PCs 1-%d: %.1f%% of total\n", umap_dims, var_df$cumulative_pct_total[umap_dims]))
cat(sprintf("Variance explained, PCs 1-%d: %.1f%% of total\n", norm_npcs, var_df$cumulative_pct_total[norm_npcs]))
cat(sprintf("Share of the variance in PCs 1-%d explained by sample: %.1f%%\n", umap_dims,
            100 * sum(var_df$pct_total_sample[1:umap_dims]) / sum(var_df$pct_total_variance[1:umap_dims])))
cat("\nVariance explained per PC (% of total) and R^2 with Sample, PCs 1-10:\n")
cat(paste(sprintf("  PC%-3d %.2f%%   R^2 %.3f", var_df$PC[1:10], var_df$pct_total_variance[1:10], var_df$sample_r2[1:10]),
          collapse = "\n"), "\n")
cat(sprintf("\nRun time: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
sink()

mqc <- c("id: 'normalisation_stats'",
         "section_name: 'Normalisation Statistics'",
         "description: 'Log normalisation statistics'",
         "plot_type: 'generalstats'",
         "",
         "data:",
         "  Normalisation:",
         "    Method: 'LogNormalize'",
         paste0("    Samples: ", length(data_layers)),
         paste0("    Total_Cells: ", ncol(merged_df)),
         paste0("    Total_Genes: ", nrow(merged_df)),
         paste0("    Variable_Features: ", length(var_features)))
writeLines(mqc, file.path(Out_dir, "normalisation_mqc.yaml"))


# Save ----
log_msg("Saving the normalised object")
saveRDS(merged_df, file.path(Out_dir, "merged_normalised_seurat.rds"))
log_msg("Done (", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min)")
