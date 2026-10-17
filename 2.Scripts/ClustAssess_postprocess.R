### Script Description ---
# Follow-up to ClustAssess.R: WHERE is the clustering unstable?
# Uses the saved ClustAssess results (no reclustering), Louvain only (the algorithm
# behind the saved RNA_snn_res.* clusters), resolutions 0.1-0.8.
#   - per-cell ECC on the UMAP
#   - ECC per saved cluster (RNA_snn_res.*) and per SingleR label
#   - mean ECS between the saved clusters and ClustAssess's 10 Louvain runs
# ECC = how consistently a cell keeps the same cluster-mates across random seeds (1 = always).
# Output: 3.Results/6.Clustering/ClustAssess/
#   ECC_per_cluster.csv, ECC_per_singler_label.csv, ECS_saved_clusters_vs_runs.csv, ECC_followup.pdf
# Usage:  Rscript 2.Scripts/ClustAssess_postprocess.R   (from any directory)

# Path configuration
# Repository root: the folder above 2.Scripts (from the script's own path with Rscript or source(),
# otherwise the working directory)
Root <- local({
  f <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), fixed = TRUE)
  if (length(f) == 0) f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (length(f) == 1) normalizePath(file.path(dirname(f), "..")) else getwd()
})
if (!dir.exists(file.path(Root, "2.Scripts"))) stop("Repository root not found (no 2.Scripts folder in ", Root, ")")
input_rds  <- file.path(Root, "3.Results", "6.Clustering", "Louvain", "seurat_multi_res_clustered.rds")
singler_md <- file.path(Root, "3.Results", "7.Annotation", "singler_zhao2024", "singler_metadata.rds")
Out_dir    <- file.path(Root, "3.Results", "6.Clustering", "ClustAssess")

library(Seurat)
library(ClustAssess)
library(ggplot2)
library(patchwork)

# 1. Load the object, SingleR labels and the ClustAssess results
seu <- readRDS(input_rds)
sr  <- readRDS(singler_md)
seu$singler_label <- droplevels(sr[colnames(seu), "singler_label"])
cs  <- readRDS(file.path(Out_dir, "clustering_stability.rds"))
louvain <- cs$split_by_resolution$Louvain              # one entry per resolution ("0.1" ... "0.8")
resolutions <- names(louvain)

# 2. Add per-cell ECC (cells are in the same order as the object)
for (res in resolutions) {
  seu[[paste0("ecc_res", res)]] <- louvain[[res]]$ecc
}
ecc_cols <- paste0("ecc_res", resolutions)

# 3. Tables: median ECC per saved cluster and per SingleR label; ECS of saved clusters vs the 10 runs
per_cluster <- do.call(rbind, lapply(resolutions, function(res) {
  cl <- seu[[paste0("RNA_snn_res.", res)]][, 1]
  ecc <- seu[[paste0("ecc_res", res)]][, 1]
  data.frame(resolution = res, cluster = levels(factor(cl)),
             n_cells = as.vector(table(cl)), median_ecc = round(tapply(ecc, cl, median), 3))
}))
write.csv(per_cluster, file.path(Out_dir, "ECC_per_cluster.csv"), row.names = FALSE)

per_label <- sapply(ecc_cols, function(e) tapply(seu[[e]][, 1], seu$singler_label, median))
colnames(per_label) <- resolutions
write.csv(round(per_label, 3), file.path(Out_dir, "ECC_per_singler_label.csv"))

saved_vs_runs <- sapply(resolutions, function(res) {
  saved <- as.integer(factor(seu[[paste0("RNA_snn_res.", res)]][, 1]))   # labels 1..k (ClustAssess expects 1-based)
  runs  <- unlist(lapply(louvain[[res]]$clusters, function(k) lapply(k$partitions, function(p) p$mb)), recursive = FALSE)
  # mean of the per-cell ECS (element_sim() itself returns invalid values at ~93k cells)
  mean(sapply(runs, function(mb) mean(element_sim_elscore(saved, mb))))
})
write.csv(data.frame(resolution = resolutions, mean_ECS_saved_vs_runs = round(saved_vs_runs, 3)),
          file.path(Out_dir, "ECS_saved_clusters_vs_runs.csv"), row.names = FALSE)
print(round(saved_vs_runs, 3))

# 4. Plots
pdf(file.path(Out_dir, "ECC_followup.pdf"), width = 16, height = 9)

# ECC on the UMAP, all resolutions
print(wrap_plots(lapply(ecc_cols, function(e)
  FeaturePlot(seu, features = e, reduction = "umap", raster = TRUE) +
    scale_colour_viridis_c(limits = c(0, 1)) + ggtitle(sub("ecc_res", "ECC, resolution ", e))), ncol = 4))

# Per resolution: saved clusters on the UMAP + ECC per cluster
for (res in resolutions) {
  cl_col <- paste0("RNA_snn_res.", res)
  d <- data.frame(cluster = seu[[cl_col]][, 1], ecc = seu[[paste0("ecc_res", res)]][, 1])
  p1 <- DimPlot(seu, group.by = cl_col, reduction = "umap", label = TRUE, raster = TRUE) +
    NoLegend() + ggtitle(paste("Saved clusters,", cl_col))
  p2 <- ggplot(d, aes(cluster, ecc, fill = cluster)) +
    geom_boxplot(outlier.size = 0.2) + coord_cartesian(ylim = c(0, 1)) + NoLegend() +
    labs(title = paste("ECC per cluster, resolution", res), x = "Cluster", y = "ECC")
  print(p1 | p2)
}

# Median ECC per SingleR label x resolution
hm <- as.data.frame(as.table(per_label)); colnames(hm) <- c("label", "resolution", "median_ecc")
print(ggplot(hm, aes(resolution, label, fill = median_ecc)) +
        geom_tile() + geom_text(aes(label = sprintf("%.2f", median_ecc)), size = 3) +
        scale_fill_viridis_c(limits = c(0, 1)) +
        labs(title = "Median ECC per SingleR label (Louvain)", x = "Resolution", y = NULL))
dev.off()
message("Outputs in ", Out_dir)
