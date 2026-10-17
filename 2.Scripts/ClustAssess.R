### Script Description ---
# Clustering stability of the blastoid object with ClustAssess
# [https://github.com/Core-Bioinformatics/ClustAssess] - the "Graph clustering" step of the
# package's "Stability-based parameter assessment" vignette (assess_clustering_stability).
#   - input: 3.Results/6.Clustering/Louvain/seurat_multi_res_clustered.rds
#     (93,185 cells) + SingleR (Zhao 2024) labels (3.Results/7.Annotation/singler_zhao2024)
#   - graph as in Clustering.R: SNN graph on 20 PCs, k = 20
#   - resolutions 0.1-0.8; each resolution x algorithm is clustered n_repetitions times
#     with different random seeds, and the stability of the partitions is measured with
#     element-centric consistency (ECC)
#   - algorithms: Louvain (1, what Clustering.R used), Louvain refined (2) and SLM (3)
#   - output: 3.Results/6.Clustering/ClustAssess/
#       clustering_stability.rds, Clustering_stability.pdf
# Usage:  Rscript 2.Scripts/ClustAssess.R          (full run, from any directory)
#         Rscript 2.Scripts/ClustAssess.R test     (quick check on 3,000 cells)
#         N_CORES sets the number of parallel workers (default 1)

args <- commandArgs(trailingOnly = TRUE)
test <- "test" %in% args

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
Out_dir    <- file.path(Root, "3.Results", "6.Clustering", "ClustAssess", if (test) "test" else "")
dir.create(Out_dir, showWarnings = FALSE, recursive = TRUE)

# Load libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ClustAssess)
  library(ggplot2)
  library(patchwork)
})
theme_set(theme_classic())
message("ClustAssess ", packageVersion("ClustAssess"))

# Parameters (graph as in Clustering.R)
n_dims        <- 20                    # PCs used for the graph
k_param       <- 20                    # neighbours for the SNN graph (FindNeighbors default)
resolutions   <- seq(0.1, 0.8, by = 0.1)
clustering_algorithms <- 1:3           # Louvain, Louvain refined, SLM
n_repetitions <- 10                    # random seeds per resolution x algorithm (vignette: 10)
ncores        <- as.integer(Sys.getenv("N_CORES", "1"))
set.seed(42)
if (test) { n_repetitions <- 2; resolutions <- c(0.3, 0.6) }

# 1. Load the object and add the SingleR annotation
seu <- readRDS(input_rds)
sr  <- readRDS(singler_md)             # data.frame: one row per cell, rownames = cell barcodes
message("Object: ", ncol(seu), " cells | SingleR table: ", nrow(sr), " cells x ", ncol(sr), " columns")
missing <- setdiff(colnames(seu), rownames(sr))
if (length(missing) > 0) {
  stop(length(missing), " cells in the object have no SingleR annotation, e.g. ",
       paste(head(missing, 3), collapse = ", "))
}
seu <- AddMetaData(seu, metadata = sr[colnames(seu), , drop = FALSE])
print(table(seu$singler_label, useNA = "ifany"))

if (test) {
  seu <- subset(seu, cells = sample(colnames(seu), 3000))
  message("TEST MODE: ", ncol(seu), " cells, ", n_repetitions, " repetitions")
}

# 2. SNN graph on the first 20 PCs (as in the clustering pipeline)
pca_emb    <- Embeddings(seu, "pca")[, seq_len(n_dims)]
adj_matrix <- FindNeighbors(pca_emb, k.param = k_param, verbose = FALSE)$snn

# 3. Parallel backend
RhpcBLASctl::blas_set_num_threads(1)
if (ncores > 1) {
  my_cluster <- parallel::makeCluster(ncores, type = "PSOCK")
  parallel::clusterCall(my_cluster, function(lib) .libPaths(lib), .libPaths())
  doParallel::registerDoParallel(cl = my_cluster)
  message("Parallel backend: ", ncores, " workers")
}

# 4. Clustering stability: resolution x algorithm, n_repetitions seeds each
t0 <- Sys.time()
clustering_stability <- assess_clustering_stability(
  graph_adjacency_matrix = adj_matrix,
  resolution = resolutions,
  n_repetitions = n_repetitions,
  clustering_algorithm = clustering_algorithms,
  ecs_thresh = 1
)
message("Clustering stability done in ",
        sprintf("%.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
saveRDS(clustering_stability, file.path(Out_dir, "clustering_stability.rds"))

# 5. Plots (as in the vignette)
legend_br <- theme(legend.position = "inside", legend.position.inside = c(1, 0), legend.justification = c(1, 0))
pdf(file.path(Out_dir, "Clustering_stability.pdf"), width = 16, height = 8)
print(wrap_plots(list(
  plot_clustering_per_value_stability(clustering_stability, "k") + legend_br,
  plot_clustering_per_value_stability(clustering_stability, "resolution") + legend_br), nrow = 1))
print(wrap_plots(list(
  plot_clustering_overall_stability(clustering_stability, "k") + legend_br,
  plot_clustering_overall_stability(clustering_stability, "resolution", mean) + legend_br), nrow = 1))
print(plot_k_resolution_corresp(clustering_stability, dodge_width = 0.7, summary_function = min) +
        plot_annotation(title = "resolution - k correspondence (colour: min ECC)"))
print(plot_k_resolution_corresp(clustering_stability, dodge_width = 0.7, colour_information = "freq_k") +
        plot_annotation(title = "resolution - k correspondence (colour: frequency of k)"))
print(plot_k_n_partitions(clustering_stability, dodge_width = 0.5) +
        plot_annotation(title = "k - # partitions correspondence (colour: ECC)"))
print(plot_k_n_partitions(clustering_stability, dodge_width = 0.5, colour_information = "freq_part") +
        plot_annotation(title = "k - # partitions correspondence (colour: frequency of the most common partition)"))
dev.off()

# Close the parallel workers
if (exists("my_cluster")) parallel::stopCluster(my_cluster)
foreach::registerDoSEQ()
message("Outputs in ", Out_dir)
