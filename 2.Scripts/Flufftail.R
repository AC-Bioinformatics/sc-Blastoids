# Flufftail [https://github.com/Core-Bioinformatics/Flufftail] at Louvain resolution 0.3
#
# Input:  3.Results/6.Clustering/Louvain/seurat_multi_res_clustered.rds (RNA_snn graph, UMAP)
#         3.Results/7.Annotation/singler_zhao2024/singler_metadata.rds (per-cell SingleR lineage)
# Steps:  1. Fuzzy clustering at resolution 0.3 (30 Louvain seeds); the configuration found by most seeds is used
#         2. Pseudotime from naive epiblast to CTB/STB (primitive endoderm removed)
#         3. DE genes at the transition between the two 24h epiblast clusters 2 and 4 (Louvain_0.3_11),
#            core markers, cell-cycle phase per state and KEGG / GO:BP enrichment (g:Profiler)
# Output: 3.Results/8.Pseudotime/flufftail/res0.3/
#           fuzzy_clust_res0.3.rds, membership_degrees_<config>.csv, ecc_boxplots.pdf, umap_ecc.pdf,
#           umap_hard_clusters.pdf, consensus_heatmap.pdf, umap_entropy_membership.pdf
#           cds_pseudotime_end_CTB_STB.rds, umap_pseudotime_end_CTB_STB.pdf, pseudotime_decomposed_end_CTB_STB.pdf
#           transition_2_vs_4/ (transition cells, DE genes, heatmap, core markers, cell-cycle phase, enrichment)
#
# Usage:  Rscript 2.Scripts/Flufftail.R   (from any directory)

suppressPackageStartupMessages({
  library(Flufftail)
  library(Seurat)
  library(SeuratWrappers)
  library(monocle3)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(gprofiler2)
})

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
Out_dir    <- file.path(Root, "3.Results", "8.Pseudotime", "flufftail", "res0.3")
dir.create(Out_dir, showWarnings = FALSE, recursive = TRUE)

# Loading data
seurat_input <- readRDS(input_rds)


# 1. Fuzzy clustering ----
set.seed(42)
seurat_object <- seurat_input
sample_names <- colnames(seurat_object)
umap <- Embeddings(seurat_object, "umap")

# SNN graph from Seurat object (this object has RNA_snn)
adj_matrix <- seurat_object@graphs$RNA_snn

# Run fuzzy clustering analysis at resolution 0.3 only
fuzzy_clust_analysis <- fuzzy_clustering_analysis(
    adj_matrix,
    sample_names = sample_names,
    n_repetitions = 30,
    algorithm = 1,
    resolution = 0.3
)

# Resolution 0.3 gives several configurations (Louvain_0.3_<k>, k = number of clusters).
# Use the k found by the most seeds; a configuration from a single seed has ECC = 1 by definition.
n_seeds <- sapply(fuzzy_clust_analysis$membership_lists, length)
print(n_seeds)
config <- names(which.max(n_seeds))
cat("Using", config, "\n")

theme_large_text <- theme(
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 15),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 15)
  )

# ECC of all configurations, sorted by median
p_ecc_box <- plot_ecc_boxplots(fuzzy_clust_analysis$ecc_list, length(n_seeds), select_highest = TRUE)
ggsave(file.path(Out_dir, "ecc_boxplots.pdf"), p_ecc_box, width = 8, height = 6)

# UMAP coloured by ECC and hard clusters
ecc <- fuzzy_clust_analysis$ecc_list[[config]]
hard_clusters <- fuzzy_clust_analysis$hard_clusters_list[[config]]

p_ecc <- plot_umap(umap, ecc, "ECC", paste(config, "- ECC")) + theme_large_text
ggsave(file.path(Out_dir, "umap_ecc.pdf"), p_ecc, width = 9, height = 7)

p_hard <- plot_umap(umap, hard_clusters, "Hard clusters", paste(config, "- Hard Clusters")) + theme_large_text
ggsave(file.path(Out_dir, "umap_hard_clusters.pdf"), p_hard, width = 9, height = 7)

# Membership degree matrix (first column is the cell name, then one column per cluster)
degree_matrix <- fuzzy_clust_analysis$degree_matrices_list[[config]]
head(degree_matrix)
write.csv(degree_matrix, file.path(Out_dir, paste0("membership_degrees_", config, ".csv")), row.names = FALSE)

# Consensus matrix on 5000 random cells (the full 93k x 93k matrix would need ~70 GB)
cells_5000 <- sample(length(sample_names), 5000)
membership_5000 <- lapply(fuzzy_clust_analysis$membership_lists[[config]], function(x) x[cells_5000])
consensus <- get_consensus_matrix(membership_5000, sample_names[cells_5000])
consensus_plot <- subsample_consensus_mat_plot(consensus, n_sample = 5000, title = paste(config, "- Consensus"))
ggsave(file.path(Out_dir, "consensus_heatmap.pdf"), consensus_plot$plot$gtable, width = 8, height = 8)

# Entropy on the membership matrix (0 = cell always in one cluster, higher = shared between clusters)
entropy_membership <- calculate_entropy(degree_matrix[, -1])
p_entropy <- plot_umap(umap, entropy_membership, "Entropy", paste(config, "- Entropy on the Membership Matrix")) + theme_large_text
ggsave(file.path(Out_dir, "umap_entropy_membership.pdf"), p_entropy, width = 9, height = 7)

# Save results for the pseudotime and transition steps
fz <- list(fuzzy_clust_analysis = fuzzy_clust_analysis, config = config)
saveRDS(fz, file.path(Out_dir, "fuzzy_clust_res0.3.rds"))


# 2. Pseudotime ----
set.seed(42)
seurat_object <- seurat_input
umap <- Embeddings(seurat_object, "umap")
hard_clusters <- fuzzy_clust_analysis$hard_clusters_list[[config]]
names(hard_clusters) <- colnames(seurat_object)

# Cell type labels for the plots: per-cell SingleR lineage (Zhao 2024 reference)
singler <- readRDS(singler_md)
seurat_object$cell_type <- as.character(singler[colnames(seurat_object), "singler_lineage"])

# Drop primitive endoderm: the PrE island (UMAP_2 > 7) and any other SingleR PrE cells; pseudotime is epiblast -> TE
keep <- colnames(seurat_object)[umap[, 2] <= 7 & seurat_object$cell_type != "Primitive endoderm"]
seurat_object <- subset(seurat_object, cells = keep)
umap <- umap[keep, ]
hard_clusters <- hard_clusters[keep]
cat("Cells after dropping PrE:", ncol(seurat_object), "\n")

# The pseudotime function reads seurat_object@assays$RNA@counts (Seurat v3 layout),
# so the v5 assay is joined and converted (scale.data is removed first, it is not needed here)
for (l in grep("^scale.data", Layers(seurat_object), value = TRUE)) LayerData(seurat_object, layer = l) <- NULL
seurat_object[["RNA"]] <- JoinLayers(seurat_object[["RNA"]])
seurat_object[["RNA"]] <- as(seurat_object[["RNA"]], "Assay")

# Pseudotime: start = naive epiblast markers, end = CTB/STB markers.
# Pan-TE markers (GATA3, GATA2, KRT7, TACSTD2) put the end in the mature TE body, and KRT7 is ~0 here.
# End genes: higher % and level in both CTB and STB than in Mature TE;
# an end cell must be high for 5 of the 7 (end_relax_ngenes = 2)
start_genes <- c("NANOG", "KLF17", "DPPA3", "DPPA5")
end_genes <- c("CAPN6", "TNFAIP3", "VGLL1", "CYP19A1", "GCM1", "SDC1", "NR2F2")
# Tag added to output file names, so runs with different end genes do not overwrite each other
run_tag <- "end_CTB_STB"

cds <- get_pseudotime_trajectory_inference(
  seurat_object,
  colnames(seurat_object),
  union(VariableFeatures(seurat_object), c(start_genes, end_genes)),
  umap,
  start_genes = start_genes,
  start_relax_ngenes = 1,
  end_genes = end_genes,
  end_relax_ngenes = 2,
  threshold = 0.8
)

# UMAP coloured by pseudotime
p_pseudotime <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,
                           label_branch_points = FALSE, label_roots = FALSE, cell_size = 0.5)
ggsave(file.path(Out_dir, paste0("umap_pseudotime_", run_tag, ".pdf")), p_pseudotime, width = 9, height = 7)

# Cell types and hard clusters along pseudotime
decomposed_plot <- get_decomposed_pseudotime_plot(cds, as.integer(as.character(hard_clusters[colnames(cds)])))
ggsave(file.path(Out_dir, paste0("pseudotime_decomposed_", run_tag, ".pdf")), decomposed_plot$combined, width = 12, height = 10)

saveRDS(cds, file.path(Out_dir, paste0("cds_pseudotime_", run_tag, ".rds")))
rm(cds, seurat_object, singler); invisible(gc())


# 3. Transition between clusters 2 and 4 ----
# Flufftail DE genes at the transition between the two 24h epiblast clusters, 2 and 4 (Louvain_0.3_11).
# Both are Early EPI by SingleR and ~90%+ 24h, but separate on the UMAP; their border cells have high entropy.
transition_config <- "Louvain_0.3_11"
clusters <- c(2, 4)   # state A = cluster 2 (pseudotime start), state B = cluster 4
if (config != transition_config) {
  message("Fuzzy clustering gave ", config, ", not ", transition_config, "; clusters ", clusters[1], " and ",
          clusters[2], " may not be the two 24h epiblast clusters, so the transition step is skipped")
  quit(save = "no", status = 0)
}

set.seed(42)
pair <- paste(clusters, collapse = "_")
Tr_dir <- file.path(Out_dir, paste0("transition_", clusters[1], "_vs_", clusters[2]))
dir.create(Tr_dir, showWarnings = FALSE, recursive = TRUE)

seurat_object <- seurat_input
rm(seurat_input); invisible(gc())
hard_clusters <- as.integer(as.character(fz$fuzzy_clust_analysis$hard_clusters_list[[config]]))
degree_matrix <- as.matrix(fz$fuzzy_clust_analysis$degree_matrices_list[[config]][, -1])
entropy <- calculate_entropy(degree_matrix)
names(hard_clusters) <- names(entropy) <- rownames(degree_matrix) <- colnames(seurat_object)

# All cells are kept (as in the vignette); get_transition_cells only uses the two clusters of interest
umap <- Embeddings(seurat_object, "umap")
seurat_object$hard_cluster <- hard_clusters
seurat_object$cell_type <- paste("Cluster", hard_clusters)   # used by the heatmap annotation
cat("Cells:", ncol(seurat_object), "\n"); print(table(seurat_object$hard_cluster))

# Seurat v3 layout for SeuratWrappers / monocle3 (as in the pseudotime step)
for (l in grep("^scale.data", Layers(seurat_object), value = TRUE)) LayerData(seurat_object, layer = l) <- NULL
seurat_object[["RNA"]] <- JoinLayers(seurat_object[["RNA"]])
seurat_object[["RNA"]] <- as(seurat_object[["RNA"]], "Assay")

# Transition cells. The midpoint between the two centroids is moved to the median position of the
# border cells (membership >= 0.1 in both clusters), where entropy is highest.
umap_df <- as.data.frame(umap)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$hard_clusters <- hard_clusters

border <- colnames(seurat_object)[degree_matrix[, as.character(clusters[1])] >= 0.1 & degree_matrix[, as.character(clusters[2])] >= 0.1]
border_median <- apply(umap[border, ], 2, median)
cat("Border cells:", length(border), "| median UMAP:", round(border_median, 2), "\n")

first_pass <- get_transition_cells(umap_df, clusters_of_interest = clusters, percentage_points_centroid = 0.4, num_points_transition = 1000)
midpoint <- unlist(first_pass$centroid_results$midpoint_centroid[1, 1:2])
translation <- as.numeric(border_median - midpoint)
cat("Midpoint:", round(midpoint, 2), "| translation:", round(translation, 2), "\n")

transition_cells_results <- get_transition_cells(umap_df, clusters_of_interest = clusters, percentage_points_centroid = 0.4,
                                                 x_y_translation = translation, num_points_transition = 1000)
cells_A <- transition_cells_results$centroids_cells[[clusters[1]]]
cells_B <- transition_cells_results$centroids_cells[[clusters[2]]]
cells_T <- transition_cells_results$cells_at_transition

# Check that the transition cells are the high-entropy ones
cat("Median entropy | cluster", clusters[1], "core:", round(median(entropy[cells_A]), 3),
    "| transition:", round(median(entropy[cells_T]), 3),
    "| cluster", clusters[2], "core:", round(median(entropy[cells_B]), 3), "\n")

umap_df$state <- "other"
umap_df[cells_A, "state"] <- paste("Cluster", clusters[1], "core")
umap_df[cells_B, "state"] <- paste("Cluster", clusters[2], "core")
umap_df[cells_T, "state"] <- "Transition"
umap_df$entropy <- entropy

p_states <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = state)) + geom_point(size = 0.4) +
  scale_color_manual(values = setNames(c("gray85", "#F8766D", "#619CFF", "purple"),
                                       c("other", paste("Cluster", clusters, "core"), "Transition"))) +
  theme_classic() + labs(color = "Cell state", title = "Cells at the centroids and at the transition")
ggsave(file.path(Tr_dir, "transition_cells.pdf"), p_states, width = 9, height = 7)

p_entropy <- plot_umap(umap, umap_df$entropy, "Entropy", "Membership entropy (Louvain_0.3_11)")
ggsave(file.path(Tr_dir, paste0("entropy_clusters_", pair, ".pdf")), p_entropy, width = 9, height = 7)

# Trajectory learned on all cells (monocle3), rooted at the cluster A core cells so pseudotime runs from A to B
cds <- as.cell_data_set(seurat_object)
cds <- cds[VariableFeatures(seurat_object), ]
fData(cds)$gene_short_name <- rownames(fData(cds))
cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition = FALSE)
cds <- order_cells(cds, root_cells = cells_A)
cds <- estimate_size_factors(cds)

p_pseudotime <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,
                           label_branch_points = FALSE, label_roots = FALSE, cell_size = 0.5)
ggsave(file.path(Tr_dir, paste0("pseudotime_", clusters[1], "_to_", clusters[2], ".pdf")), p_pseudotime, width = 9, height = 7)

# DE genes at the transition (Moran's I along the trajectory, transition cells only)
DE_genes_results <- identify_DE_genes_at_the_transition(cds, cells_T, q_value_thresh = 0.05, morans_I_thresh = 0.02)
de_genes <- DE_genes_results$test_results %>% mutate(morans_I = as.numeric(morans_I)) %>% arrange(desc(morans_I))
cat("DE genes at the transition:", nrow(de_genes), "\n"); print(head(de_genes[, c("gene_short_name", "morans_I", "q_value")], 30))
write.csv(de_genes, file.path(Tr_dir, paste0("DE_genes_transition_", pair, ".csv")), row.names = FALSE)

# The heatmap matches cells by name after data.frame() has run make.names() on them; our names start with
# a digit ("24h-..."), which make.names prefixes with "X", so "X" is added here to make the names match
# (also to the stored pseudotime, which monocle3 looks up by cell name)
cds_heatmap <- DE_genes_results$cds_sub
pt_aux <- principal_graph_aux(cds_heatmap)$UMAP$pseudotime
names(pt_aux) <- paste0("X", names(pt_aux))
colnames(cds_heatmap) <- paste0("X", colnames(cds_heatmap))
principal_graph_aux(cds_heatmap)$UMAP$pseudotime <- pt_aux
meta_heatmap <- seurat_object@meta.data
rownames(meta_heatmap) <- paste0("X", rownames(meta_heatmap))
p_heatmap <- plot_DE_genes_heatmap(cds_heatmap, meta_heatmap, de_genes$gene_short_name[1:min(20, nrow(de_genes))], 10)
ggsave(file.path(Tr_dir, "DE_genes_transition_heatmap.pdf"), p_heatmap, width = 10, height = 8)

# Direct comparison of the two cluster cores (which side of the border each gene is higher on)
Idents(seurat_object) <- "hard_cluster"
seurat_object <- NormalizeData(seurat_object, verbose = FALSE)
core_markers <- FindMarkers(seurat_object, ident.1 = cells_A, ident.2 = cells_B, logfc.threshold = 0.25, min.pct = 0.1)
core_markers$gene <- rownames(core_markers)
core_markers$higher_in <- ifelse(core_markers$avg_log2FC > 0, paste("Cluster", clusters[1]), paste("Cluster", clusters[2]))
core_markers <- core_markers %>% filter(p_val_adj < 0.05) %>% arrange(desc(abs(avg_log2FC)))

# Cell-cycle phase per state (the top transition genes are S and G2/M genes)
seurat_object <- CellCycleScoring(seurat_object, s.features = cc.genes.updated.2019$s.genes,
                                  g2m.features = cc.genes.updated.2019$g2m.genes)
phase_table <- table(umap_df[colnames(seurat_object), "state"], seurat_object$Phase)
cat("\nCell-cycle phase per state (% of cells):\n"); print(round(100 * prop.table(phase_table, 1), 1))
write.csv(as.data.frame.matrix(phase_table), file.path(Tr_dir, "cell_cycle_phase_per_state.csv"))
p_phase <- DimPlot(seurat_object, group.by = "Phase", raster = TRUE) + ggtitle("Cell-cycle phase")
ggsave(file.path(Tr_dir, "umap_cell_cycle_phase.pdf"), p_phase, width = 9, height = 7)
write.csv(core_markers, file.path(Tr_dir, paste0("markers_core_", clusters[1], "_vs_", clusters[2], ".csv")), row.names = FALSE)
cat("\nTop markers, cluster", clusters[1], "core vs cluster", clusters[2], "core:\n")
print(head(core_markers[core_markers$avg_log2FC > 0, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")], 20))
print(head(core_markers[core_markers$avg_log2FC < 0, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")], 20))

# Transition DE genes with the side they are higher on
de_genes <- left_join(de_genes, core_markers[, c("gene", "avg_log2FC", "pct.1", "pct.2", "higher_in")], by = c("gene_short_name" = "gene"))
write.csv(de_genes, file.path(Tr_dir, paste0("DE_genes_transition_", pair, ".csv")), row.names = FALSE)

saveRDS(list(cds = cds, transition_cells_results = transition_cells_results, de_genes = de_genes, core_markers = core_markers),
        file.path(Tr_dir, paste0("transition_", pair, "_objects.rds")))

# KEGG and GO:BP enrichment of the transition DE genes (g:Profiler). The vignette queries all tested genes
# (DE_genes_results$full_results); the DE genes are used here. Background = the genes that
# were tested (the variable features in the trajectory), since the DE genes could only come from those.
tested_genes <- DE_genes_results$full_results$gene_short_name

results <- gost(
  query = de_genes$gene_short_name,
  sources = c("KEGG", "GO:BP"),
  organism = "hsapiens",
  custom_bg = tested_genes,
  as_short_link = FALSE,
  evcodes = TRUE
)

if (is.null(results)) {
  cat("\nNo significant KEGG or GO:BP terms\n")
} else {
  gsea_results <- results$result %>% arrange(desc(intersection_size))
  write.csv(gsea_results[, c("source", "term_id", "term_name", "p_value", "term_size", "query_size", "intersection_size", "intersection")],
            file.path(Tr_dir, paste0("enrichment_DE_genes_transition_", pair, ".csv")), row.names = FALSE)

  # One plot per source; top 5 by intersection size highlighted
  for (src in unique(gsea_results$source)) {
    src_results <- gsea_results %>% filter(source == src)
    cat("\n", src, ":", nrow(src_results), "significant terms\n")
    print(head(src_results[, c("term_name", "p_value", "term_size", "query_size", "intersection_size")], 10))
    top5 <- src_results %>% dplyr::slice(1:5)

    p_enrich <- ggplot(src_results, aes(x = reorder(term_name, -log10(p_value)),
                                        y = -log10(p_value),
                                        size = intersection_size)) +
      geom_point(color = "grey70", alpha = 0.6) +
      geom_point(data = top5, aes(color = term_name), size = 4) +
      geom_text_repel(data = top5,
                      aes(label = term_name),
                      size = 5,
                      max.overlaps = Inf,
                      box.padding = 0.6,
                      segment.color = "grey60") +
      theme_classic() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)
      ) +
      expand_limits(
        y = max(-log10(src_results$p_value)) + 1,
        x = c(0.5, nrow(src_results) + 5)
      ) +
      labs(x = paste(src, "terms"), y = "-log10(p-adj)") +
      guides(color = "none")
    ggsave(file.path(Tr_dir, paste0(gsub(":", "", src), "_DE_genes_transition_", pair, ".pdf")), p_enrich, width = 10, height = 7)
  }
}
