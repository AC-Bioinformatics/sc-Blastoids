### SingleR annotation of the APAP blastoid data against the Zhao et al. 2024
### human embryo reference (Nat Methods, doi 10.1038/s41592-024-02493-2), per cell
### and per cluster at every clustering resolution.
###
### Reference: 1.Data/references/zhao2024_human_embryo/psd.R3.6.em.seurat.ob.rds,
### 3,304 Smart-seq2 cells from six integrated embryo datasets (zygote to CS7),
### log-normalised data and the authors' sub-lineage labels (sub_rename_EML).
###   - CS7 gastrula cells (Tyser 2021) removed, except amnion (keep_amnion): the
###     embryo datasets have no other amnion annotation, and amnion-like cells are
###     a known off-target population in blastoids.
###   - Labels regrouped (ref_map): Zygote / 2-4 cell / 8 cell -> Cleavage;
###     Early/Late_EPI -> Early/Late EPI; Early/Late_Hypoblast -> Early/Late PrE;
###     TE at E5-E6 -> Early TE, TE at E7-E8 -> Mature TE (in the reference NR2F2,
###     CGA, CCR7 and VGLL1 rise and CDX2 falls from E6 to E7); Ambiguous and
###     Unknown dropped. Labels with fewer than min_ref_cells cells are dropped.
###     Each label also has a coarse lineage (lineage_map).
###   The prepared reference is saved next to the original for reuse.
### SingleR (trainSingleR / classifySingleR; markers from pairwise Wilcoxon tests
### between reference labels; Spearman correlation on those markers, 0.8 quantile
### per label, fine-tuning on the closest labels):
###   - per cell: label, pruned label (NA = low confidence: delta to the next label
###     more than 3 MADs below that label's median), scores, delta;
###   - per cluster (every resolution): SingleR on the mean log-normalised profile
###     of the cluster, the runner-up label, and the majority per-cell label.
###
### Input: 3.Results/6.Clustering/Louvain/seurat_multi_res_clustered.rds
###
### Usage: Rscript 2.Scripts/Annotation_SingleR.R [--test]   (from any directory)
###   --test  5,000 random cells and two resolutions, written to <out_dir>/test
###
### Outputs (3.Results/7.Annotation/singler_zhao2024/):
###   singler_annotation.pdf        reference labels by stage, per-cell labels and
###                                 lineages, scores and deltas, labels by sample, cluster
###                                 labels at every resolution, one page per resolution
###                                 (clusters, cluster label, per-cell labels per cluster,
###                                 sample composition), labels by sample
###   singler_reference_labels.csv  reference label x stage x dataset (cell counts)
###   singler_label_umaps.pdf       one UMAP per per-cell label (its cells in colour over all
###                                 cells in grey), with the number of cells
###   singler_cell_labels_by_sample.csv  per-cell labels and lineages per sample
###   singler_cluster_labels.csv    per resolution / cluster: SingleR label, runner-up,
###                                 delta, majority cell label and share, % pruned,
###                                 sample composition
###   singler_summary.txt           all of the above as text
###   singler_metadata.rds          per-cell labels, scores, delta and the cluster label
###                                 at each resolution (add with AddMetaData)
library(Seurat)
library(SingleR)
library(SingleCellExperiment)
library(BiocParallel)
library(ggplot2)
library(patchwork)

# ============================================================
# 1. Path configuration
# ============================================================
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

input_rds <- file.path(Results, "6.Clustering", "Louvain", "seurat_multi_res_clustered.rds")
ref_rds   <- file.path(Root, "1.Data", "references", "zhao2024_human_embryo", "psd.R3.6.em.seurat.ob.rds")
out_dir   <- file.path(Results, "7.Annotation", "singler_zhao2024")

# ============================================================
# 2. Parameters
# ============================================================
res_use       <- NULL     # resolutions to annotate, e.g. c(0.3, 0.8); NULL = all in the object
res_prefix    <- "RNA_snn_res."
umap_name     <- "umap"
sample_col    <- "Sample"
keep_amnion   <- TRUE     # keep the CS7 amnion cells when the gastrula is removed
min_ref_cells <- 10       # drop reference labels with fewer cells
de_method     <- "wilcox" # recommended for single-cell references
de_n          <- 20       # markers per pairwise label comparison
n_workers     <- as.integer(Sys.getenv("N_CORES", "4"))
test_cells    <- 0        # > 0: run on a random subset of this many cells
title_theme   <- theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
pt_size       <- 3            # UMAP point size, in pixels of the rasterised layer
raster_dpi    <- c(800, 800)   # resolution of the rasterised UMAP layer

ref_map <- c("Zygote" = "Cleavage", "2-4 cell" = "Cleavage", "8 cell" = "Cleavage",
             Morula = "Morula", Prelineage = "Prelineage", ICM = "ICM",
             "EPI.PrE.INT" = "EPI/PrE intermediate",
             Early_EPI = "Early EPI", Late_EPI = "Late EPI",
             Early_Hypoblast = "Early PrE", Late_Hypoblast = "Late PrE",
             TE = "TE", CTB = "CTB", STB = "STB", EVT = "EVT", ExE_Mes = "ExE mesoderm",
             Amnion.Ecto = "Amnion", Late_Amnion = "Amnion")
lineage_map <- c(Cleavage = "Pre-lineage", Morula = "Pre-lineage", Prelineage = "Pre-lineage",
                 ICM = "ICM", "EPI/PrE intermediate" = "ICM",
                 "Early EPI" = "Epiblast", "Late EPI" = "Epiblast",
                 "Early PrE" = "Primitive endoderm", "Late PrE" = "Primitive endoderm",
                 "Early TE" = "Trophectoderm", "Mature TE" = "Trophectoderm",
                 CTB = "Post-implantation trophoblast", STB = "Post-implantation trophoblast",
                 EVT = "Post-implantation trophoblast",
                 "ExE mesoderm" = "Extraembryonic mesoderm", Amnion = "Amnion")
label_cols <- c(Cleavage = "#D9C2A7", Morula = "#B08968", Prelineage = "#7F5539", ICM = "#F4A261",
                "EPI/PrE intermediate" = "#E9C46A", "Early EPI" = "#D62728", "Late EPI" = "#8B0000",
                "Early PrE" = "#2CA02C", "Late PrE" = "#145A32", "Early TE" = "#6BAED6",
                "Mature TE" = "#1F4E9A", CTB = "#17BECF", STB = "#9467BD", EVT = "#E377C2",
                "ExE mesoderm" = "#BCBD22", Amnion = "#FF7F0E", Pruned = "grey80")
lineage_cols <- c("Pre-lineage" = "#B08968", ICM = "#F4A261", Epiblast = "#D62728",
                  "Primitive endoderm" = "#2CA02C", Trophectoderm = "#1F77B4",
                  "Post-implantation trophoblast" = "#9467BD", "Extraembryonic mesoderm" = "#BCBD22",
                  Amnion = "#FF7F0E", Pruned = "grey80")

args <- commandArgs(trailingOnly = TRUE)
if ("--test" %in% args) {
  test_cells <- 5000
  res_use    <- c(0.2, 0.8)
  out_dir    <- file.path(out_dir, "test")
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
bp <- if (n_workers > 1 && .Platform$OS.type != "windows") MulticoreParam(n_workers) else SerialParam()

# ============================================================
# 3. Reference
# ============================================================
message("Loading reference: ", ref_rds)
ref_seu <- UpdateSeuratObject(readRDS(ref_rds))
rmd <- ref_seu[[]]
ref_lab <- unname(ref_map[rmd$sub_rename_EML])
te <- which(rmd$sub_rename_EML == "TE")
ref_lab[te] <- ifelse(rmd$devTime[te] %in% c("E5", "E6"), "Early TE", "Mature TE")
keep <- !is.na(ref_lab) & (rmd$devTime != "CS7" | (keep_amnion & ref_lab %in% "Amnion"))
n_lab <- table(ref_lab[keep])
small <- names(n_lab)[n_lab < min_ref_cells]
keep  <- keep & !ref_lab %in% small
ref_levels <- intersect(names(lineage_map), unique(ref_lab[keep]))
stage_levels <- intersect(c("Zygote", paste0("E", 1:14), "CS7"), unique(rmd$devTime[keep]))
ref <- SingleCellExperiment(
  assays  = list(logcounts = LayerData(ref_seu, layer = "data")[, keep]),
  colData = DataFrame(label = ref_lab[keep], lineage = unname(lineage_map[ref_lab[keep]]),
                      stage = rmd$devTime[keep], publication = rmd$publication[keep],
                      original_label = rmd$sub_rename_EML[keep]))
saveRDS(ref, file.path(dirname(ref_rds), "zhao2024_singler_reference.rds"))
rm(ref_seu); invisible(gc())
ref_tab <- as.data.frame(table(label = factor(ref$label, levels = ref_levels), stage = ref$stage,
                               publication = ref$publication))
ref_tab <- ref_tab[ref_tab$Freq > 0, ]
write.csv(ref_tab, file.path(out_dir, "singler_reference_labels.csv"), row.names = FALSE)
message("Reference: ", ncol(ref), " cells, ", length(ref_levels), " labels",
        if (length(small)) paste0(" (dropped, < ", min_ref_cells, " cells: ", paste(small, collapse = ", "), ")") else "")

# ============================================================
# 4. Query
# ============================================================
message("Loading Seurat object: ", input_rds)
seu <- readRDS(input_rds)
if (test_cells > 0 && test_cells < ncol(seu)) {
  set.seed(1)
  seu <- subset(seu, cells = sample(colnames(seu), test_cells))
}
DefaultAssay(seu) <- "RNA"
seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
test_mat <- LayerData(seu, assay = "RNA", layer = "data")
n_cells  <- ncol(seu)
res_cols <- grep(paste0("^", gsub(".", "\\.", res_prefix, fixed = TRUE)), colnames(seu[[]]), value = TRUE)
res_grid <- sort(as.numeric(sub(res_prefix, "", res_cols, fixed = TRUE)))
if (!is.null(res_use)) res_grid <- res_grid[round(res_grid, 3) %in% round(res_use, 3)]
res_cols <- paste0(res_prefix, res_grid)
smp      <- seu[[sample_col, drop = TRUE]]
samples  <- sort(unique(smp))
sample_cols <- setNames(c("#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B2")[seq_along(samples)], samples)
common <- intersect(rownames(test_mat), rownames(ref))
message("Loaded object: ", n_cells, " cells; ", length(res_grid), " resolutions; ", length(common),
        " genes shared with the reference")

# ============================================================
# 5. SingleR per cell
# ============================================================
t0 <- Sys.time()
trained <- trainSingleR(ref, labels = ref$label, de.method = de_method, de.n = de_n,
                        test.genes = rownames(test_mat), BPPARAM = bp)
pred <- classifySingleR(test_mat, trained, BPPARAM = bp)
message(sprintf("Per-cell SingleR: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

cell_lab    <- pred$labels
cell_pruned <- pred$pruned.labels
seu$singler_label        <- factor(cell_lab, levels = ref_levels)
seu$singler_pruned       <- factor(ifelse(is.na(cell_pruned), "Pruned", cell_pruned), levels = c(ref_levels, "Pruned"))
seu$singler_lineage      <- factor(ifelse(is.na(cell_pruned), "Pruned", lineage_map[cell_pruned]),
                                   levels = c(unique(lineage_map[ref_levels]), "Pruned"))
seu$singler_score        <- apply(pred$scores, 1, max)
seu$singler_delta_next   <- pred$delta.next
for (l in colnames(pred$scores)) seu[[paste0("singler_score_", make.names(l))]] <- setNames(pred$scores[, l], colnames(seu))

# ============================================================
# 6. SingleR per cluster at every resolution
# ============================================================
cl_list <- list()
for (k in seq_along(res_cols)) {
  col <- res_cols[k]; r <- res_grid[k]
  cl  <- droplevels(factor(seu[[col, drop = TRUE]]))
  ncl <- as.numeric(table(cl)[levels(cl)])
  mm  <- Matrix::sparse.model.matrix(~ 0 + cl)
  agg <- as.matrix(test_mat %*% mm) %*% diag(1 / ncl, nrow = length(ncl))   # mean log expression per cluster
  dimnames(agg) <- list(rownames(test_mat), levels(cl))
  cp  <- classifySingleR(agg, trained, prune = FALSE)
  second <- sapply(seq_len(nrow(cp$scores)), function(i) {   # best-scoring label other than the final one
    v <- cp$scores[i, ]; v <- v[names(v) != cp$labels[i]]; names(v)[which.max(v)] })
  tab <- table(cl, factor(cell_pruned, levels = ref_levels))       # pruned cells (NA) not counted
  maj <- ifelse(apply(tab, 1, max) > 0, colnames(tab)[apply(tab, 1, which.max)], NA)
  comp <- prop.table(table(cl, factor(smp, levels = samples)), 1)
  cl_list[[as.character(r)]] <- data.frame(
    resolution = r, cluster = levels(cl), cells = ncl,
    label = cp$labels, lineage = unname(lineage_map[cp$labels]), second_label = second,
    delta_next = round(cp$delta.next, 3),
    majority_cell_label = maj, majority_share = round(100 * apply(tab, 1, max) / ncl, 1),
    pct_pruned = round(100 * as.numeric(tapply(is.na(cell_pruned), cl, mean)), 1),
    matrix(round(100 * comp, 1), ncol = length(samples), dimnames = list(NULL, paste0("pct_", samples))),
    check.names = FALSE, row.names = NULL)
  seu[[paste0("singler_res.", r)]] <- setNames(factor(cp$labels[match(as.character(cl), levels(cl))], levels = ref_levels),
                                               colnames(seu))
  message("Resolution ", r, ": ", nlevels(cl), " clusters; labels: ", paste(names(table(cp$labels)), table(cp$labels),
                                                                           sep = " x", collapse = ", "))
}
cl_all <- do.call(rbind, cl_list)
write.csv(cl_all, file.path(out_dir, "singler_cluster_labels.csv"), row.names = FALSE)
lab_cols <- paste0("singler_res.", res_grid)
md <- seu[[]]
saveRDS(md[, c("singler_label", "singler_pruned", "singler_lineage", "singler_score", "singler_delta_next",
               grep("^singler_score_", colnames(md), value = TRUE), lab_cols)],
        file.path(out_dir, "singler_metadata.rds"))

by_sample <- rbind(
  data.frame(level = "label", as.data.frame(table(call = md$singler_pruned, sample = md[[sample_col]]))),
  data.frame(level = "lineage", as.data.frame(table(call = md$singler_lineage, sample = md[[sample_col]]))))
by_sample$pct <- round(100 * by_sample$Freq / as.numeric(table(md[[sample_col]])[as.character(by_sample$sample)]), 2)
write.csv(by_sample, file.path(out_dir, "singler_cell_labels_by_sample.csv"), row.names = FALSE)

# ============================================================
# 7. Summary text
# ============================================================
sm <- file(file.path(out_dir, "singler_summary.txt"), open = "wt")
cat("==================================================\n", file = sm)
cat("SingleR annotation against the Zhao 2024 human embryo reference\n", file = sm)
cat("==================================================\n", file = sm)
cat("Input:", input_rds, "\n", file = sm)
cat("Reference:", ref_rds, "\n", file = sm)
cat("Generated on:", as.character(Sys.time()), "\n", file = sm)
cat("Cells:", n_cells, "| samples:", paste(names(table(smp)), table(smp), sep = " = ", collapse = ", "), "\n", file = sm)
cat("Reference cells:", ncol(ref), "(CS7 gastrula removed", if (keep_amnion) "except amnion" else "", ") | genes shared:",
    length(common), "\n", file = sm)
cat("SingleR: de.method", de_method, ", de.n", de_n, ", quantile 0.8, fine-tuning on, pruning on (per cell)\n", file = sm)
cat("Cluster label: SingleR on the mean log-normalised profile of each cluster\n\n", file = sm)
cat("Reference labels by stage:\n", file = sm)
capture.output(print(table(factor(ref$label, levels = ref_levels), factor(ref$stage, levels = stage_levels))),
               file = sm, append = TRUE)
cat("\nPer-cell labels (Pruned = low confidence), % of each sample:\n", file = sm)
capture.output(print(round(100 * prop.table(table(md$singler_pruned, md[[sample_col]]), 2), 1)), file = sm, append = TRUE)
cat("\nPer-cell lineages, % of each sample:\n", file = sm)
capture.output(print(round(100 * prop.table(table(md$singler_lineage, md[[sample_col]]), 2), 1)), file = sm, append = TRUE)
cat("\n--------------------------------------------------\n", file = sm)
cat("Cluster labels per resolution\n", file = sm)
cat("--------------------------------------------------\n", file = sm)
for (r in names(cl_list)) {
  cat("\n--- Resolution", r, "---\n", file = sm)
  capture.output(print(cl_list[[r]][, c("cluster", "cells", "label", "second_label", "delta_next", "majority_cell_label",
                                         "majority_share", "pct_pruned", paste0("pct_", samples))],
                       row.names = FALSE), file = sm, append = TRUE)
}
close(sm)

# ============================================================
# 8. Plots
# ============================================================
rt <- as.data.frame(table(label = factor(ref$label, levels = rev(ref_levels)), stage = factor(ref$stage, levels = stage_levels)))
p_ref <- ggplot(rt[rt$Freq > 0, ], aes(stage, label, fill = Freq)) + geom_tile(colour = "white") +
  geom_text(aes(label = Freq), size = 3) + scale_fill_gradient(low = "#DEEBF7", high = "#08519C", name = "cells") +
  theme_minimal() + labs(title = paste0("Reference labels by stage (", ncol(ref), " cells; CS7 gastrula removed",
                                        if (keep_amnion) ", amnion kept" else "", ")"), x = "Stage", y = NULL)

p_cells <- (DimPlot(seu, reduction = umap_name, group.by = "singler_pruned", cols = label_cols, raster = TRUE, raster.dpi = raster_dpi, pt.size = pt_size) +
              labs(title = "Per-cell SingleR label (grey = pruned)")) |
  (DimPlot(seu, reduction = umap_name, group.by = "singler_lineage", cols = lineage_cols, raster = TRUE, raster.dpi = raster_dpi, pt.size = pt_size) +
     labs(title = "Lineage"))
p_cells2 <- (DimPlot(seu, reduction = umap_name, group.by = sample_col, cols = sample_cols, raster = TRUE, raster.dpi = raster_dpi, pt.size = pt_size) +
               labs(title = sample_col)) |
  (FeaturePlot(seu, features = "singler_delta_next", reduction = umap_name, raster = TRUE, raster.dpi = raster_dpi, pt.size = pt_size) +
     scale_colour_viridis_c() + labs(title = "Delta to the next label (low = ambiguous)"))

# One UMAP per per-cell label (pruned labels, "Pruned" included): all cells in grey, the cells
# called as that label in colour on top (Pruned in black, as its own colour is grey).
um    <- Embeddings(seu, umap_name)
lab_n <- table(md$singler_pruned)
lab_n <- lab_n[lab_n > 0]
bg    <- data.frame(x = um[, 1], y = um[, 2])
label_page <- function(L) {
  sel <- which(md$singler_pruned == L)
  ggplot(bg, aes(x, y)) +
    ggrastr::rasterise(geom_point(colour = "grey85", size = 0.3), dpi = 300) +
    ggrastr::rasterise(geom_point(data = bg[sel, ], colour = if (L == "Pruned") "black" else label_cols[[L]],
                                  size = 0.8), dpi = 300) +
    theme_void() + theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
    labs(title = sprintf("%s: %s cells (%.1f%% of all cells)", L, format(length(sel), big.mark = ","),
                         100 * length(sel) / n_cells))
}

sc_df <- do.call(rbind, lapply(ref_levels[ref_levels %in% cell_lab], function(l)
  data.frame(assigned = l, reference = colnames(pred$scores),
             score = colMeans(pred$scores[cell_lab == l, , drop = FALSE]), n = sum(cell_lab == l))))
sc_df$assigned  <- factor(paste0(sc_df$assigned, " (", sc_df$n, ")"),
                          levels = rev(unique(paste0(sc_df$assigned, " (", sc_df$n, ")"))))
sc_df$reference <- factor(sc_df$reference, levels = ref_levels)
p_scores <- ggplot(sc_df, aes(reference, assigned, fill = score)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", score)), size = 2.4) + scale_fill_viridis_c(name = "mean\nscore") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Mean SingleR score (Spearman) of the cells given each label, against every reference label",
       x = "Reference label", y = "Assigned label (cells)")
dd <- data.frame(label = factor(cell_lab, levels = ref_levels), delta = pred$delta.next, pruned = is.na(cell_pruned))
p_delta <- ggplot(dd, aes(label, delta, fill = label)) + geom_violin(scale = "width", colour = NA) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = label_cols, guide = "none") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = paste0("Delta to the next label per cell (", round(100 * mean(dd$pruned), 1), "% pruned)"), x = NULL, y = "delta")

bs <- by_sample
bs$call <- factor(bs$call, levels = unique(c(ref_levels, lineage_map, "Pruned")))
p_by_sample <- ggplot(bs[bs$level == "label", ], aes(sample, pct, fill = call)) + geom_col() +
  scale_fill_manual(values = label_cols) + theme_minimal() + labs(title = "Per-cell labels by sample", x = NULL, y = "%", fill = NULL) |
  ggplot(bs[bs$level == "lineage", ], aes(sample, pct, fill = call)) + geom_col() +
  scale_fill_manual(values = lineage_cols) + theme_minimal() + labs(title = "Per-cell lineages by sample", x = NULL, y = "%", fill = NULL)

grid_plots <- lapply(res_grid, function(r)
  DimPlot(seu, reduction = umap_name, group.by = paste0("singler_res.", r), cols = label_cols, raster = TRUE, raster.dpi = raster_dpi, pt.size = pt_size) +
    NoAxes() + labs(title = paste0("res ", r, " (", nrow(cl_list[[as.character(r)]]), " clusters)")) +
    theme(legend.text = element_text(size = 6), legend.key.size = unit(0.3, "cm")))

res_page <- function(r) {
  cs <- cl_list[[as.character(r)]]
  cl <- factor(seu[[paste0(res_prefix, r), drop = TRUE]], levels = cs$cluster)
  p_cl <- DimPlot(seu, reduction = umap_name, group.by = paste0(res_prefix, r), raster = TRUE, raster.dpi = raster_dpi, pt.size = pt_size, label = TRUE, repel = TRUE) +
    NoLegend() + labs(title = paste0("Clusters (", nrow(cs), ")"))
  p_lab <- DimPlot(seu, reduction = umap_name, group.by = paste0("singler_res.", r), cols = label_cols, raster = TRUE, raster.dpi = raster_dpi, pt.size = pt_size) +
    labs(title = "SingleR cluster label")
  ylab <- setNames(paste0(cs$cluster, " (", cs$cells, ") ", cs$label, " | 2nd: ", cs$second_label), cs$cluster)
  cc <- as.data.frame(prop.table(table(cluster = cl, label = seu$singler_pruned), 1))
  cc$cluster <- factor(cc$cluster, levels = rev(cs$cluster))
  p_cc <- ggplot(cc, aes(100 * Freq, cluster, fill = label)) + geom_col() + scale_fill_manual(values = label_cols) +
    scale_y_discrete(labels = ylab) + theme_minimal() + theme(axis.text.y = element_text(size = 6)) +
    labs(title = "Per-cell labels within each cluster (axis: cluster label, runner-up)", x = "% of cells", y = NULL, fill = NULL)
  comp <- do.call(rbind, lapply(samples, function(s) data.frame(cluster = cs$cluster, sample = s, pct = cs[[paste0("pct_", s)]])))
  comp$cluster <- factor(comp$cluster, levels = rev(cs$cluster))
  p_comp <- ggplot(comp, aes(pct, cluster, fill = sample)) + geom_col() + scale_fill_manual(values = sample_cols) +
    theme_minimal() + theme(axis.text.y = element_text(size = 6)) + labs(title = "Sample composition", x = "%", y = NULL, fill = NULL)
  ((p_cl | p_lab) / ((p_cc | p_comp) + plot_layout(widths = c(3, 1)))) +
    plot_annotation(title = paste("Resolution", r), theme = title_theme)
}

lab_long <- do.call(rbind, lapply(seq_along(res_grid), function(k)
  data.frame(resolution = factor(res_grid[k], levels = res_grid), as.data.frame(table(label = md[[lab_cols[k]]])))))
p_lab_res <- ggplot(lab_long[lab_long$Freq > 0, ], aes(resolution, Freq, fill = label)) +
  geom_col(position = "fill") + scale_fill_manual(values = label_cols) +
  scale_y_continuous(labels = function(v) paste0(100 * v, "%")) + theme_minimal() +
  labs(title = "Share of cells per SingleR cluster label at each resolution", x = "Resolution", y = "Cells", fill = NULL)
sample_long <- do.call(rbind, lapply(seq_along(res_grid), function(k) {
  d <- as.data.frame(prop.table(table(label = md[[lab_cols[k]]], sample = smp), 2))
  d$resolution <- factor(res_grid[k], levels = res_grid); d
}))
p_lab_sample <- ggplot(sample_long, aes(sample, Freq, fill = label)) + geom_col() +
  scale_fill_manual(values = label_cols) + scale_y_continuous(labels = function(v) paste0(100 * v, "%")) +
  facet_wrap(~ resolution, ncol = 5, labeller = label_both) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6)) +
  labs(title = "SingleR cluster label composition by sample at each resolution", x = NULL, y = "Cells", fill = NULL)

pdf(file.path(out_dir, "singler_annotation.pdf"), width = 17, height = 11)
print(p_ref)
print(p_cells)
print(p_cells2)
print(p_scores)
print(p_delta)
print(p_by_sample)
print(wrap_plots(grid_plots, ncol = 5) + plot_annotation(title = "SingleR cluster labels at each resolution", theme = title_theme))
for (r in res_grid) print(res_page(r))
print(p_lab_res)
print(p_lab_sample)
invisible(dev.off())
pdf(file.path(out_dir, "singler_label_umaps.pdf"), width = 11, height = 9)
for (L in names(lab_n)) print(label_page(L))
invisible(dev.off())
message("All outputs written to ", out_dir)
