### Version 2 of the Slingshot analysis. Changes from version 1:
###   - Morula and Prelineage cells are removed (drop_labels): Prelineage looked like noise
###     leaking between states, and the blastoids show no naive / morula signal.
###   - slingshot(shrink = FALSE): lineages that share a start are not pulled towards a
###     common average curve before they branch, so the weak Early PrE branch (gain 0.018
###     in version 1) cannot bend the trophoblast lineage. Each lineage is fitted on its
###     own cells.
###   Everything else is as in version 1.
###
### Slingshot pseudotime of the APAP blastoid data, using the per-cell SingleR labels
### (Zhao et al. 2024 human embryo reference, with the Early / Mature TE split) as the
### groups for the lineage tree.
###
### Sources:
###   - slingshot 2.14.0 vignette and help pages (vignette.Rmd, ?slingshot,
###     ?slingPseudotime, ?embedCurves) and TrajectoryUtils ?createClusterMST.
###     Street K, et al. (2018) BMC Genomics 19:477.
###   - Gene trends: the idea of tradeSeq (Van den Berge K, et al. (2020) Nat Commun
###     11:1201), a smooth GAM of expression along pseudotime per gene. tradeSeq is not
###     installed, so mgcv is used directly: Gaussian GAM on log-normalised expression,
###     expr ~ s(pseudotime, k = gam_k), per lineage, on a subsample of cells.
###
### Steps:
###   1. All cells with a SingleR label (pruned.labels); Pruned cells and labels with
###      fewer than min_label_cells cells are left out (a group needs cells to have a
###      centre and spread).
###   2. PCA recomputed on these cells as in the normalisation
###      (log-normalised data, n_features variable genes (vst), ScaleData, no
###      regression), n_pcs PCs; UMAP on the same PCs, for plots only.
###   3. slingshot(PCs 1:n_pcs, SingleR labels, start.clus = start_clus) fitted on a
###      random subset of fit_cells cells (fitting on all ~93,000 cells took over an hour
###      and did not finish); the other cells are projected onto the fitted curves with
###      predict() (slingshot vignette, "projecting cells onto existing trajectories").
###      End points are not fixed. Defaults otherwise (dist.method "slingshot",
###      approx_points 150).
###   4. Checks: pseudotime by sample (24h / 48h / 72h; not used to build the
###      trajectory, so an independent check) and by label; marker genes along
###      pseudotime: a hand-picked panel (marker_genes) and the Zhao et al. 2024 top 25
###      markers per cell type (Supplementary Data 3), summarised per cell type as the
###      mean z-scored trend and where it peaks.
###   5. Genes changing along each lineage: GAM per variable gene, BH-adjusted
###      p-values, deviance explained; heatmap of the top genes ordered by peak, with
###      dashed lines where pseudotime passes from one SingleR label to the next
###      (halfway between the labels' median pseudotimes).
###      Fitted serially (about 0.1 s per gene with 5,000 cells): forked parallel
###      workers (BiocParallel MulticoreParam) hung after the PCA / UMAP steps.
###
### Input: 3.Results/6.Clustering/Louvain/seurat_multi_res_clustered.rds,
###        3.Results/7.Annotation/singler_zhao2024/singler_metadata.rds,
###        1.Data/references/zhao2024_human_embryo/zhao2024_celltype_markers_top25.csv
###
### Usage: Rscript 2.Scripts/Slingshot.R [--test]   (from any directory)
###   --test  5,000 random cells (2,000 fitted, 3,000 projected), written to <out_dir>/test
###
### Outputs (3.Results/8.Pseudotime/slingshot/):
###   slingshot_pseudotime.pdf      tree, curves on UMAP and PCA, pseudotime by sample and
###                                 label, marker gene trends, heatmaps of top genes
###   slingshot_summary.txt         lineages, tree edges, checks, top genes per lineage
###   lineage_genes_<lineage>.csv   GAM results per gene for each lineage
###   zhao_marker_trends.csv        mean z-scored trend of each Zhao cell type's markers
###   zhao_marker_peaks.csv         pseudotime at which each Zhao cell type's trend peaks
###   pseudotime_metadata.rds       per-cell label, pseudotime and weight per lineage,
###                                 average pseudotime, subset PCA/UMAP coordinates
###   slingshot_pto.rds             the slingshot PseudotimeOrdering object
library(Seurat)
library(slingshot)
library(mgcv)
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

input_rds  <- file.path(Results, "6.Clustering", "Louvain", "seurat_multi_res_clustered.rds")
singler_md <- file.path(Results, "7.Annotation", "singler_zhao2024", "singler_metadata.rds")
zhao_top25 <- file.path(Root, "1.Data", "references", "zhao2024_human_embryo", "zhao2024_celltype_markers_top25.csv")
out_dir    <- file.path(Results, "8.Pseudotime", "slingshot")

# ============================================================
# 2. Parameters
# ============================================================
label_col       <- "singler_pruned"  # per-cell SingleR label, "Pruned" = low confidence
start_clus      <- "Early EPI"       # root of the tree (blastoid TE is induced from naive PSCs)
drop_labels     <- c("Morula", "Prelineage")  # removed entirely (not fitted, not projected)
sling_shrink    <- FALSE             # slingshot(shrink): FALSE = no averaging of shared starts
min_label_cells <- 30                # leave out labels with fewer cells (must exceed n_pcs, see below)
sample_col      <- "Sample"
n_features      <- 3000              # variable genes, as in the normalisation
n_pcs           <- 20                # PCs for slingshot and UMAP, as in the clustering
fit_cells       <- 20000             # random cells used to fit the curves; the rest are projected
min_fit_label   <- 200               # each label keeps at least this many cells in the fit (or all)
gam_cells       <- 5000              # cells per lineage for the gene GAMs
gam_k           <- 6                 # basis size of the smooth
gam_min_cells   <- 50                # skip genes detected in fewer cells of the subsample
top_genes       <- 40                # genes per lineage in the heatmap
n_grid          <- 100               # pseudotime points for predicted trends
test_cells      <- 0                 # > 0: run on a random subset of this many cells
pt_size         <- 3
raster_dpi      <- c(800, 800)
title_theme     <- theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

marker_genes <- list(
  "Pluripotency / epiblast" = c("POU5F1", "NANOG", "SOX2", "KLF17", "DPPA3"),
  "Pan-TE / nascent TE"     = c("GATA3", "GATA2", "TFAP2C", "KRT18", "CDX2"),
  "Early TE"                = c("ESRRB", "ENPEP"),
  "Mature TE"               = c("NR2F2", "CGA", "GCM1", "VGLL1", "OVOL1"),
  "CTB"                     = c("TP63", "ELF5", "TEAD4", "PAGE4"),
  "STB"                     = c("CGB3", "CGB5", "CGB8", "SDC1", "CYP19A1", "HSD3B1", "KISS1"),
  "PrE / hypoblast"         = c("GATA6", "SOX17", "GATA4", "PDGFRA"))

# Zhao et al. 2024 cell types whose top 25 markers are followed along pseudotime: those on
# the path. Left out: gastrula types, Z4cell / 8 cell, EVT, ExE_Mes, Amnion (no cells here).
zhao_types <- c("Morula", "Prelineage", "ICM", "Epiblast", "Hypoblast", "TE", "CTB", "STB")
zhao_cols  <- c(Morula = "#B08968", Prelineage = "#7F5539", ICM = "#F4A261", Epiblast = "#D62728",
                Hypoblast = "#2CA02C", TE = "#3182BD", CTB = "#17BECF", STB = "#9467BD")

label_cols <- c(Cleavage = "#D9C2A7", Morula = "#B08968", Prelineage = "#7F5539", ICM = "#F4A261",
                "EPI/PrE intermediate" = "#E9C46A", "Early EPI" = "#D62728", "Late EPI" = "#8B0000",
                "Early PrE" = "#2CA02C", "Late PrE" = "#145A32", "Early TE" = "#6BAED6",
                "Mature TE" = "#1F4E9A", CTB = "#17BECF", STB = "#9467BD", EVT = "#E377C2",
                "ExE mesoderm" = "#BCBD22", Amnion = "#FF7F0E")

args <- commandArgs(trailingOnly = TRUE)
if ("--test" %in% args) {
  test_cells <- 5000
  fit_cells  <- 2000     # so that the projection step is also tested
  gam_cells  <- 2000
  out_dir    <- file.path(out_dir, "test")
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
t_start <- Sys.time()

# ============================================================
# 3. Cells and labels
# ============================================================
message("Loading Seurat object: ", input_rds)
seu <- readRDS(input_rds)
sr  <- readRDS(singler_md)
seu$label <- as.character(sr[colnames(seu), label_col])

n_all <- ncol(seu)
cells <- colnames(seu)[!is.na(seu$label) & !seu$label %in% c("Pruned", drop_labels)]
if (test_cells > 0 && test_cells < length(cells)) {
  set.seed(1)
  cells <- sample(cells, test_cells)
}
n_lab    <- table(seu$label[match(cells, colnames(seu))])
small    <- names(n_lab)[n_lab < min_label_cells]
keep     <- colnames(seu) %in% cells & !seu$label %in% small   # Pruned and drop_labels are not in cells
excluded <- table(seu$label[!keep])
if (test_cells > 0) excluded <- table(seu$label[colnames(seu) %in% cells & !keep])
seu <- subset(seu, cells = colnames(seu)[keep])
invisible(gc())
lab_levels <- intersect(names(label_cols), unique(seu$label))
seu$label  <- factor(seu$label, levels = lab_levels)
smp   <- seu[[sample_col, drop = TRUE]]
hours <- as.numeric(sub("^([0-9]+)h.*", "\\1", smp))   # "24h-blastoids" -> 24
message("Cells used: ", ncol(seu), " of ", n_all, "; left out: ",
        paste(names(excluded), excluded, sep = " = ", collapse = ", "))

# ============================================================
# 4. PCA and UMAP on these cells
# ============================================================
DefaultAssay(seu) <- "RNA"
seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = n_features, verbose = FALSE)
seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = n_pcs, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:n_pcs, reduction.name = "umap_sub", reduction.key = "umapsub_",
               seed.use = 42, verbose = FALSE)
message("PCA and UMAP done: ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min")

# ============================================================
# 5. Slingshot
# ============================================================
pcs <- Embeddings(seu, "pca")[, 1:n_pcs]
# Random subset, topped up so every label keeps at least min_fit_label cells: Slingshot's
# distance between labels uses each label's covariance in PC space, which is singular
# ("system is computationally singular") when a label has fewer cells than PCs
lab_all <- as.character(seu$label)
set.seed(1)
if (nrow(pcs) > fit_cells) {
  fit_idx <- sample(nrow(pcs), fit_cells)
  for (idx in split(seq_len(nrow(pcs)), lab_all)) {
    rest <- setdiff(idx, fit_idx)
    need <- min(length(idx), min_fit_label) - (length(idx) - length(rest))
    if (need > 0) fit_idx <- c(fit_idx, rest[sample.int(length(rest), need)])
  }
  fit_idx <- sort(fit_idx)
} else fit_idx <- seq_len(nrow(pcs))
fit_lab <- lab_all[fit_idx]
message("Fitting Slingshot on ", length(fit_idx), " cells: ",
        paste(names(table(fit_lab)), table(fit_lab), sep = " = ", collapse = ", "))
t0  <- Sys.time()
pto <- slingshot(pcs[fit_idx, ], clusterLabels = fit_lab, start.clus = start_clus, shrink = sling_shrink)
t_sling <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
saveRDS(pto, file.path(out_dir, "slingshot_pto.rds"))

lineages <- slingLineages(pto)
lin_names <- names(lineages)
lin_desc  <- setNames(sapply(lineages, paste, collapse = " > "), lin_names)
message("Slingshot: ", length(lineages), " lineages in ", t_sling, " min")
for (l in lin_names) message("  ", l, ": ", lin_desc[l])

# Pseudotime and weights for every cell: fitted cells from the fit, the others projected
pt  <- matrix(NA_real_, nrow(pcs), length(lin_names), dimnames = list(rownames(pcs), lin_names))
wts <- pt
pt[fit_idx, ]  <- slingPseudotime(pto, na = TRUE)[, lin_names]   # NA = not on that lineage
wts[fit_idx, ] <- slingCurveWeights(pto)[, lin_names]
t0 <- Sys.time()
if (length(fit_idx) < nrow(pcs)) {
  prj <- predict(pto, pcs[-fit_idx, , drop = FALSE])
  pt[-fit_idx, ]  <- slingPseudotime(prj, na = TRUE)[, lin_names]
  wts[-fit_idx, ] <- slingCurveWeights(prj)[, lin_names]
}
t_pred <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
avg <- rowSums(pt * wts, na.rm = TRUE) / rowSums(wts)         # weighted average, as slingAvgPseudotime
avg[!is.finite(avg)] <- NA
fitted <- seq_len(nrow(pcs)) %in% fit_idx
message("Projected ", sum(!fitted), " cells in ", t_pred, " min")

mst   <- slingMST(pto)
edges <- igraph::as_edgelist(mst)
edge_attr <- function(a) { v <- igraph::edge_attr(mst, a); if (is.null(v)) NA_real_ else round(v, 3) }
edges <- data.frame(from = edges[, 1], to = edges[, 2], length = edge_attr("weight"), gain = edge_attr("gain"))
edges <- edges[!edges$from %in% ".OMEGA" & !edges$to %in% ".OMEGA", ]

for (l in lin_names) {
  seu[[paste0("pt_", l)]] <- pt[, l]
  seu[[paste0("w_", l)]]  <- wts[, l]
}
seu$pt_avg <- avg
seu$fitted <- fitted

md_out <- data.frame(label = seu$label, sample = smp, fitted = fitted, pt_avg = avg, pt, check.names = FALSE)
colnames(md_out)[-(1:4)] <- paste0("pt_", lin_names)
md_out <- cbind(md_out, setNames(as.data.frame(wts), paste0("w_", lin_names)),
                Embeddings(seu, "umap_sub"), Embeddings(seu, "pca")[, 1:n_pcs])
saveRDS(md_out, file.path(out_dir, "pseudotime_metadata.rds"))

# ============================================================
# 6. Checks: pseudotime by sample and by label
# ============================================================
by_sample <- do.call(rbind, lapply(lin_names, function(l) {
  ok <- !is.na(pt[, l])
  data.frame(lineage = l, sample = sort(unique(smp)),
             cells = as.numeric(table(factor(smp[ok], levels = sort(unique(smp))))),
             median_pt = round(sapply(sort(unique(smp)), function(s) median(pt[ok & smp == s, l])), 2))
}))
rho <- sapply(lin_names, function(l) {
  ok <- !is.na(pt[, l])
  round(cor(pt[ok, l], hours[ok], method = "spearman"), 3)
})
by_label <- do.call(rbind, lapply(lin_names, function(l) {
  ok <- !is.na(pt[, l])
  tb <- table(seu$label[ok])
  d  <- data.frame(lineage = l, label = names(tb), cells = as.numeric(tb))
  d  <- d[d$cells > 0, ]
  d$median_pt <- round(sapply(d$label, function(x) median(pt[ok & seu$label == x, l])), 2)
  d[order(d$median_pt), ]
}))
# Same label, fitted vs projected cells: similar medians mean the subset did not bias the curves
fit_vs_prj <- do.call(rbind, lapply(lin_names, function(l) {
  ok <- !is.na(pt[, l])
  lv <- intersect(lab_levels, unique(as.character(seu$label[ok])))
  med <- function(sel) round(sapply(lv, function(x) median(pt[ok & sel & seu$label == x, l])), 2)
  n   <- function(sel) sapply(lv, function(x) sum(ok & sel & seu$label == x))
  data.frame(lineage = l, label = lv, fitted_cells = n(fitted), fitted_median = med(fitted),
             projected_cells = n(!fitted), projected_median = med(!fitted))
}))

# ============================================================
# 7. Gene trends along each lineage (GAM per gene)
# ============================================================
data_mat <- LayerData(seu, assay = "RNA", layer = "data")
genes    <- VariableFeatures(seu)

fit_genes <- function(E, x, grid) {
  res <- lapply(seq_len(nrow(E)), function(i) {
    y <- E[i, ]
    if (sum(y > 0) < gam_min_cells) return(rep(NA_real_, 3 + length(grid)))
    f <- gam(y ~ s(x, k = gam_k), method = "REML")
    s <- summary(f)
    c(s$s.table[1, "p-value"], s$dev.expl, s$s.table[1, "F"],
      as.numeric(predict(f, newdata = data.frame(x = grid))))
  })
  do.call(rbind, res)
}

gene_res <- list()
gene_pred <- list()
grids <- list()
for (l in lin_names) {
  t0 <- Sys.time()
  on_l <- which(!is.na(pt[, l]))
  set.seed(1)
  sub  <- if (length(on_l) > gam_cells) sort(sample(on_l, gam_cells)) else on_l
  x    <- pt[sub, l]
  grid <- seq(min(x), max(x), length.out = n_grid)
  E    <- as.matrix(data_mat[genes, sub, drop = FALSE])
  r    <- fit_genes(E, x, grid)
  pred <- r[, -(1:3), drop = FALSE]
  rownames(pred) <- genes
  d <- data.frame(gene = genes, p_value = r[, 1], dev_explained = round(r[, 2], 4), F = round(r[, 3], 2),
                  start = round(pred[, 1], 3), end = round(pred[, n_grid], 3),
                  peak_pt = round(grid[max.col(replace(pred, is.na(pred), -Inf), ties.method = "first")], 2))
  d$peak_pt[is.na(d$p_value)] <- NA
  d$p_adj  <- p.adjust(d$p_value, method = "BH")
  d$change <- round(d$end - d$start, 3)
  d <- d[order(-d$dev_explained), c("gene", "p_value", "p_adj", "dev_explained", "F", "start", "end", "change", "peak_pt")]
  write.csv(d, file.path(out_dir, paste0("lineage_genes_", l, ".csv")), row.names = FALSE)
  gene_res[[l]] <- d; gene_pred[[l]] <- pred; grids[[l]] <- grid
  message("GAMs ", l, ": ", length(sub), " cells, ", sum(!is.na(d$p_value)), " genes, ",
          sum(d$p_adj < 0.05, na.rm = TRUE), " with adj. p < 0.05 (",
          round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min)")
}

# Predicted trend (GAM) of a set of genes along one lineage, on the same cell subsample
# as above; genes detected in fewer than 10 cells of the subsample give NA
trend_matrix <- function(gs, l) {
  on_l <- which(!is.na(pt[, l]))
  set.seed(1)
  sub  <- if (length(on_l) > gam_cells) sort(sample(on_l, gam_cells)) else on_l
  x    <- pt[sub, l]
  grid <- seq(min(x), max(x), length.out = n_grid)
  E    <- as.matrix(data_mat[gs, sub, drop = FALSE])
  m <- do.call(rbind, lapply(gs, function(g) {
    y <- E[g, ]
    if (sum(y > 0) < 10) return(rep(NA_real_, n_grid))
    as.numeric(predict(gam(y ~ s(x, k = gam_k), method = "REML"), newdata = data.frame(x = grid)))
  }))
  rownames(m) <- gs
  list(grid = grid, m = m)
}

# Hand-picked marker panel (not limited to variable genes)
mk_all <- unique(unname(unlist(marker_genes)))
mk_all <- mk_all[mk_all %in% rownames(data_mat)]
mk_trend <- do.call(rbind, lapply(lin_names, function(l) {
  tm <- trend_matrix(mk_all, l)
  data.frame(lineage = l, gene = rep(mk_all, each = n_grid), pt = rep(tm$grid, length(mk_all)),
             expr = as.numeric(t(tm$m)))
}))
mk_trend$gene <- factor(mk_trend$gene, levels = mk_all)

# Zhao et al. 2024 top 25 markers per cell type: each gene's trend is z-scored along the
# lineage, then averaged per cell type; the peak of that average places the cell type
zt <- read.csv(zhao_top25, check.names = FALSE)
zhao_markers <- lapply(setNames(zhao_types, zhao_types), function(t) {
  g <- zt[[t]]
  g[!is.na(g) & g != "" & g %in% rownames(data_mat)]
})
zhao_genes <- unique(unlist(zhao_markers, use.names = FALSE))
z_type <- data.frame(type = rep(names(zhao_markers), lengths(zhao_markers)),
                     gene = unlist(zhao_markers, use.names = FALSE))
zhao_trend  <- list()
zhao_module <- list()
for (l in lin_names) {
  tm <- trend_matrix(zhao_genes, l)
  zm <- t(scale(t(tm$m)))                       # z-score per gene along pseudotime
  zhao_trend[[l]] <- list(grid = tm$grid, z = zm)
  zhao_module[[l]] <- do.call(rbind, lapply(zhao_types, function(t) {
    g <- zhao_markers[[t]]
    data.frame(lineage = l, type = t, pt = tm$grid,
               score = colMeans(zm[g, , drop = FALSE], na.rm = TRUE),
               genes_fitted = sum(!is.na(zm[g, 1])))
  }))
}
zhao_module <- do.call(rbind, zhao_module)
zhao_module$type <- factor(zhao_module$type, levels = zhao_types)
zhao_peak <- do.call(rbind, lapply(split(zhao_module, list(zhao_module$lineage, zhao_module$type), drop = TRUE), function(d) {
  if (all(is.na(d$score))) return(NULL)
  k <- which.max(d$score)
  data.frame(lineage = d$lineage[1], type = d$type[1], genes_fitted = d$genes_fitted[1],
             peak_pt = round(d$pt[k], 2), peak_pct = round(100 * (d$pt[k] - min(d$pt)) / diff(range(d$pt)), 1))
}))
zhao_peak <- zhao_peak[order(zhao_peak$lineage, zhao_peak$peak_pt), ]
rownames(zhao_peak) <- NULL
write.csv(zhao_module, file.path(out_dir, "zhao_marker_trends.csv"), row.names = FALSE)
write.csv(zhao_peak, file.path(out_dir, "zhao_marker_peaks.csv"), row.names = FALSE)
message("Zhao marker trends: ", length(zhao_genes), " genes, ", length(zhao_types), " cell types")

# ============================================================
# 8. Summary text
# ============================================================
sm <- file(file.path(out_dir, "slingshot_summary.txt"), open = "wt")
cat("==================================================\n", file = sm)
cat("Slingshot pseudotime with SingleR (Zhao 2024) labels, version 2\n", file = sm)
cat("==================================================\n", file = sm)
cat("Generated on:", as.character(Sys.time()), "\n", file = sm)
cat("Cells used:", ncol(seu), "of", n_all, if (test_cells > 0) "(test subset)" else "", "\n", file = sm)
cat("Left out:", paste(names(excluded), excluded, sep = " = ", collapse = ", "), "\n", file = sm)
cat("PCA on these cells:", n_features, "variable genes, PCs 1 -", n_pcs, "used\n", file = sm)
cat("Version 2: removed", paste(drop_labels, collapse = ", "), "; shrink =", sling_shrink, "\n", file = sm)
cat("Slingshot: start", start_clus, ", end points free, defaults otherwise; fitted on", length(fit_idx),
    "random cells in", t_sling, "min; other", sum(!fitted), "cells projected with predict() in", t_pred, "min\n\n", file = sm)
cat("Cells per label:\n", file = sm)
capture.output(print(table(seu$label, smp)), file = sm, append = TRUE)
cat("\nLineages:\n", file = sm)
for (l in lin_names) cat("  ", l, ": ", lin_desc[l], " (", sum(!is.na(pt[, l])), " cells)\n", sep = "", file = sm)
cat("\nTree edges (length = distance between label centres; gain = confidence, higher = clearer):\n", file = sm)
capture.output(print(edges, row.names = FALSE), file = sm, append = TRUE)
cat("\nSpearman correlation of pseudotime with sample time (24 / 48 / 72h):\n", file = sm)
capture.output(print(rho), file = sm, append = TRUE)
cat("\nMedian pseudotime by sample:\n", file = sm)
capture.output(print(by_sample, row.names = FALSE), file = sm, append = TRUE)
cat("\nMedian pseudotime by label (ordered), per lineage:\n", file = sm)
capture.output(print(by_label, row.names = FALSE), file = sm, append = TRUE)
cat("\nFitted vs projected cells, median pseudotime per label:\n", file = sm)
capture.output(print(fit_vs_prj, row.names = FALSE), file = sm, append = TRUE)
for (l in lin_names) {
  d <- gene_res[[l]]
  d <- d[!is.na(d$p_adj) & d$p_adj < 0.05, ]
  cat("\n", l, " (", lin_desc[l], "): ", nrow(d), " of ", length(genes),
      " variable genes change along pseudotime (adj. p < 0.05). Top 20 by deviance explained:\n", sep = "", file = sm)
  capture.output(print(head(d[, c("gene", "dev_explained", "start", "end", "change", "peak_pt")], 20), row.names = FALSE),
                 file = sm, append = TRUE)
}
cat("\nZhao et al. 2024 top 25 markers per cell type: where the mean z-scored trend peaks\n",
    "(peak_pct = % of the way along the lineage; genes_fitted = markers detected in >= 10 cells):\n", file = sm)
capture.output(print(zhao_peak, row.names = FALSE), file = sm, append = TRUE)
cat("\nTotal runtime:", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), "min\n", file = sm)
close(sm)

# ============================================================
# 9. Plots
# ============================================================
curve_df <- function(obj, dims = 1:2) {
  cvs <- slingCurves(obj)                     # unnamed list, same order as slingLineages()
  do.call(rbind, lapply(seq_along(cvs), function(i) {
    s <- cvs[[i]]$s[cvs[[i]]$ord, dims, drop = FALSE]
    data.frame(lineage = lin_names[i], x = s[, 1], y = s[, 2])
  }))
}
add_curves <- function(p, cd) {
  p + geom_path(data = cd, aes(x, y, group = lineage), colour = "black", linewidth = 0.8, inherit.aes = FALSE)
}
umap_sub <- Embeddings(seu, "umap_sub")
cd_umap  <- curve_df(embedCurves(pto, umap_sub[fit_idx, ]))
cd_pca   <- curve_df(pto, 1:2)
lin_cols <- setNames(scales::hue_pal()(length(lin_names)), lin_names)
dim_args <- list(raster = TRUE, raster.dpi = raster_dpi, pt.size = pt_size)

# Tree: label centres on the UMAP joined by the tree edges
cent <- aggregate(as.data.frame(umap_sub), list(label = seu$label), median)
colnames(cent)[2:3] <- c("x", "y")
ed <- edges
ed[, c("x", "y")]       <- cent[match(ed$from, cent$label), c("x", "y")]
ed[, c("xend", "yend")] <- cent[match(ed$to, cent$label), c("x", "y")]
p_tree <- do.call(DimPlot, c(list(seu, reduction = "umap_sub", group.by = "label", cols = label_cols), dim_args)) +
  geom_segment(data = ed, aes(x, y, xend = xend, yend = yend), inherit.aes = FALSE, linewidth = 1) +
  geom_label(data = cent, aes(x, y, label = label), inherit.aes = FALSE, size = 3) +
  labs(title = paste0("Slingshot tree between SingleR labels (start: ", start_clus, ")"))

p_curves <- add_curves(do.call(DimPlot, c(list(seu, reduction = "umap_sub", group.by = "label", cols = label_cols), dim_args)), cd_umap) +
  labs(title = "Slingshot curves on the UMAP of these cells") |
  do.call(DimPlot, c(list(seu, reduction = "umap_sub", group.by = sample_col), dim_args)) + labs(title = "Sample")
p_avg <- add_curves(do.call(FeaturePlot, c(list(seu, features = "pt_avg", reduction = "umap_sub"), dim_args)), cd_umap) +
  scale_colour_viridis_c(option = "plasma") + labs(title = "Average pseudotime") |
  add_curves(do.call(DimPlot, c(list(seu, reduction = "pca", group.by = "label", cols = label_cols), dim_args)), cd_pca) +
  labs(title = "Curves in PCA space (PC1, PC2)")

p_lin <- lapply(lin_names, function(l)
  add_curves(do.call(FeaturePlot, c(list(seu, features = paste0("pt_", l), reduction = "umap_sub"), dim_args)),
             cd_umap[cd_umap$lineage == l, ]) +
    scale_colour_viridis_c(option = "plasma", na.value = "grey90") +
    labs(title = l, subtitle = lin_desc[l]) + theme(plot.subtitle = element_text(size = 8)))

pt_long <- do.call(rbind, lapply(lin_names, function(l) {
  ok <- !is.na(pt[, l])
  data.frame(lineage = l, pt = pt[ok, l], sample = smp[ok], label = seu$label[ok])
}))
p_by_sample <- ggplot(pt_long, aes(pt, fill = sample)) + geom_density(alpha = 0.5, colour = NA) +
  facet_wrap(~ lineage, scales = "free", labeller = labeller(lineage = function(v) paste0(v, ": rho = ", rho[v]))) +
  theme_minimal() + labs(title = "Pseudotime by sample (rho = Spearman correlation with sample time)",
                         x = "Pseudotime", y = "Density", fill = NULL)
p_by_label <- ggplot(pt_long, aes(pt, reorder(label, pt, FUN = median), fill = label)) +
  geom_boxplot(outlier.shape = NA) + scale_fill_manual(values = label_cols, guide = "none") +
  facet_wrap(~ lineage, scales = "free_x") + theme_minimal() +
  labs(title = "Pseudotime by SingleR label (ordered by median)", x = "Pseudotime", y = NULL)

p_markers <- ggplot(mk_trend[!is.na(mk_trend$expr), ], aes(pt, expr, colour = lineage)) + geom_line(linewidth = 0.8) +
  scale_colour_manual(values = lin_cols) + facet_wrap(~ gene, scales = "free_y", ncol = 8) +
  theme_minimal() + theme(strip.text = element_text(size = 8)) +
  labs(title = "Marker genes along pseudotime (GAM fit, log-normalised expression)", x = "Pseudotime", y = "Expression")

p_zhao_module <- ggplot(zhao_module[!is.na(zhao_module$score), ], aes(pt, score, colour = type)) +
  geom_line(linewidth = 1) + scale_colour_manual(values = zhao_cols) +
  facet_wrap(~ lineage, scales = "free_x", labeller = labeller(lineage = function(v) paste0(v, ": ", lin_desc[v]))) +
  theme_minimal() + theme(strip.text = element_text(size = 7)) +
  labs(title = "Zhao et al. 2024 top 25 markers per cell type: mean z-scored trend along each lineage",
       x = "Pseudotime", y = "Mean z-score of the cell type's markers", colour = "Zhao cell type")

# Where pseudotime passes from one SingleR label to the next along a lineage: labels of the
# lineage placed at their median pseudotime (all cells on the lineage); boundaries halfway
# between consecutive medians. Drawn as dashed lines, with the labels on a top axis.
label_bounds <- function(l) {
  ok  <- !is.na(pt[, l])
  med <- sapply(lineages[[l]], function(x) median(pt[ok & seu$label == x, l]))
  med <- sort(med[!is.na(med)])
  list(pos = data.frame(label = names(med), pt = as.numeric(med)),
       bounds = (head(med, -1) + tail(med, -1)) / 2)
}
add_bounds <- function(p, l, xr) {
  b   <- label_bounds(l)
  bd  <- b$bounds[b$bounds >= xr[1] & b$bounds <= xr[2]]
  pos <- b$pos[b$pos$pt >= xr[1] & b$pos$pt <= xr[2], ]
  p + geom_vline(xintercept = bd, linetype = "dashed", colour = "grey20", linewidth = 0.4) +
    scale_x_continuous(sec.axis = dup_axis(breaks = pos$pt, labels = pos$label, name = NULL)) +
    theme(axis.text.x.top = element_text(angle = 40, hjust = 0, vjust = 0, size = 9, face = "bold"),
          axis.ticks.x.top = element_line())
}

zhao_heat <- function(l) {
  zl <- zhao_trend[[l]]
  d  <- z_type[!is.na(zl$z[z_type$gene, 1]), ]
  if (nrow(d) < 2) return(NULL)
  d$peak <- max.col(zl$z[d$gene, , drop = FALSE], ties.method = "first")
  d  <- d[order(match(d$type, zhao_types), d$peak), ]
  d$row <- paste(d$type, d$gene, sep = ": ")
  hd <- data.frame(row = factor(rep(d$row, each = n_grid), levels = rev(d$row)),
                   type = factor(rep(d$type, each = n_grid), levels = zhao_types),
                   pt = rep(zl$grid, nrow(d)), z = as.numeric(t(zl$z[d$gene, , drop = FALSE])))
  ggplot(hd, aes(pt, row, fill = z)) + geom_tile() +
    facet_grid(type ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = function(v) sub("^.*: ", "", v)) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", name = "z") + theme_minimal() +
    theme(axis.text.y = element_text(size = 3.5), strip.text.y = element_text(angle = 0)) +
    labs(title = paste0(l, ": Zhao et al. 2024 top 25 markers per cell type along pseudotime (z-score per gene)"),
         subtitle = lin_desc[l], x = "Pseudotime", y = NULL) -> p
  add_bounds(p, l, range(zl$grid))
}

heat_page <- function(l) {
  d <- gene_res[[l]]
  d <- d[!is.na(d$p_adj) & d$p_adj < 0.05, ]
  g <- head(d$gene, top_genes)
  if (length(g) < 2) return(NULL)
  m <- gene_pred[[l]][g, , drop = FALSE]
  m <- t(scale(t(m)))
  g <- g[order(max.col(m, ties.method = "first"))]
  hd <- data.frame(gene = factor(rep(g, each = n_grid), levels = rev(g)), pt = rep(grids[[l]], length(g)),
                   z = as.numeric(t(m[g, ])))
  ggplot(hd, aes(pt, gene, fill = z)) + geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", name = "z") + theme_minimal() +
    theme(axis.text.y = element_text(size = 6)) +
    labs(title = paste0(l, ": top ", length(g), " genes changing along pseudotime (ordered by peak)"),
         subtitle = lin_desc[l], x = "Pseudotime", y = NULL) -> p
  add_bounds(p, l, range(grids[[l]]))
}

pdf(file.path(out_dir, "slingshot_pseudotime.pdf"), width = 17, height = 11)
print(p_tree)
print(p_curves)
print(p_avg)
print(wrap_plots(p_lin, ncol = min(3, length(p_lin))) + plot_annotation(title = "Pseudotime per lineage (grey = not on the lineage)", theme = title_theme))
print(p_by_sample)
print(p_by_label)
print(p_markers)
print(p_zhao_module)
for (l in lin_names) { h <- zhao_heat(l); if (!is.null(h)) print(h) }
for (l in lin_names) { h <- heat_page(l); if (!is.null(h)) print(h) }
invisible(dev.off())
message("All outputs written to ", out_dir, " (", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min)")
