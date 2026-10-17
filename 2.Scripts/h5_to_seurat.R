# This script converts the cellbender output .h5 files to a merged v5 seurat object for downstream processing. 
# Inputs -- 3.Results/2.CellBender/<Sample>_cellbender_filtered.h5 for each sample (from cellbender.sh). 
# Outputs -- 3.Results/2.CellBender/merged_filtered_seurat.rds, merged seurat v5 object with each sample as a different assay / layer,
#            plus seurat_processing_log.txt and seurat_layer_stats.txt. 
# Usage -- Rscript 2.Scripts/h5_to_seurat.R (from any directory)

# 1. SET PATHS RELATIVE TO THE REPOSITORY ROOT
# Repository root: the folder above 2.Scripts (from the script's own path with Rscript or source(),
# otherwise the working directory)
Root <- local({
  f <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), fixed = TRUE)
  if (length(f) == 0) f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (length(f) == 1) normalizePath(file.path(dirname(f), "..")) else getwd()
})
if (!dir.exists(file.path(Root, "2.Scripts"))) stop("Repository root not found (no 2.Scripts folder in ", Root, ")")

out_dir <- file.path(Root, "3.Results", "2.CellBender")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

output_rds <- file.path(out_dir, "merged_filtered_seurat.rds")
log_filepath <- file.path(out_dir, "seurat_processing_log.txt")

log_file <- file(log_filepath, open = "wt")
sink(log_file, type = "message")
sink(log_file, type = "output")

library(Seurat)
library(scCustomize)  # Read_CellBender_h5_Multi_File / Read_CellBender_h5_Mat
library(hdf5r)        # needed by scCustomize to read the h5 files
library(Matrix)

cat("==================================================\n")
cat("Starting Seurat Processing Workflow\n")
cat("Date/Time:", as.character(Sys.time()), "\n")
cat("==================================================\n\n")

# 2. CONFIGURATION
data_dir <- out_dir

# 3. READ ALL H5 FILES AT ONCE
# NOTE: adjust custom_name to match your actual filename suffix exactly
# (whatever gsub needs to strip to leave the bare sample ID).
raw_data_list <- Read_CellBender_h5_Multi_File(
  data_dir    = data_dir,
  custom_name = "_filtered.h5",   # <-- confirm this matches your real suffix
  merge       = FALSE              # keep per-sample matrices, not one pooled matrix
)

if (length(raw_data_list) == 0) {
  stop("Error: No matching h5 files found in directory: ", data_dir)
}

cat("Found", length(raw_data_list), "files to process.\n")
names(raw_data_list)

seurat_list <- list()

# 4. BUILD A SEURAT OBJECT + QC METRICS FOR EACH SAMPLE, WITH PROGRESS REPORTING
n_samples  <- length(raw_data_list)
start_time <- Sys.time()

for (i in seq_along(raw_data_list)) {

  sample_id    <- names(raw_data_list)[i]
  sample_start <- Sys.time()

  cat("\n----------------------------------------\n")
  cat("[", i, "/", n_samples, "] Processing sample:", sample_id, "\n")
  cat("  Started:", format(sample_start, "%H:%M:%S"), "\n")

  seurat_obj <- CreateSeuratObject(
    counts  = raw_data_list[[sample_id]],
    project = sample_id
  )

  seurat_obj$Sample <- sample_id

  n_before <- ncol(seurat_obj)

  # Filter cells: keep droplets with at least 200 genes
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= 200)

  n_after  <- ncol(seurat_obj)
  pct_kept <- round(100 * n_after / n_before, 1)

  # QC metrics
  seurat_obj$log1p_total_counts <- log1p(seurat_obj$nCount_RNA)
  seurat_obj[["pct_counts_mt"]]   <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["pct_counts_ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL")
  seurat_obj[["pct_counts_hb"]]   <- PercentageFeatureSet(seurat_obj, pattern = "^HB[^(P)]")

  seurat_list[[sample_id]] <- seurat_obj

  sample_elapsed <- round(difftime(Sys.time(), sample_start, units = "secs"), 1)
  total_elapsed  <- round(difftime(Sys.time(), start_time, units = "mins"), 1)

  cat("  Cells:", n_before, "->", n_after, "(", pct_kept, "% retained )\n")
  cat("  Median genes/cell:", median(seurat_obj$nFeature_RNA),
      "| Median UMIs/cell:", median(seurat_obj$nCount_RNA), "\n")
  cat("  Mean %MT:", round(mean(seurat_obj$pct_counts_mt), 2), "\n")
  cat("  Sample time:", sample_elapsed, "sec | Total elapsed:", total_elapsed, "min\n")
  cat("  Progress:", round(100 * i / n_samples, 1), "% complete\n")
}

cat("\n==================================================\n")
cat("All samples processed. Total time:",
    round(difftime(Sys.time(), start_time, units = "mins"), 1), "min\n")
cat("==================================================\n")

# 5. MERGE AND SAVE
cat("\n----------------------------------------\n")
cat("Merging all datasets into one object...\n")

merged_seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[2:length(seurat_list)],
  add.cell.ids = names(seurat_list),
  project = "Merged_data"
)

cat("Saving merged object to:", output_rds, "\n")
saveRDS(merged_seurat, file = output_rds)

cat("\nExtracting layer structures and metadata statistics...\n")

# Define the statistics log path
stats_filepath <- file.path(out_dir, "seurat_layer_stats.txt")
stats_conn <- file(stats_filepath, open = "wt")

cat("==================================================\n", file = stats_conn)
cat("SEURAT OBJECT LAYER METADATA & STATISTICS\n", file = stats_conn)
cat("Generated on:", as.character(Sys.time()), "\n", file = stats_conn)
cat("==================================================\n\n", file = stats_conn)

# --- Overall object class ---
cat("Object class:", paste(class(merged_seurat), collapse = ", "), "\n", file = stats_conn)
cat("Seurat object version:", as.character(merged_seurat@version), "\n\n", file = stats_conn)

# 1. Identify all assays present in the object (usually just "RNA")
assays_present <- Assays(merged_seurat)
cat("Assays present in object:", paste(assays_present, collapse = ", "), "\n\n", file = stats_conn)

for (assay_name in assays_present) {
  cat("----------------------------------------\n", file = stats_conn)
  cat("Assay Name:", assay_name, "\n", file = stats_conn)
  cat("Assay class:", paste(class(merged_seurat[[assay_name]]), collapse = ", "), "\n", file = stats_conn)
  cat("----------------------------------------\n", file = stats_conn)

  # Feature-level metadata for this assay (e.g. annotation columns in meta.features)
  feature_meta <- merged_seurat[[assay_name]][[]]
  feature_meta_cols <- colnames(feature_meta)
  cat("Feature-level metadata columns:",
      if (length(feature_meta_cols) > 0) paste(feature_meta_cols, collapse = ", ") else "(none)",
      "\n\n", file = stats_conn)

  # 2. Get all layer names inside this specific assay
  layers_present <- Layers(merged_seurat[[assay_name]])
  cat("Layers found:", paste(layers_present, collapse = ", "), "\n\n", file = stats_conn)

  # 3. Loop through each layer to pull dimensions, class, and zero-inflated statistics
  for (layer_name in layers_present) {
    cat("  ▶ Layer:", layer_name, "\n", file = stats_conn)

    layer_mat <- GetAssayData(merged_seurat, assay = assay_name, layer = layer_name)

    # Calculate dimensional traits
    n_features <- nrow(layer_mat)
    n_cells    <- ncol(layer_mat)

    # Calculate basic matrix fill statistics safely using sparse matrix behaviors
    total_elements <- as.numeric(n_features) * as.numeric(n_cells)
    non_zero_count <- length(layer_mat@x) # Quick access to non-zero slots in dgCMatrix
    sparsity_pct   <- (1 - (non_zero_count / total_elements)) * 100

    cat("    - Class              :", paste(class(layer_mat), collapse = ", "), "\n", file = stats_conn)
    cat("    - Features (Genes)   :", n_features, "\n", file = stats_conn)
    cat("    - Cells (Barcodes)   :", n_cells, "\n", file = stats_conn)
    cat("    - Sparsity (Zeros)   :", sprintf("%.2f%%", sparsity_pct), "\n\n", file = stats_conn)
  }
}

# 4. Save Global Cell Metadata Summary Checklist
cat("----------------------------------------\n", file = stats_conn)
cat("Global Metadata Column Summaries\n", file = stats_conn)
cat("----------------------------------------\n", file = stats_conn)

meta_cols <- colnames(merged_seurat[[]])
cat("Cell-level metadata columns (", length(meta_cols), "):\n", file = stats_conn)
cat(paste0("  - ", meta_cols, collapse = "\n"), "\n\n", file = stats_conn)

cat("Total Cells in Master Object:", ncol(merged_seurat), "\n", file = stats_conn)
cat("Total Genes in Master Object:", nrow(merged_seurat), "\n\n", file = stats_conn)

# Print a snapshot breakdown of cell allocations across your custom Sample column
sample_breakdown <- table(merged_seurat$Sample)
cat("Cell Count Breakdown Per Sample Group:\n", file = stats_conn)
for (s_name in names(sample_breakdown)) {
  cat("  -", s_name, ":", sample_breakdown[[s_name]], "cells\n", file = stats_conn)
}

close(stats_conn)
cat("Layer statistics successfully saved to:", stats_filepath, "\n")

cat("Finished! All samples successfully merged into a seurat object.\n")

# Turn off log file recording and close the file connection safely
sink(type = "message")
sink(type = "output")
close(log_file)