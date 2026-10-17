# Helper functions for the sc-Blastoids scripts (sourced by QC.R)

log_msg <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", ...)

or_default <- function(a, b) if (is.null(a)) b else a


# Thresholds ----
mad_bounds <- function(x, n_mads, constant = 1.4826) {
  if (is.null(n_mads)) return(c(lower = -Inf, upper = Inf))
  med <- median(x, na.rm = TRUE)
  m <- mad(x, na.rm = TRUE, constant = constant)
  c(lower = med - n_mads * m, upper = med + n_mads * m)
}

# TRUE when a metric takes more than one value across cells
has_variation <- function(x) length(unique(x[!is.na(x)])) > 1

hard_bounds <- function(lower, upper) c(lower = or_default(lower, -Inf), upper = or_default(upper, Inf))

# Open interval for MAD and nFeature bounds; upper caps (nCount, mito) keep values <= cap
in_bounds <- function(x, bounds) x > bounds["lower"] & x < bounds["upper"]
below_cap <- function(x, cap) if (is.null(cap)) rep(TRUE, length(x)) else x <= cap


# QC metrics ----
# Percentage of counts in the n most expressed genes of each cell
pct_top_n <- function(mat, n = 20) {
  mat <- as(mat, "CsparseMatrix")
  p <- mat@p
  x <- mat@x
  top_sum <- vapply(seq_len(ncol(mat)), function(j) {
    if (p[j + 1] == p[j]) return(0)
    vals <- x[(p[j] + 1):p[j + 1]]
    sum(sort(vals, decreasing = TRUE)[seq_len(min(n, length(vals)))])
  }, numeric(1))
  setNames(100 * top_sum / Matrix::colSums(mat), colnames(mat))
}

# QC metrics from one counts matrix (genes x cells)
qc_metrics <- function(mat, top_n = 20) {
  total <- Matrix::colSums(mat)
  pct <- function(pattern) {
    idx <- grep(pattern, rownames(mat))
    100 * Matrix::colSums(mat[idx, , drop = FALSE]) / total
  }
  data.frame(
    nCount_RNA = total,
    nFeature_RNA = Matrix::colSums(mat > 0),
    pct_counts_mt = pct("^MT-"),
    pct_counts_ribo = pct("^RPS|^RPL"),
    pct_counts_hb = pct("^HB[^(P)]"),
    pct_top = pct_top_n(mat, top_n),
    row.names = colnames(mat)
  )
}


# Per-sample histograms ----
# Histogram of one metric, kept (green) vs removed (red) cells, with dashed threshold lines
qc_histogram <- function(df, var, vlines = list(), bins = 100, title = var, xlab = var, xlim = NULL) {
  plot_df <- data.frame(value = df[[var]], status = df$qc_status)
  constant <- !has_variation(plot_df$value)
  hist_layer <- if (constant) {
    geom_histogram(binwidth = 0.01, boundary = plot_df$value[1], alpha = 0.85, color = NA, position = "stack")
  } else if (!is.null(xlim)) {
    geom_histogram(binwidth = diff(xlim) / bins, boundary = xlim[1], alpha = 0.85, color = NA, position = "stack")
  } else {
    geom_histogram(bins = bins, alpha = 0.85, color = NA, position = "stack")
  }
  p <- ggplot(plot_df, aes(x = value, fill = status)) +
    hist_layer +
    scale_fill_manual(values = c(Kept = "#2ca02c", Removed = "#d62728"), name = "QC status") +
    theme_minimal() +
    labs(title = title, x = xlab, y = "Number of cells")
  # A constant metric is drawn as a narrow bar on a 0-1 axis
  if (constant) {
    v <- plot_df$value[1]
    p <- p + coord_cartesian(xlim = c(min(0, v), max(0, v) + 1)) +
      labs(subtitle = paste0("All cells: ", xlab, " = ", signif(v, 4)))
  } else if (!is.null(xlim)) {
    p <- p + coord_cartesian(xlim = xlim)
  }
  line_labels <- c()
  for (nm in names(vlines)) {
    vl <- vlines[[nm]]
    if (!is.null(vl$x) && is.finite(vl$x)) {
      p <- p + geom_vline(xintercept = vl$x, color = vl$color, linewidth = 0.9, linetype = "dashed")
      line_labels <- c(line_labels, paste0(nm, " = ", signif(vl$x, 4)))
    }
  }
  if (length(line_labels) > 0) p <- p + labs(caption = paste(line_labels, collapse = "   |   "))
  p
}


# Cross-sample QC plots ----
# Exponential saturation fit per sample on the log10 scale:
#   log10_nFeature = Asym + (R0 - Asym) * exp(-exp(lrc) * log10_nCount)
fit_saturation <- function(md) {
  fits <- lapply(split(md, md$Sample), function(d) {
    fit <- tryCatch(
      nls(log10_nFeature ~ SSasymp(log10_nCount, Asym, R0, lrc), data = d),
      error = function(e) { log_msg("  Saturation fit failed for ", d$Sample[1], ": ", conditionMessage(e)); NULL }
    )
    if (is.null(fit)) return(NULL)
    grid <- data.frame(log10_nCount = seq(min(d$log10_nCount), max(d$log10_nCount), length.out = 200))
    grid$log10_nFeature <- predict(fit, newdata = grid)
    grid$Sample <- d$Sample[1]
    list(curve = grid, coef = coef(fit))
  })
  fits[!vapply(fits, is.null, logical(1))]
}

ridge_plot <- function(md, var, label) {
  ggplot(md, aes(x = .data[[var]], y = Sample, fill = Sample)) +
    geom_density_ridges(alpha = 0.7, scale = 1.2) +
    geom_vline(xintercept = median(md[[var]], na.rm = TRUE), color = "red", linewidth = 0.8, linetype = "dashed") +
    theme_ridges() + theme(legend.position = "none") +
    labs(title = paste(label, "distribution per sample"), x = var, y = NULL)
}

# Cross-sample QC plots (violin, ridge, scatter, saturation) and the saturation fit coefficients
qc_combined_plots <- function(obj) {
  md <- obj[[]]
  md$log10_nCount <- log10(md$nCount_RNA + 1)
  md$log10_nFeature <- log10(md$nFeature_RNA + 1)

  vln_features <- intersect(c("nFeature_RNA", "nCount_RNA", "pct_counts_mt", "pct_counts_hb", "pct_counts_ribo"), colnames(md))

  scatter_feat_count <- ggplot(md, aes(x = nCount_RNA, y = nFeature_RNA, color = Sample)) +
    geom_point(size = 0.4, alpha = 0.4) + facet_wrap(~ Sample) +
    theme_minimal() + theme(legend.position = "none") +
    labs(title = "nFeature_RNA vs nCount_RNA")

  fits <- fit_saturation(md)
  curves <- do.call(rbind, lapply(fits, `[[`, "curve"))
  fit_caption <- paste(vapply(names(fits), function(s) {
    cf <- fits[[s]]$coef
    sprintf("%s: Asym = %.3f, R0 = %.3f, lrc = %.3f", s, cf["Asym"], cf["R0"], cf["lrc"])
  }, character(1)), collapse = "   |   ")
  scatter_log_log <- ggplot(md, aes(x = log10_nCount, y = log10_nFeature, color = Sample)) +
    geom_point(size = 0.4, alpha = 0.4) + facet_wrap(~ Sample) +
    theme_minimal() + theme(legend.position = "none") +
    labs(title = "log10(nFeature_RNA) vs log10(nCount_RNA)",
         subtitle = "Exponential saturation fit: y = Asym + (R0 - Asym) * exp(-exp(lrc) * x)",
         caption = fit_caption)
  if (!is.null(curves)) {
    scatter_log_log <- scatter_log_log +
      geom_line(data = curves, aes(x = log10_nCount, y = log10_nFeature), color = "black", linewidth = 0.5, inherit.aes = FALSE)
  }

  list(
    plots = list(
      violin = VlnPlot(obj, features = vln_features, group.by = "Sample", ncol = 3, pt.size = 0.1, raster = TRUE, alpha = 0.01),
      ridge_nCount = ridge_plot(md, "nCount_RNA", "nCount_RNA"),
      ridge_mito = ridge_plot(md, "pct_counts_mt", "Mitochondrial %"),
      ridge_nFeature = ridge_plot(md, "nFeature_RNA", "nFeature_RNA"),
      scatter_feat_count = scatter_feat_count,
      scatter_log_log = scatter_log_log
    ),
    fits = lapply(fits, `[[`, "coef")
  )
}

# All cross-sample QC plots in one PDF, after a title page with cell counts
write_combined_qc_pdf <- function(qc, pdf_path, title, n_before, n_after) {
  pdf(pdf_path, width = 17, height = 11)
  plot.new()
  text(0.5, 0.7, title, cex = 2.5, font = 2)
  text(0.5, 0.6, paste("Generated:", Sys.time()), cex = 1.2)
  text(0.5, 0.5, paste("Total cells before:", n_before), cex = 1.2)
  text(0.5, 0.45, paste("Total cells after:", n_after), cex = 1.2)
  text(0.5, 0.4, paste("Retention rate:", round(100 * n_after / n_before, 1), "%"), cex = 1.2)
  for (p in qc$plots) print(p)
  invisible(dev.off())
}

# Four-panel QC figure: A violin plots, B nCount_RNA ridges, C nFeature_RNA ridges, D saturation fit
write_qc_figure <- function(qc, pdf_path) {
  # The violin grid is one panel; its margin leaves room for the tag
  violin <- wrap_elements(full = qc$plots$violin & theme(plot.margin = margin(40, 5, 5, 60)))
  panels <- c(list(violin), qc$plots[c("ridge_nCount", "ridge_nFeature", "scatter_log_log")])
  pdf(pdf_path, width = 34, height = 22)
  print(wrap_plots(panels, ncol = 2) +
    plot_annotation(tag_levels = "A", theme = theme(plot.margin = margin(20, 20, 20, 20))) &
    theme(plot.tag = element_text(size = 40, face = "bold")))
  invisible(dev.off())
}
