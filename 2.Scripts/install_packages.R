# Installs the R packages used by the sc-Blastoids scripts (only those not installed yet)
#
# Sources: CRAN, Bioconductor and GitHub (monocle3, SeuratWrappers, Flufftail)
# Requires: R >= 4.4, a C/C++ compiler (Rtools on Windows, Xcode command line tools on macOS) and the
#           HDF5 library for hdf5r (e.g. brew install hdf5, or apt install libhdf5-dev)
# Usage:    Rscript 2.Scripts/install_packages.R   (from any directory)
#
# Versions the results in 3.Results were made with (R 4.4.2) are listed in `versions_used` and
# printed next to the installed versions at the end.

options(repos = c(CRAN = "https://cloud.r-project.org"))

have <- function(p) nzchar(system.file(package = p))
install_missing <- function(pkgs, installer) {
  missing <- pkgs[!vapply(pkgs, have, logical(1))]
  if (length(missing) > 0) {
    message("Installing: ", paste(missing, collapse = ", "))
    installer(missing)
  }
}

cran_pkgs <- c("remotes", "BiocManager", "Seurat", "scCustomize", "hdf5r", "Matrix", "ggplot2", "ggridges",
               "patchwork", "scales", "dplyr", "clustree", "ggrastr", "ggrepel", "gprofiler2", "ClustAssess",
               "RhpcBLASctl", "doParallel", "foreach", "mgcv")
bioc_pkgs <- c("SingleR", "SingleCellExperiment", "BiocParallel", "slingshot", "TrajectoryUtils",
               # monocle3 dependencies (monocle3 installation instructions)
               "BiocGenerics", "DelayedArray", "DelayedMatrixStats", "limma", "lme4", "S4Vectors",
               "SummarizedExperiment", "batchelor", "HDF5Array", "terra")
github_pkgs <- c(monocle3 = "cole-trapnell-lab/monocle3",
                 SeuratWrappers = "satijalab/seurat-wrappers",
                 Flufftail = "Core-Bioinformatics/Flufftail@9fd1d8b")

versions_used <- c(Seurat = "5.1.0", SeuratObject = "5.0.2", scCustomize = "3.3.0", hdf5r = "1.3.11",
                   Matrix = "1.7.1", ggplot2 = "3.5.1", ggridges = "0.5.6", patchwork = "1.3.0", scales = "1.3.0",
                   dplyr = "1.1.4", clustree = "0.5.1", ggrastr = "1.0.2", ggrepel = "0.9.6", gprofiler2 = "0.2.4",
                   ClustAssess = "1.2.0", RhpcBLASctl = "0.23.42", doParallel = "1.0.17", foreach = "1.5.2",
                   mgcv = "1.9.1", uwot = "0.2.2", SingleR = "2.8.0", SingleCellExperiment = "1.28.1",
                   BiocParallel = "1.40.0", slingshot = "2.14.0", TrajectoryUtils = "1.14.0",
                   monocle3 = "1.4.27", SeuratWrappers = "0.4.0", Flufftail = "0.0.1")

install_missing(c("remotes", "BiocManager"), install.packages)
install_missing(setdiff(cran_pkgs, c("remotes", "BiocManager")), install.packages)
install_missing(bioc_pkgs, function(p) BiocManager::install(p, update = FALSE, ask = FALSE))
for (p in names(github_pkgs)) {
  if (!have(p)) {
    message("Installing ", p, " from GitHub (", github_pkgs[[p]], ")")
    remotes::install_github(github_pkgs[[p]], upgrade = "never", dependencies = NA)
  }
}

# Installed versions next to the versions used for 3.Results
installed <- vapply(names(versions_used), function(p)
  if (have(p)) as.character(packageVersion(p)) else "MISSING", character(1))
print(data.frame(package = names(versions_used), installed = installed, used = versions_used, row.names = NULL),
      row.names = FALSE)
if (any(installed == "MISSING")) stop("Some packages could not be installed: ",
                                      paste(names(installed)[installed == "MISSING"], collapse = ", "))
message("All packages installed")
