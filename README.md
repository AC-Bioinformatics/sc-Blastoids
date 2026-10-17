# Single-cell transcriptomics analysis of naive hPSC derived blastoids

Single-cell RNA-seq analysis of human blastoids (stem cell-derived models of the blastocyst) at three time points:
24 h, 48 h and 72 h. The samples were profiled with 10x Genomics Flex (probe-based fixed RNA profiling) and pooled in a
single capture with probe sample barcodes, so there is no technical batch between the time points; there is one sample
per time point.

The analysis goes from the raw 10X matrices to cell types and trajectories:

- ambient RNA removal (CellBender) and quality control per sample
- log normalisation, PCA and UMAP
- Louvain clustering at resolutions 0.1–1.5 and clustering stability (ClustAssess)
- cell type annotation with SingleR against the Zhao et al. 2024 human embryo reference
- fuzzy clustering, pseudotime and the transition between the two 24 h epiblast clusters (Flufftail)
- lineage trajectories from the epiblast (Slingshot)

## Repository layout

```
1.Data/       input data (not tracked by git; see 1.Data/README.md)
2.Scripts/    one script per step, run in order by run_pipeline.sh; shared helpers in functions.R
3.Results/    results, one folder per step (figures, tables and reports; large objects are not tracked)
```

## Requirements

- R ≥ 4.4 with a C/C++ compiler, and the HDF5 library (e.g. `brew install hdf5`, `apt install libhdf5-dev`)
- For the full dataset (~105,000 cells), about 64 GB of RAM; see [Runtime and memory](#runtime-and-memory)
- Only to start from the raw data: CellBender (see step 3 below); a GPU is strongly recommended

## How to run

1. **Get the data.** The data are shared on request. Put the shared files in `1.Data/`, check and unzip them:

   ```
   cd 1.Data
   md5sum -c md5sums.txt        # macOS: md5 *.zip, and compare with md5sums.txt
   unzip '*.zip'
   cd ..
   ```

   `1.Data/README.md` lists the files and the layout the scripts expect.

2. **Install the R packages** (only those not installed yet):

   ```
   Rscript 2.Scripts/install_packages.R
   ```

3. **Install CellBender** (only to start from the raw data); recommended in a conda environment:

   ```
   conda create -n cellbender python=3.12
   conda activate cellbender
   pip install cellbender
   ```

   Alternatively, run it in its container with `2.Scripts/cellbender_container.sh` (Apptainer or Docker).

4. **Run the pipeline**, either

   - from the merged object after CellBender (`1.Data/merged_filtered_seurat.rds`), which skips CellBender:

     ```
     bash 2.Scripts/run_pipeline.sh --from qc
     ```

   - or from the raw 10X matrices, including CellBender:

     ```
     bash 2.Scripts/run_pipeline.sh                          # with the cellbender conda environment active
     CELLBENDER=apptainer bash 2.Scripts/run_pipeline.sh     # or in the container (apptainer or docker)
     ```

   Each step writes its results to `3.Results/` and a log to `3.Results/logs/<step>.log`. A stopped run can be continued
   from any step with `--from <step>`. Single scripts can also be run on their own, from any directory, e.g.
   `Rscript 2.Scripts/QC.R`.

## Pipeline steps

| Step | Script | Output in `3.Results/` |
|---|---|---|
| `cellbender` | `cellbender.sh` (or `cellbender_container.sh`) | `2.CellBender/` |
| `h5_to_seurat` | `h5_to_seurat.R` | `2.CellBender/merged_filtered_seurat.rds` |
| `qc` | `QC.R` | `3.QC_filtering/` |
| `normalisation` | `Normalisation.R` | `5.Normalisation/` |
| `clustering` | `Clustering.R` | `6.Clustering/Louvain/` |
| `annotation` | `Annotation_SingleR.R` | `7.Annotation/singler_zhao2024/` |
| `clustassess` | `ClustAssess.R`, `ClustAssess_postprocess.R` | `6.Clustering/ClustAssess/` |
| `flufftail` | `Flufftail.R` | `8.Pseudotime/flufftail/res0.3/` |
| `slingshot` | `Slingshot.R` | `8.Pseudotime/slingshot/` |

The Louvain clusters are saved in the repository (`3.Results/6.Clustering/Louvain/cluster_metadata.rds`); `Clustering.R`
reuses them when they cover every cell and only remakes the plots.

## Options

Set as environment variables in front of the command, e.g. `N_CORES=8 bash 2.Scripts/run_pipeline.sh --from qc`:

| Variable | Default | Effect |
|---|---|---|
| `N_CORES` | 4 | parallel workers for SingleR and ClustAssess |
| `CELLBENDER` | `local` | `local` (pip install), `apptainer` or `docker` |
| `SC_RESULTS` | `3.Results` | results folder, e.g. `SC_RESULTS=test` for a separate run |
| `SC_RECLUSTER` | `false` | `true` recomputes the Louvain clusters instead of reusing the saved ones |

## Runtime and memory

Measured on the full dataset on a 16-core server:

| Step | Time | Peak memory |
|---|---|---|
| CellBender | ~35 min per sample on a GPU (many hours on a CPU) | |
| h5_to_seurat | 5 min | 24 GB |
| QC | 5 min | 26 GB |
| Normalisation | 11 min | 38 GB |
| Clustering | 10 min with the saved clusters (19 min reclustering) | 24 GB |
| SingleR annotation | 5 min | 46 GB |
| ClustAssess | 7 min | 27 GB |
| Flufftail | 30 min | 63 GB |
| Slingshot | 19 min | 46 GB |

After CellBender, the whole pipeline takes about 1.5–2 hours.

## Reference

Zhao C, et al. A comprehensive human embryo reference tool using single-cell RNA-sequencing data.
Nat Methods (2024). doi:10.1038/s41592-024-02493-2
