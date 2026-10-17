# Data is downloaded here

The input data are not tracked by git; they are shared on request as a Google Drive folder with these files:

| File | Contents | Used by |
|---|---|---|
| `24h-blastoids-raw_feature_bc_matrix.zip` | raw 10X matrix, 24h blastoids | `cellbender.sh` |
| `48h-blastoids-raw_feature_bc_matrix.zip` | raw 10X matrix, 48h blastoids | `cellbender.sh` |
| `72h-blastoids-raw_feature_bc_matrix.zip` | raw 10X matrix, 72h blastoids | `cellbender.sh` |
| `merged_filtered_seurat.zip` | merged Seurat object after CellBender and before QC (`3.Results/2.CellBender/merged_filtered_seurat.rds`), to start at QC without running CellBender | `QC.R` |
| `zhao2024_reference.zip` | Zhao et al. 2024 human embryo reference (Nat Methods, doi 10.1038/s41592-024-02493-2): `psd.R3.6.em.seurat.ob.rds` and the top 25 markers per cell type, `zhao2024_celltype_markers_top25.csv` (Supplementary Data 3) | `Annotation_SingleR.R`, `Slingshot.R` |
| `md5sums.txt` | checksums of the zip files | |

Check the downloads (`md5sum -c md5sums.txt`; on macOS `md5 *.zip` and compare with `md5sums.txt`), then unzip
the raw matrices and the reference here, in `1.Data/`:

```
unzip '*-raw_feature_bc_matrix.zip'
unzip zhao2024_reference.zip
```

This gives the layout the scripts expect:

```
1.Data/
  24h-blastoids-raw_feature_bc_matrix/   barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
  48h-blastoids-raw_feature_bc_matrix/
  72h-blastoids-raw_feature_bc_matrix/
  references/zhao2024_human_embryo/      psd.R3.6.em.seurat.ob.rds, zhao2024_celltype_markers_top25.csv
```

## Where to start

- From the raw data (runs CellBender, which needs a GPU or many hours on a CPU):
  `bash 2.Scripts/run_pipeline.sh`
- From the merged object after CellBender, skipping CellBender and `h5_to_seurat.R`: unzip
  `merged_filtered_seurat.zip` in the repository root (it contains `3.Results/2.CellBender/merged_filtered_seurat.rds`), then
  `bash 2.Scripts/run_pipeline.sh --from qc`

The Louvain clusters are in the repository (`3.Results/6.Clustering/Louvain/cluster_metadata.rds`) and are reused by
`Clustering.R`, so they are not recomputed.
