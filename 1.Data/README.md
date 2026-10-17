# Data is downloaded here

Input data for 2.Scripts/run_pipeline.sh (not tracked by git):

- `<Sample>-raw_feature_bc_matrix/`: raw 10X matrix per sample (`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`),
  for 24h-blastoids, 48h-blastoids and 72h-blastoids
- `references/zhao2024_human_embryo/`: Zhao et al. 2024 human embryo reference (Nat Methods, doi 10.1038/s41592-024-02493-2)
  - `psd.R3.6.em.seurat.ob.rds`: integrated reference Seurat object, used by Annotation_SingleR.R
  - `zhao2024_celltype_markers_top25.csv`: top 25 markers per cell type (Supplementary Data 3), used by Slingshot.R
