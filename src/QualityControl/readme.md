Run Quality Control using Tallulah's pipeline.

Two major QC performed are:
1. EmptyDrops - remove empty droplets
2. SoupX - remove background RNA

After each QC, we run standard Seurat pipeline + basic cell type annotation with scmap.

The Seurat files saved are
Cleaned_output_EmptyOnly.rds
Cleaned_output_SoupX.rds
where the final QC output should be Cleaned_output_SoupX.rds

However, for intergration, Cleaned_output_EmptyOnly.rds is used since it still provides integer for its count matrix.
