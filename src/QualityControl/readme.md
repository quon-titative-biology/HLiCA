# Quality control



## Set up environment

```
bash create_R_env.sh
```

## Move files to QC directory

Run from reformat_for_sharing in DataManagement

```
python 01_cellranger_aggr_create_csv.py
python 03_copy_h5_to_share.py
```

## Run Tallulah's pipeline
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

```
R TallulahPipeline/modified_tallulah_pipeline.R
R TallulahPipeline/combine_qc_outputs.R
python TallulahPipeline/combine_metadata.py
```

For upgrading the metadata files with alignment and QC (GoogleSheetMetadata_sample_qced.csv)
```
python src/DataManagement/reformat_for_sharing
/reformat_sample_metadata.py
```

For additional stats, run following script

```
python plot_stats.py
```
