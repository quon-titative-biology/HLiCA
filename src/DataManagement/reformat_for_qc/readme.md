# Upload to box

This module will be used to prepare the aligned read files for sharing (h5ad)

## Copy metadata

```
python 01_cellranger_aggr_create_csv.py
```

Create {dataset_name}.csv that contains sample_id and molecule_h5 columns for cellranger aggr

## Manually copy files for Mullen dataset

```
bash 04_copy_mullen.sh
```

## Copy the count matrices to shared folder

```
python 03_copy_h5_to_share.py
```

# Run QC pipeline

```
R modified_tallulah_pipeline.R
R combine_qc_outputs.R
```
