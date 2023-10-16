# Upload to box

This module will be used to prepare the aligned read files for sharing (h5ad)

## Update the metadata files to include the h5ad information

```
python 01_cellranger_aggr_create_csv.py
```

Create {dataset_name}.csv that contains sample_id and molecule_h5 columns for cellranger aggr

## Aggregate the outputs

This script uses the .csv file create from the previous step to run cellranger aggr function.

```
02_cellranger_aggr.sh
```

## Manually modify aggr csv for GuilliamsScott

ADT should not be combined with others
