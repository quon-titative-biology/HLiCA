import os, glob

import numpy as np
import pandas as pd

FASTQ_DIR = 'fastqfiles'

DS_DIR = glob.glob(os.path.join(FASTQ_DIR,"*/GoogleSheetMetadata_file.csv"))

DS_names = [ds.split('/')[-2] for ds in DS_DIR]

DS_metadatas = [pd.read_csv(ds) for ds in DS_DIR]

def create_column_chart(DS_metadatas, DS_names):
    columns_per_ds = [df_meta.columns for df_meta in DS_metadatas]
    columns = np.unique(np.concatenate(columns_per_ds))
    columns_bool = {ds_name: np.isin(columns,col) for ds_name,col in zip(DS_names,columns_per_ds)}
    df_meta_stat = pd.DataFrame(columns_bool, index=columns).transpose()
    return df_meta_stat

def count_unique_columns(DS_metadatas,column):
    return [len(df_meta[column].unique()) for df_meta in DS_metadatas]

df_meta_stat = create_column_chart(DS_metadatas,DS_names)

for col in ['tissue_alias', 'library_alias', 'SAMPLE']:
    df_meta_stat['# unique ' + col] = count_unique_columns(DS_metadatas,col)

df_meta_stat.transpose()
