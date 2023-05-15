import os,glob

import numpy as np
import pandas as pd
import shutil

ALIGNMENT_DIR = 'alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13'
SHARE_DIR = os.path.join(ALIGNMENT_DIR,'share')

files_to_move = ['raw_feature_bc_matrix.h5',
                 'filtered_feature_bc_matrix.h5',
                 'metrics_summary.csv',
                 'web_summary.html']

for csv_dir in glob.glob(os.path.join(ALIGNMENT_DIR,'*.csv')):
    ds_name = os.path.basename(csv_dir).split('.')[0]
    DS_DIR = os.path.join(SHARE_DIR,ds_name)
    os.makedirs(DS_DIR,exist_ok=True)
    #
    df_ds = pd.read_csv(csv_dir)
    assert np.all(df_ds.sample_id.value_counts() == 1)
    for _,(sample_id,file) in df_ds.iterrows():
        dest = os.path.join(SHARE_DIR,ds_name,sample_id)
        src_prefix = os.path.dirname(file)
        os.makedirs(dest,exist_ok=True)
        for file_mv in files_to_move:
            src  = os.path.join(src_prefix,file_mv)
            shutil.copy(src,dest)
