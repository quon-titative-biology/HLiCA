import os,glob

import numpy as np
import pandas as pd
import shutil
from distutils.dir_util import copy_tree


# Source folder that contains metadata csv
ALIGNMENT_DIR = 'alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13'
# Desitnation folder that will be shared
SHARE_DIR = os.path.join(ALIGNMENT_DIR,'share_qc')

datasets = ['DasGupta',
            'Gruen',
            'Henderson',
            'GuilliamsScott',
            'Toronto']

files_to_move = ['raw_feature_bc_matrix.h5',
                 'filtered_feature_bc_matrix.h5',
                 'metrics_summary.csv',
                 'web_summary.html']

files_to_move = []

folder_to_move = ['raw_feature_bc_matrix']


for ds_name in datasets:
    csv_dir = os.path.join(ALIGNMENT_DIR,f'{ds}.csv')
    DS_DIR = os.path.join(SHARE_DIR,ds_name)
    os.makedirs(DS_DIR,exist_ok=True)
    #
    df_ds = pd.read_csv(csv_dir)
    assert np.all(df_ds.sample_id.value_counts() == 1)
    for _,series in df_ds.iterrows():
        sample_id = series.sample_id
        file      = series.molecule_h5
        dest = os.path.join(SHARE_DIR,ds_name,sample_id)
        src_prefix = os.path.dirname(file)
        os.makedirs(dest,exist_ok=True)
        for file_mv in files_to_move:
            src  = os.path.join(src_prefix,file_mv)
            shutil.copy(src,dest)
        for file_mv in folder_to_move:
            src  = os.path.join(src_prefix,file_mv)
            dest_folder = os.path.join(dest,file_mv)
            os.makedirs(dest_folder)
            copy_tree(src,dest_folder)
