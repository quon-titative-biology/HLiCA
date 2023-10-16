"""
User to manually combine h5 files from cellranger count into a single h5ad file
"""

import os,glob
import pandas as pd
import scanpy as sc

import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--METADATA_DIR', type=str, default='',
                    help='Directory where GoogleSheetMetadata_sample.csv')
parser.add_argument('--dataset_name', type=str, default='',
                    help='Name of the dataset')
parser.add_argument('--ALIGNMENT_DIR', type=str, default='',
                    help='Directory where alignment outputs are stored')
parser.add_argument('--OUT_DIR', type=str, default='',
                    help='Directory where the combined read in h5 form will be stored')



# ==============================================================================
# Set directory
# ==============================================================================

METADATA_DIR = "fastqfiles/Gruen/GoogleSheetMetadata_sample.csv"
dataset_name = 'Gruen'
METADATA_DIR = "fastqfiles/Henderson/GoogleSheetMetadata_sample.csv"
dataset_name = 'Henderson'
METADATA_DIR = "fastqfiles/Guiliams/GoogleSheetMetadata_sample.csv"
dataset_name = 'GuilliamsScott'
METADATA_DIR = "fastqfiles/Dasgupta/GoogleSheetMetadata_sample.csv"
dataset_name = 'DasGupta'

ALIGNMENT_DIR = "alignment/ref_GRCh38p13_gencode_v42"
OUT_DIR = os.path.join(ALIGNMENT_DIR,"share_2023_05_13")
os.makedirs(OUT_DIR,exist_ok=True)

# ==============================================================================
#
# ==============================================================================

df_meta = pd.read_csv(METADATA_DIR)

# Check for h5ad

df_meta['result_dir'] = 'NA'

adata_list = []

# Load the data
for i,sample in df_meta.iterrows():
    result_dir = os.path.join(ALIGNMENT_DIR,
                            sample.FASTQS_DIR,
                            sample.SAMPLE)
    sample.result_dir = result_dir
    #
    try:
        os.listdir(os.path.join(result_dir,'outs'))
    except:
        print('doesnt exist')
    #
    if os.path.isfile(os.path.join(result_dir,'outs','raw_feature_bc_matrix.h5')):
        adata = sc.read_10x_h5(os.path.join(result_dir,'outs','raw_feature_bc_matrix.h5'))
        adata_filtered = sc.read_10x_h5(os.path.join(result_dir,'outs','filtered_feature_bc_matrix.h5'))
        adata.obs['filtered'] = False
        adata.obs['filtered'][adata_filtered.obs.index] = True
        adata.obs['Dataset'] = dataset_name
        # #
        # for k,v in sample.items():
        #     adata.obs[k] = v
        #
        adata_list.append(adata)
    else:
        pass

# Make names unique
for adata in adata_list:
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

# Combine into one anndata to save as h5ad
adata_combined = sc.concat(adata_list, axis=0)
adata_combined.write_h5ad(os.path.join(OUT_DIR,f"{dataset_name}.h5"))
