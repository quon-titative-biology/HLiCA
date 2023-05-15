"""
Create csv format for combining results of multiple experiments

Reference
https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate
https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/aggr-outputs#overview
"""

"""
Update metadata
"""

import os,glob
import numpy as np
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

# METADATA_DIR = "fastqfiles/Gruen/GoogleSheetMetadata_sample.csv"
# dataset_name = 'Gruen'
# METADATA_DIR = "fastqfiles/Henderson/GoogleSheetMetadata_sample.csv"
# dataset_name = 'Henderson'
# METADATA_DIR = "fastqfiles/Guiliams/GoogleSheetMetadata_sample.csv"
# dataset_name = 'GuilliamsScott'
METADATA_DIR = "fastqfiles/Dasgupta/GoogleSheetMetadata_sample.csv"
dataset_name = 'DasGupta'

ALIGNMENT_DIR = "alignment/ref_GRCh38p13_gencode_v42"
OUT_DIR = os.path.join(ALIGNMENT_DIR,"share_2023_05_13")
os.makedirs(OUT_DIR,exist_ok=True)

# ==============================================================================
# Save CSV for cellranger aggr
# ==============================================================================

df_aggr = pd.read_csv(METADATA_DIR)

df_aggr['molecule_h5'] = ALIGNMENT_DIR + "/" + df_aggr.FASTQS_DIR + '/' + df_aggr.SAMPLE + '/outs/molecule_info.h5'
df_aggr['sample_id'] = df_aggr.SAMPLE

df_aggr = df_aggr[df_aggr.molecule_h5.apply(os.path.isfile)]

assert np.all(df_aggr.molecule_h5.apply(os.path.isfile)), 'Not all specified molecule_h5.info does not exist'

df_aggr.to_csv(os.path.join(OUT_DIR,f'{dataset_name}.csv'),index=False)
