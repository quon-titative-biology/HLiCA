import os,glob

import pandas as pd

ALIGNMENT_DIR = 'alignment/ref_GRCh38p13_gencode_v42'

sample_meta_dir = "fastqfiles/Toronto/GoogleSheetMetadata_sample.csv"

df_meta = pd.read_csv(sample_meta_dir)

os.listdir(ALIGNMENT_DIR)
