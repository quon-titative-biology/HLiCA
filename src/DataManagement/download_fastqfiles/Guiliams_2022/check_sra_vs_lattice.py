"""
Check the download files from SRA match with metadata from lattice
"""

import os,glob

import pandas as pd

DATA_DIR = 'fastqfiles/sra_downloads'
METADATA_DIR = 'fastqfiles/Guiliams/GoogleSheetMetadata.csv'

df_metadata = pd.read_csv(METADATA_DIR,skiprows=1)

SRAfiles = glob.glob(os.path.join(DATA_DIR,'*.fastq.gz'))
df_sra = pd.DataFrame({'FILE_DIR': SRAfiles})

def parse_sra(sra):
    filename = os.path.basename(sra)
    SRAname = filename.replace(".fastq.gz",'')
    SRA_id,file_type = SRAname.split('_')
    return SRA_id, file_type
