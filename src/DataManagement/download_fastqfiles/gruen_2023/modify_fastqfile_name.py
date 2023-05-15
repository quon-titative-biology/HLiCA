"""
Check the download files from SRA match with metadata from lattice
"""

import os,glob
import re

import pandas as pd

FASTQ_DIR    = 'fastqfiles'
METADATA_DIR = 'fastqfiles/Gruen/GoogleSheetMetadata.csv'

df_metadata = pd.read_csv(METADATA_DIR)

def format_fastq_name(s):
    d2 = {"R1.fq.gz": "S1_L001_R1_001.fastq.gz",
          "R2.fq.gz": "S1_L001_R2_001.fastq.gz",
          "R1.fastq.gz": "S1_L001_R1_001.fastq.gz",
          "R2.fastq.gz": "S1_L001_R2_001.fastq.gz"}
    for k,v in d2.items():
        s = s.replace(k,v)
    return s

pattern = "^s3://"
FILE_DIR = df_metadata.s3_uri.map(lambda s: re.sub(pattern,"",s))

FILE_DIR_NEW = FILE_DIR.map(format_fastq_name)

df_metadata['FILE_DIR_MODIFIED'] = FILE_DIR_NEW
df_metadata.to_csv(METADATA_DIR, index=False)

# Rename
for file_old, file_new in zip(FILE_DIR,FILE_DIR_NEW):
    file_old_exist = os.path.isfile(os.path.join(FASTQ_DIR,file_old))
    file_new_exist = os.path.isfile(os.path.join(FASTQ_DIR,file_new))
    assert file_old_exist or file_new_exist, f"Neither new or old file exist for {file_old}"
    if file_old_exist and not file_new_exist:
        os.rename(os.path.join(FASTQ_DIR,file_old),
                  os.path.join(FASTQ_DIR,file_new))
