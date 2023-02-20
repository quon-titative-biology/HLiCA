"""
Take in
"""

import os,glob
import re

import pandas as pd

METADATA_DIR = "fastqfiles/GoogleSheetMetadata.csv"

OUT_DIR = ""
pattern = "_S*_L*_R*_*.fastq.gz"
pattern = "_S\d+_L\d+_R\d+_\d+.fastq.gz"

string = df_metadata.s3_uri[0]
re.sub(pattern,"",string)
re.sub(pattern,"",string)

df_metadata = pd.read_csv(METADATA_DIR)

# ==============================================================================
# Save metadata per fastq file
# ==============================================================================

# Set FILE_DIR, the relative path where file is stored
# 's3://submissions-czi004liv/andrews_2021/SRR7276476/McGilvery_Sonya__TLH_June_29_MissingLibrary_1_CBC03ANXX/SRR7276476_S1_L007_R2_001.fastq.gz
# 'submissions-czi004liv/andrews_2021/SRR7276476/McGilvery_Sonya__TLH_June_29_MissingLibrary_1_CBC03ANXX/SRR7276476_S1_L007_R2_001.fastq.gz
pattern = "^s3://"
col_to_add = df_metadata.s3_uri.map(lambda s: re.sub(pattern,"",s))
df_metadata.insert(1,'FILE_DIR',col_to_add)

# Get the fastqs dir, input to cellranger count
# 'submissions-czi004liv/andrews_2021/SRR7276476/McGilvery_Sonya__TLH_June_29_MissingLibrary_1_CBC03ANXX/SRR7276476_S1_L007_R2_001.fastq.gz
# 'submissions-czi004liv/andrews_2021/SRR7276476/McGilvery_Sonya__TLH_June_29_MissingLibrary_1_CBC03ANXX'
col_to_add = df_metadata.FILE_DIR.map(lambda s: os.path.dirname(s))
df_metadata.insert(1,'FASTQS_DIR',col_to_add)

# Get the sample name by stripping the naming scheme
# 'submissions-czi004liv/andrews_2021/SRR7276476/McGilvery_Sonya__TLH_June_29_MissingLibrary_1_CBC03ANXX/SRR7276476_S1_L007_R2_001.fastq.gz
# 'SRR7276476'
pattern = "_S\d+_L\d+_\D\d+_\d+.fastq.gz"
col_to_add = df_metadata.FILE_DIR.map(lambda s: re.sub(pattern,"",os.path.basename(s)))
df_metadata.insert(1,'SAMPLE',col_to_add)

# Save the file
filedir = os.path.splitext(METADATA_DIR)[0]+'_file.csv'
df_metadata.to_csv(filedir,index=False)

# ==============================================================================
# Save the metadata per sample
# ==============================================================================

# Drop duplicate sample
df_metadata_sample = df_metadata.drop_duplicates(['FASTQS_DIR','SAMPLE'])

filedir = os.path.splitext(METADATA_DIR)[0]+'_sample.csv'
df_metadata_sample.to_csv(filedir,index=False)
