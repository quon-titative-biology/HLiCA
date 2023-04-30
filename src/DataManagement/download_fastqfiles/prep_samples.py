"""
Take in GoogleSheetMetadata and prepare the metadata for alignment pipeline
Applicable for the lattice datasets where GoogleSheetMetadata were provided
"""

import os,glob
import re

import pandas as pd


# ==============================================================================
# Specify metadata directory
# ==============================================================================

# Set METADATA_DIR
# METADATA_DIR = "fastqfiles/GoogleSheetMetadata.csv"
METADATA_DIR = "fastqfiles/Gruen/GoogleSheetMetadata.csv"
METADATA_DIR = "fastqfiles/Henderson/GoogleSheetMetadata.csv"
METADATA_DIR = "fastqfiles/Guiliams/GoogleSheetMetadata.csv"


# ==============================================================================
# Save metadata per fastq file
# ==============================================================================

df_metadata = pd.read_csv(METADATA_DIR)

# Set FILE_DIR, the relative path where file is stored
# 's3://submissions-czi004liv/andrews_2021/SRR7276476/McGilvery_Sonya__TLH_June_29_MissingLibrary_1_CBC03ANXX/SRR7276476_S1_L007_R2_001.fastq.gz
# 'submissions-czi004liv/andrews_2021/SRR7276476/McGilvery_Sonya__TLH_June_29_MissingLibrary_1_CBC03ANXX/SRR7276476_S1_L007_R2_001.fastq.gz

if 'FILE_DIR_MODIFIED' in df_metadata.columns:

    col_to_add = df_metadata.FILE_DIR_MODIFIED

else:
    pattern = "^s3://"
    col_to_add = df_metadata.s3_uri.map(lambda s: re.sub(pattern,"",s))

    # Modify FILE_DIR to
    def format_fastq_name(s):
        d2 = {"R1.fq.gz": "S1_L001_R1_001.fastq.gz",
              "R2.fq.gz": "S1_L001_R2_001.fastq.gz"}
        for k,v in d2.items():
            s = s.replace(k,v)
        return s

    col_to_add = col_to_add.map(format_fastq_name)

df_metadata.insert(1,'FILE_DIR',col_to_add)

# Get the fastqs dir, input to cellranger count
# 'submissions-czi004liv/andrews_2021/SRR7276476/McGilvery_Sonya__TLH_June_29_MissingLibrary_1_CBC03ANXX/SRR7276476_S1_L007_R2_001.fastq.gz
# 'submissions-czi004liv/andrews_2021/SRR7276476/McGilvery_Sonya__TLH_June_29_MissingLibrary_1_CBC03ANXX'
col_to_add = df_metadata.FILE_DIR.map(lambda s: os.path.dirname(s))
df_metadata.insert(1,'FASTQS_DIR',col_to_add)

# Get the sample name by stripping the naming scheme
# 'submissions-czi004liv/andrews_2021/SRR7276476/McGilvery_Sonya__TLH_June_29_MissingLibrary_1_CBC03ANXX/SRR7276476_S1_L007_R2_001.fastq.gz
# 'SRR7276476_S1_L007_R2_001.fastq.gz'
# _S\d+ = _S1
# _L\d+ =  _L007
# _\D\d+ = _R2
# _\d+ = _001
# .f[a-zA-z]*q.gz$ = .fastq.gz / .fq.gz
# 'SRR7276476'

#
pattern = "_S\d+_L\d+_\D\d+_\d+.f[a-zA-z]*q.gz$"
col_to_add = df_metadata.FILE_DIR.map(lambda s: re.sub(pattern,"",os.path.basename(s)))

#
pattern = "_S\d+_L\d+_\D\d+_\d+.fq.gz$"
col_to_add = col_to_add.map(lambda s: re.sub(pattern,"",os.path.basename(s)))

# col_to_add = df_metadata.FILE_DIR.map(lambda s: "_".join(os.path.basename(s).split("_")[:3]))
df_metadata.insert(1,'SAMPLE',col_to_add)

# Save the file
filedir = os.path.splitext(METADATA_DIR)[0]+'_file.csv'
df_metadata.to_csv(filedir,index=False)

# ==============================================================================
# Save the metadata per sample
# ==============================================================================

# Drop duplicate sample
df_metadata_sample = df_metadata.drop_duplicates(['FASTQS_DIR','SAMPLE'])
df_metadata_sample.value_counts(['SAMPLE','FASTQS_DIR'])
df_metadata.value_counts(['SAMPLE','FASTQS_DIR'])

filedir = os.path.splitext(METADATA_DIR)[0]+'_sample.csv'
df_metadata_sample.to_csv(filedir,index=False)
