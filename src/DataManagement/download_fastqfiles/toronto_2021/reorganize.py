"""
Re-organize fastqfiles so all the samples are in the same folder
* Run after downloading

/andrew_2021
    /PRJNA769141
        /SRR#
        /SRR#
    /SRR#...
    /SRR#...

SRR#s under /PRJNA769141 should be grouped by their sample name.
"""

import os,glob

import numpy as np
import pandas as pd

fastqfiles="fastqfiles" # relative
fastqfiles="/share/quonlab/workspaces/czi_liver_atlas/data/LiverNetworkData/fastqfiles" #absolute

GroupName=f"{fastqfiles}/Toronto"
GoogleSheetFile=f"{GroupName}/GoogleSheetMetadata.csv"

df_meta = pd.read_csv(GoogleSheetFile)

df_meta['s3_uri_orig'] = df_meta['s3_uri']


# Some files are contained under additional folder /PRJNA769141
# Remove so all SRR# level directories are removed at andrew_2021/PRJNA769141

import re
s='s3://submissions-czi004liv/andrews_2021/SRR7276476/McGilvery_Sonya__TLH_June_29_MissingLibrary_1_CBC03ANXX/SRR7276476_S1_L007_R2_001.fastq.gz'
s = df_meta.s3_uri[50]
pattern = '/PRJNA769141/SRR\d+/'

df_meta.s3_uri = df_meta.s3_uri.map(lambda s: re.sub(pattern,'/PRJNA769141/',s))

df_meta.to_csv(GoogleSheetFile,index=False)

def convert2file(f):
    return os.path.join(fastqfiles,f.replace('s3://',''))

for i,file in df_meta.iterrows():
    print(i)
    if file.s3_uri != file.s3_uri_orig:
        os.rename(convert2file(file.s3_uri_orig),
                  convert2file(file.s3_uri))
