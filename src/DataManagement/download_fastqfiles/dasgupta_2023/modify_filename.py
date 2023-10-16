import os, glob
import re

import pandas as pd

METADATA_DIR = "fastqfiles/DasGupta/GoogleSheetMetadata.csv"

df_metadata = pd.read_csv(METADATA_DIR)

df_metadata['SAMPLE'] = df_metadata.library_alias.str.replace('ramanuj-dasgupta:','').str.replace('_lib1','')
df_metadata['SAMPLE_FOLDER'] = df_metadata.s3_uri.map(lambda s: s.split("/")[-2])

df_metadata.groupby(['library_alias','SAMPLE','tissue']).size()
df_metadata.groupby(['SAMPLE','SAMPLE_FOLDER','tissue']).size()
df_metadata.groupby(['SAMPLE_FOLDER','SAMPLE','tissue']).size()

def add_to_filename(DIR,str_add):
    path,ext = os.path.splitext(DIR)
    return path+str_add+ext

df_metadata.to_csv(add_to_filename(METADATA_DIR,"_orig"),index=False)

df_metadata = df_metadata[df_metadata.tissue == 'liver']

df_metadata.groupby(['library_alias','SAMPLE','tissue']).size()
df_metadata.groupby(['SAMPLE','SAMPLE_FOLDER','tissue']).size()
df_metadata.groupby(['SAMPLE_FOLDER','SAMPLE','tissue']).size()

df_metadata.to_csv(METADATA_DIR,index=False)

def reformat(row):
    s = row.s3_uri
    sample = row.SAMPLE
    s = s.replace("_unaligned_","_")
    s = s.replace("_HJVHFDSXY_","_")
    head,end = s.split("_L0")
    # library = head.split('/')[-2]
    file = sample + "_S1_L0" + end
    return os.path.join(os.path.dirname(os.path.dirname(head)), sample, file)

df_metadata['FILE_DIR_MODIFIED'] = df_metadata.apply(reformat, axis=1)

def replace_L(s,i):
    L = re.search("_L\d+_",s)[0]
    L_new = list(L)
    L_new[3] = str(i)
    L_new = "".join(L_new)
    return s.replace(L,L_new)

def rename_duplicates(df_metadata):
    counts = df_metadata['FILE_DIR_MODIFIED'].value_counts()
    duplicates = counts[counts > 1]
    files = df_metadata['FILE_DIR_MODIFIED'].values
    for d in duplicates.sort_index().index:
        df_meta_dup = df_metadata[df_metadata.FILE_DIR_MODIFIED == d]
        new_files = [replace_L(s,i) for i,s in enumerate(df_meta_dup['FILE_DIR_MODIFIED'])]
        files[df_meta_dup.index] = new_files
        print(df_meta_dup.index)
    df_metadata['FILE_DIR_MODIFIED'] = files
    return df_metadata

df_metadata = rename_duplicates(df_metadata)

assert df_metadata.FILE_DIR_MODIFIED.value_counts().max() == 1, "Files with same name exists"
df_metadata.to_csv(METADATA_DIR,index=False)

#### After download
df_metadata = pd.read_csv(METADATA_DIR)
import numpy as np
assert np.all(df_metadata.tissue == 'liver')

# Move
for _,row in df_metadata.iterrows():
    src   = row.s3_uri.replace('s3://','fastqfiles/')
    dest  = row.FILE_DIR_MODIFIED.replace('s3://','fastqfiles/')
    if os.path.isfile(src) and not os.path.isfile(dest):
        print('Moving')
        os.rename(src,dest)
    else:
        print('Not moving')
        print(os.path.isfile(src), not os.path.isfile(dest))
        print(src, dest)

assert np.all([os.path.isfile(dest) for dest in df_metadata.FILE_DIR_MODIFIED.str.replace('s3://','fastqfiles/')])

# # Move
# for _,row in df_metadata.iterrows():
#     src   = row.s3_uri.replace('s3://','fastqfiles/')
#     dest  = row.FILE_DIR_MODIFIED.replace('s3://','fastqfiles/')
#     print(os.path.isfile(src), not os.path.isfile(dest))
#     # Move
#     if os.path.isfile(src) and not os.path.isfile(dest):
#         print('Moving')
#         os.rename(src,dest)
#     # # Reverse
#     # if os.path.isfile(dest):
#     #     print('Moving')
#     #     os.rename(dest,src)

df_metadata.to_csv(METADATA_DIR,index=False)
