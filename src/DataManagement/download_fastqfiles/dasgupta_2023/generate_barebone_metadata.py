import os,glob

import pandas as pd

FASTQ_DIR = 'fastqfiles/submissions-czi004liv/dasgupta_2023'


# Search for all fastqfiles
fastqfiles = glob.glob(os.path.join(FASTQ_DIR,"**","*.fastq.gz"),recursive=True)

# Initialize metadata data frame
df_meta = pd.DataFrame({"s3_uri":fastqfiles})

# Remove s3:// from the filedir to get s3_uri
df_meta.s3_uri = df_meta.s3_uri.apply(lambda s: 's3://'+s.replace('fastqfiles/',''))

df_meta['s3_uri_orig'] = df_meta.s3_uri
df_meta['filename_orig'] = df_meta.s3_uri.map(lambda s: os.path.basename(s))

# Format the filenames to match 10x naming scheme

def reformat(s):
    if '_unaligned' in s:
        s = s.replace('_unaligned',"")
        prefix = '_unaligned_S1'
        s1,s2 = s.split('_L')
        s = s1 + prefix + '_L' + s2
    if '_HJVHFDSXY' in s:
        s = s.replace('_HJVHFDSXY',"")
        prefix = '_HJVHFDSXY'
        s1,s2 = s.split('_L')
        s = s1 + prefix + '_L' + s2
    return s

df_meta['filename'] = df_meta.filename_orig.apply(reformat)
s3_uri_new = df_meta.s3_uri.map(lambda s: os.path.dirname(s)) + "/" + df_meta.filename
df_meta.s3_uri = s3_uri_new

# Rename to new format
for s3_new, s3_old in zip(df_meta.s3_uri, df_meta.s3_uri_orig):
    old_dir = 'fastqfiles/'+s3_old.replace('s3://','')
    new_dir = 'fastqfiles/'+s3_new.replace('s3://','')
    os.path.isfile(old_dir)
    os.path.isfile(new_dir)
    os.rename(old_dir,
              new_dir)

os.makedirs('fastqfiles/Dasgupta',exist_ok=True)
df_meta.to_csv('fastqfiles/Dasgupta/GoogleSheetMetadata.csv',index=False)
