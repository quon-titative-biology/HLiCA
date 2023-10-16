import os
import pandas as pd

DATA_DIR = 'alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13/share/Mullen'
OUT_DIR  = 'fastqfiles/Mullen'
os.makedirs(OUT_DIR,exist_ok=True)

df_meta_in = pd.read_csv(os.path.join(DATA_DIR,'metadata_for_analysis.txt'),sep='\t')

TEMPLATE_META_DIR = 'fastqfiles/Gruen/GoogleSheetMetadata_sample.csv'
df_meta_template = pd.read_csv(TEMPLATE_META_DIR)
columns = df_meta_template.columns

df_meta = pd.DataFrame(index=df_meta_in.index,
                       columns=columns)

df_meta.donor_sex = df_meta_in.Sex
df_meta.SAMPLE = df_meta_in.SampleID.str.replace('-','_')
df_meta.tissue = df_meta_in.SampleName.apply(lambda s: s.split('-')[0].lower())
df_meta.donor_age = df_meta_in.Age

df_meta.to_csv(os.path.join(OUT_DIR,'GoogleSheetMetadata_sample.csv'),index=None)
