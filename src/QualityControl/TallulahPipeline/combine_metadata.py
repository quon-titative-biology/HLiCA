import os, glob
import data

# Might need to modify to grab all the {dataset}_metadata.csv
df_dir_list = glob.glob('*.csv')
df_dir_list = [df_dir for df_dir in df_dir_list if not df_dir.startswith('Combined')]

df_list = []

for df_dir in df_dir_list:
    df_ = pd.read_csv(df_dir)
    df_['Dataset'] = df_dir.split("_")[0]
    df_list.append(df_)

df = pd.concat(df_list)
df = df.iloc[:,1:]
df = df.iloc[:,1:]
col_to_move = df.pop('Dataset')
df.insert(0, "Dataset", col_to_move)

df.to_csv('Combined_metadata.csv')

df_qc = df[df.QC_Pipeline]
df_qc = df_qc.fillna('NaN')

df_qc.to_csv('Combined_metadata_qconly.csv')

df_qc.groupby(['Dataset','assay']).size()
df_qc.Dataset.value_counts()
df_qc.assay.value_counts()
