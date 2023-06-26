import os,glob

import pandas as pd

filenames = glob.glob(os.path.join('alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13/share_qc/share_qc',"*",'stats.csv'))

df_stat_list = []
for file in filenames:
    df_stat = pd.read_csv(file)
    df_stat = df_stat.mean(1)
    dataset = os.path.basename(os.path.dirname(file))
    df_stat['dataset'] = dataset
    df_stat_list.append(df_stat)

pd.concat(df_stat_list,axis=1)

df_stat_list = []
for file in filenames:
    df_stat = pd.read_csv(file)
    df_stat = df_stat.reset_index().melt(id_vars='index')
    dataset = os.path.basename(os.path.dirname(file))
    df_stat['dataset'] = dataset
    df_stat_list.append(df_stat)

df_stat = pd.concat(df_stat_list)

import seaborn as sns
import matplotlib.pyplot as plt

g = sns.catplot(data=df_stat,col='index',col_wrap=5,
            x='dataset',y='value',kind='violin',sharey=False)
plt.savefig('alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13/share_qc/share_qc/stat_violin_plot.pdf')
plt.close()
