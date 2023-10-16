import os,glob

import pandas as pd

from matplotlib import pyplot as plt
import seaborn as sns

BASE_DIR = 'alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13/share_qc'

filenames = glob.glob(os.path.join('alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13/share_qc/share_qc',"*",'stats.csv'))
filenames = glob.glob(os.path.join(BASE_DIR,"*",'stats.csv'))

def analyze_qc(DS_DIR):
    FILES = glob.glob(os.path.join(DS_DIR,"*",'raw_feature_bc_matrix/qc_output/Cleaned_output_status.csv'))
    df_status_all = []
    for f in FILES:
        sample    = f.split('/')[-4]
        df_status = pd.read_csv(f)
        df_status['SAMPLE'] = sample
        df_status_all.append(df_status)
    df_status = pd.concat(df_status_all)
    df_status['status'] = 'NaN'
    for i in ['soupx','emptydrops']:
        if i in df_status.columns:
            df_status.status[~df_status[i].isna()] = i+'-'+ df_status[~df_status[i].isna()][i]
    return df_status

#### Create summarized table where each row is dataset
df_stat_list = []
for file in filenames:
    df_stat = pd.read_csv(file)
    n_sample = df_stat.shape[1]
    df_stat = df_stat.mean(1)
    dataset = os.path.basename(os.path.dirname(file))
    df_stat['dataset'] = dataset
    df_stat['n_sample'] = n_sample
    df_stat_list.append(df_stat)
    df_status = analyze_qc(os.path.dirname(file))
    df_status.to_csv(os.path.join(os.path.dirname(file),'status.csv'))

df_stat = pd.concat(df_stat_list,axis=1)
df_stat.columns= df_stat.loc['dataset']
df_stat = df_stat.transpose()
df_stat['Total cell empty drop'] = (df_stat['EmptyDrop # cell'] * df_stat['n_sample']).astype('int')
df_stat = df_stat.transpose()

df_plot = df_stat.reset_index().melt(id_vars='index',var_name='dataset')
df_plot = df_plot.rename(columns={'index':'Stat'})
df_plot = df_plot[df_plot.value.map(lambda v: not isinstance(v,str))]

g = sns.catplot(data=df_plot,col='Stat',col_wrap=5,
                x='dataset',y='value',kind='bar', sharex=False, hue='dataset',
                sharey=False,dodge=False)
g.set_xticklabels(rotation=90)
plt.tight_layout()
plt.savefig(os.path.join(BASE_DIR,'stat_bar_plot.pdf'))
plt.close()

#### Create table where row is sample for plotting
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
plt.savefig(os.path.join(BASE_DIR,'stat_violin_plot.pdf'))
plt.close()
