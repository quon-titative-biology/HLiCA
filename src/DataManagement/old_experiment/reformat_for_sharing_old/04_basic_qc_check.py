import os,glob

import pandas as pd

import scanpy as sc

from matplotlib import pyplot as plt
import seaborn as sns

SHARE_DIR = 'alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13/share'

adata_files = glob.glob(os.path.join(SHARE_DIR,"*/filtered_feature_bc_matrix.h5"))

for adata_file in adata_files:
    adata = sc.read_10x_h5(adata_file)
    # Basic filtering
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    # Mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    #
    # sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    #          jitter=0.4, multi_panel=True)
    # plt.savefig(os.path.join(os.path.dirname(adata_file),'basic_qc.pdf'))
    # plt.close()
    #
    g = sns.violinplot(data=adata.obs, y='pct_counts_mt')
    plt.savefig(os.path.join(os.path.dirname(adata_file),'basic_qc_pct_counts_mt.pdf'))
    plt.close()


####

dataset_names = os.listdir(SHARE_DIR)

df_total = []

for ds_name in dataset_names:
    DS_DIR = os.path.join(SHARE_DIR, ds_name)
    h5_files = glob.glob(os.path.join(DS_DIR,"*/filtered_feature_bc_matrix.h5",))
    samples = [f.split('/')[-2] for f in h5_files]
    sizes = []
    for adata_file in h5_files:
        adata = sc.read_10x_h5(adata_file)
        sizes.append(adata.shape[0])
    df_meta = pd.DataFrame({'Sample': samples,
                           'n_cells': sizes,
                           'dataset': ds_name})
    df_total.append(df_meta)

df_total = pd.concat(df_total)

# of samples
df_total.dataset.value_counts()
# of cells
df_total.groupby('dataset').sum()

g = sns.violinplot(data=df_total, x='dataset', y='n_cells', cut=0)
plt.savefig(os.path.join(SHARE_DIR,'n_cells.pdf'))
plt.close()
