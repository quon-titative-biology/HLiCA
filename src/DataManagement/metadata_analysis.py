"""
Perform basic metadata anaylsis
"""

import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
import seaborn as sns

df_metadata = pd.read_csv("Liver Network data - Andrews_2021_hepatol_commun.csv")

# Get name of columns that are categorical
col_barplots = [col for col in df_metadata.columns if len(df_metadata[col].value_counts()) < 20]
len(col_barplots)

df_plot = df_metadata[col_barplots].melt(var_name='Metadata',value_name='Value')

g = sns.FacetGrid(data=df_plot, col='Metadata', col_wrap=4,
                  sharey=False, sharex=False)
g.map_dataframe(sns.countplot, x='Value')

for axes in g.axes.flat:
    max_length = 10
    shortened_label = [item.get_text()[:max_length] for item in axes.get_xticklabels()]
    _ = axes.set_xticklabels(shortened_label, rotation=90)

plt.tight_layout()

plt.savefig('metadata_distribution.pdf',
            bbox_inches='tight')
plt.close()
