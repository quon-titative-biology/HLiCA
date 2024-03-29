"""
Analyze the metrics summary from cell ranger count
"""

import os,glob

import pandas as pd

from matplotlib import pyplot as plt
import seaborn as sns

ALIGNMENT_DIR = "alignments"
ALIGNMENT_DIR = "alignment/ref_GRCh38p13_gencode_v42"
ALIGNMENT_DIR = "fastqfiles/submissions-czi004liv/gruen_2023"

PLOT_DIR = ALIGNMENT_DIR

RESULTS_DIR = glob.glob(os.path.join(ALIGNMENT_DIR,"*counts_and_metrics"))
RESULTS_DIR = glob.glob(os.path.join(ALIGNMENT_DIR,"*/outs"))
RESULTS_DIR = [os.path.dirname(g) for g in RESULTS_DIR]

def load_metric(RESULT_DIR):
    sample_name = os.path.basename(RESULT_DIR)
    def parse_sample_name(sample_name):
        sample,ref = sample_name.split("-")
        return sample,ref
    sample,ref = parse_sample_name(sample_name)
    try:
        df_metric = pd.read_csv(os.path.join(RESULT_DIR,'outs/metrics_summary.csv'))
    except:
        df_metric = pd.read_csv(os.path.join(RESULT_DIR,'metrics_summary.csv'))
    df_metric['sample'] = sample
    df_metric['ref'] = ref
    return df_metric

def load_metric(RESULT_DIR):
    sample_name = os.path.basename(RESULT_DIR)
    try:
        df_metric = pd.read_csv(os.path.join(RESULT_DIR,'outs/metrics_summary.csv'))
    except:
        df_metric = pd.read_csv(os.path.join(RESULT_DIR,'metrics_summary.csv'))
    df_metric['sample'] = sample_name
    return df_metric

df_metric = pd.concat([load_metric(r) for r in RESULTS_DIR])
df_metric['ref'] = 'GRCh38p13_gencode_v42'

df_meta = pd.read_csv('fastqfiles/Gruen/GoogleSheetMetadata_sample.csv')

metric_to_plot = []
df_plot = df_metric.copy()
df_plot = df_plot.melt(id_vars = ['sample','ref'], value_name='Value', var_name='Metric')
df_plot.Value = df_plot.Value.apply(lambda s: float(str(s).replace(",","").replace("%","")))

# g = sns.FacetGrid(data=df_plot, columns='Metric')
g = sns.catplot(data=df_plot,
            col='Metric', col_wrap=5,
            kind="bar",
            x='sample',y='Value',hue='ref',
            sharey=False,sharex=False)
g.set_xticklabels(rotation=90)
plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR,"alignments_metric.pdf"),
            bbox_inches='tight')
plt.close()

metric_to_plot = ['Reads Mapped to Genome',
                  'Reads Mapped Confidently to Genome',
                  'Reads Mapped Confidently to Transcriptome',
                  'Reads Mapped Confidently to Intergenic Regions',
                  'Reads Mapped Confidently to Intronic Regions',
                  'Reads Mapped Confidently to Exonic Regions',
                  'Reads Mapped Antisense to Gene'] \
                  + ['sample','ref']
df_plot = df_metric[metric_to_plot]
df_plot = df_plot.melt(id_vars = ['sample','ref'], value_name='Value', var_name='Metric')
df_plot.Value = df_plot.Value.apply(lambda s: float(str(s).replace(",","").replace("%","")))

g = sns.catplot(data=df_plot,
            col='Metric',
            kind="bar",
            x='sample',y='Value',hue='ref',
            sharey=False,sharex=False)
g.set_xticklabels(rotation=90)
plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR,"alignments_metric-ReadsMapped.pdf"),
            bbox_inches='tight')
plt.close()
