import os
import pandas as pd

ref="GRCh38p13_gencode_v42"
OUT_DIR_BASE=f"alignment/ref_{ref}"

#### Add column to see if sample passed through cellranger pipeline
for ds in ['Gruen','Henderson','Toronto','Guilliams','DasGupta']:
    print(ds)
    df_meta = pd.read_csv(os.path.join('fastqfiles',ds,'GoogleSheetMetadata_sample.csv'))
    CellRanger = []
    for i,row in df_meta.iterrows():
        OUT_DIR = os.path.join(OUT_DIR_BASE,row.FASTQS_DIR,
                               row.SAMPLE,
                               'outs/filtered_feature_bc_matrix.h5')
        if os.path.exists(OUT_DIR):
            CellRanger.append(True)
        else:
            print(OUT_DIR)
            print(row.SAMPLE)
            CellRanger.append(False)
    df_meta['AlignmentCellRanger'] = CellRanger
    df_meta.to_csv(os.path.join('fastqfiles',ds,'GoogleSheetMetadata_sample_aligned.csv'),index=None)


ds = 'Mullen'
df_meta = pd.read_csv(os.path.join('fastqfiles',ds,'GoogleSheetMetadata_sample.csv'))
df_meta['AlignmentCellRanger'] = True
df_meta.to_csv(os.path.join('fastqfiles',ds,'GoogleSheetMetadata_sample_aligned.csv'),index=None)

#### Add column to see if sample passed through QC pipeline
datasets = ['Gruen','Henderson','Toronto','GuilliamsScott','DasGupta','Mullen']
OUT_DIR_BASE = 'alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13/share_qc'
# OUT_DIR_BASE = 'alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13/share_qc/Mullen/6854_24/raw_feature_bc_matrix/qc_output/Cleaned_output_SoupX.rds'
# OUT_DIR_BASE = 'alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13/share_qc/Mullen/6854-24/raw_feature_bc_matrix/qc_output/Cleaned_output_SoupX.rds'

datasets = ['DasGupta']
for ds in datasets:
    print(ds)
    df_meta = pd.read_csv(os.path.join('fastqfiles',ds,'GoogleSheetMetadata_sample_aligned.csv'))
    passed = []
    for i,row in df_meta.iterrows():
        OUT_DIR = os.path.join(OUT_DIR_BASE,
                               ds,
                               row.SAMPLE,
                               'raw_feature_bc_matrix/qc_output/Cleaned_output_SoupX.rds'
                               )
        if os.path.exists(OUT_DIR):
            passed.append(True)
        else:
            print(OUT_DIR)
            print(row.SAMPLE)
            passed.append(False)
    df_meta['QC_Pipeline'] = passed
    df_meta.to_csv(os.path.join('fastqfiles',ds,'GoogleSheetMetadata_sample_qced.csv'),index=None)
