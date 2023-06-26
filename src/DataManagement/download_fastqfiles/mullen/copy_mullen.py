import os,glob
import shutil
from distutils.dir_util import copy_tree

DATA_DIR = 'alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13/share/Mullen'
OUT_DIR = 'alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13/share_qc/Mullen'
os.makedirs(OUT_DIR,exist_ok=True)

FILES_DIR = glob.glob(os.path.join(DATA_DIR,"*",'outs/raw_feature_bc_matrix'))

for file in FILES_DIR:
    outs_dir = os.path.dirname(file)
    sample   = os.path.basename(os.path.dirname(outs_dir))
    SAMPLE_DIR_OUT = os.path.join(OUT_DIR,sample,'raw_feature_bc_matrix')
    os.makedirs(SAMPLE_DIR_OUT,exist_ok=True)
    copy_tree(file,SAMPLE_DIR_OUT)
    # outs_10x_dir = os.path.dirname(outs_dir)
