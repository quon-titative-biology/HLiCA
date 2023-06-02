SHARE_DIR="alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13/share_qc"

source="/share/quonlab/workspaces/czi_liver_atlas/data/LiverNetworkData/${SHARE_DIR}"

dest="/group/gquongrp/collaborations/HLiCA/${SHARE_DIR}"
mkdir -p $dest

scp -r yonchoi@esb.genomecenter.ucdavis.edu:$source $dest
