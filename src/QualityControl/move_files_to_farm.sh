####

sample_name="C70_Caudate_3pr_v2"
sample_name="Human_Caudate_C72_3pr_v2"
sample_name="TLH_18_09_20_CITE_3"
sample_prefix="submissions-czi004liv/andrews_2021/PRJNA769141/${sample_name}"

src="yonchoi@esb.genomecenter.ucdavis.edu:/share/quonlab/workspaces/czi_liver_atlas/data/LiverNetworkData/alignment/ref_GRCh38p13_gencode_v42/${sample_prefix}/outs/raw_feature_bc_matrix"
dest="/group/gquongrp/collaborations/HLiCA/alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13/share_qc/share_qc/Toronto/${sample_name}"

mkdir -p $dest

scp -r $src $dest

#
src="yonchoi@esb.genomecenter.ucdavis.edu:/share/quonlab/workspaces/czi_liver_atlas/data/LiverNetworkData/alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13/share_qc"
dest="/group/gquongrp/collaborations/HLiCA/alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13/share_qc3"

mkdir -p $dest

scp -r $src $dest
