module load cellranger/7.0.0

# dataset_name="GuilliamsScott"
dataset_name="Toronto"
# dataset_name="Henderson"
# dataset_name="DasGupta"
# dataset_name="Gruen"

SHARE_DIR="alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13"

csv="$SHARE_DIR/${dataset_name}_aggr.csv"

cellranger aggr --id=${dataset_name} \
                --csv=${csv}

mv $dataset_name $SHARE_DIR

SHARE_FOLDER="${SHARE_DIR}/share"
mkdir -p $SHARE_FOLDER

dataset_list=("GuilliamsScott" "Henderson" "DasGupta" "Gruen" "Toronto")

# Copy to share folder
for dataset_name in "${dataset_list[@]}"
do

  dest="$SHARE_FOLDER/$dataset_name"
  mkdir $dest
  echo $dest

  src="$SHARE_DIR/$dataset_name/outs/count/filtered_feature_bc_matrix.h5"
  cp $src $dest

  src="$SHARE_DIR/$dataset_name.csv"
  cp $src $dest

done
