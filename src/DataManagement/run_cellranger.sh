"""
3. Perform sequence alignment with cellranger
https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count#cr-mkfastq
"""

module load cellranger/7.0.0

#### Apply on all

ref="GRCh38p13_gencode_v42"
transcriptome="genome/GRCh38/GRCh38_gencode/v42/GRCh38p13_gencode_v42"

# Set base OUT_DIR
OUT_DIR="alignment/ref_${ref}"
mkdir -p $OUT_DIR

# Set path to dataset specific metadata file
METADATA_SAMPLE_DIR="fastqfiles/Gruen/GoogleSheetMetadata_sample.csv"
METADATA_SAMPLE_DIR="fastqfiles/Henderson/GoogleSheetMetadata_sample.csv"
METADATA_SAMPLE_DIR="fastqfiles/Guiliams/GoogleSheetMetadata_sample.csv"
METADATA_SAMPLE_DIR="fastqfiles/Dasgupta/GoogleSheetMetadata_sample.csv"

dos2unix $METADATA_SAMPLE_DIR

mem="128" #GB
time="2" #hours
ntasks="1"
cpus_per_task="16"
partition="production"

sed 1d $METADATA_SAMPLE_DIR |
while IFS=, read -r s3_uri sample FASTQS_DIR FILE_DIR rest; do

  echo "Starting alignment: sample info:"
  echo $sample
  echo $FASTQS_DIR

  OUT_DIR_SAMPLE="$OUT_DIR/$FASTQS_DIR"
  mkdir -p $OUT_DIR_SAMPLE

  if [ -f "${OUT_DIR_SAMPLE}/${sample}/outs/filtered_feature_bc_matrix.h5" ]

  then
    echo 'matrix exists'

  else
    echo 'matrix does not exist, align'

    ls "fastqfiles/${FASTQS_DIR}"

    time cellranger count --id="${sample}" \
                          --transcriptome=${transcriptome} \
                          --fastqs="fastqfiles/${FASTQS_DIR}" \
                          --localcores=32 \
                          --localmem=256 \
                          --sample=${sample} \
                          --no-bam

    # Run through job manager
    # sbatch --partition=${partition} \
    #        --gres=gpu:1 \
    #        --time="${time}:00:00" \
    #        --mem="${mem}G" \
    #        --error="${OUT_DIR_SAMPLE}/slurm_msg/error_FA.%j.out" \
    #        --output="${OUT_DIR_SAMPLE}/slurm_msg/output_FA.%j.out" \
    #        --cpus-per-task="${cpus_per_task}" \
    #        --ntasks="${ntasks}" \
    #        alignment/cellranger.sh $sample $transcriptome "fastqfiles/${FASTQS_DIR}" 8 256 ${sample} ${OUT_DIR_SAMPLE}

    mv $sample $OUT_DIR_SAMPLE

  fi

  echo $OUT_DIR_SAMPLE


  echo "Alignment finished"
  echo ""

  sleep 0.1

done
