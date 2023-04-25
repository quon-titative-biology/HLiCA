#!/usr/bin/env bash

echo " "
echo "Running cellranger.sh with params:"
echo "$@"

# aklog
hostname

# Load necessary modules
module load cellranger/7.0.0

#
## Gather args
args=("$@")

# Directories
id=${args[0]}
transcriptome=${args[1]}
fastqs=${args[2]}
localcores=${args[3]}
localmem=${args[4]}
sample=${args[5]}
OUT_DIR_SAMPLE=${args[6]}

echo "id=${id}"
echo "transcriptome=${transcriptome}"
echo "fastqs=${fastqs}"
echo "localcores=${localcores}"
echo "localmem=${localmem}"
echo "sample=${sample}"
echo "OUT_DIR_SAMPLE=${OUT_DIR_SAMPLE}"

time cellranger count --id="${sample}" \
                      --transcriptome=${transcriptome} \
                      --fastqs="fastqfiles/${FASTQS_DIR}" \
                      --localcores="${localcores}" \
                      --localmem="${localmem}" \
                      --sample=${sample}

mkdir -p $OUT_DIR_SAMPLE
mv $sample $OUT_DIR_SAMPLE
