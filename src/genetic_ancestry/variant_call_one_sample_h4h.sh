#!/bin/bash
#SBATCH --mem=10G
#SBATCH -p all
#SBATCH -c 32
#SBATCH -N 1
#SBATCH -t 1-00:00:00
#SBATCH --output=%x.out



module load singularity/3.5.2

echo "$1"



### Set up paths to files ###
BAM=/cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/${3}/${2}_possorted_genome_bam.bam
FASTA=/cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/hg38.p13.fa
OUT=/cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/${3}/${2}
SIF=/cluster/projects/macparland/RE/ancestry_calling/scadmix/scadmix.sif
BIND=/cluster/projects/macparland/RE


### Define the bed file for the specific chromosome being used ###
BED="/cluster/projects/macparland/RE/ancestry_calling/scadmix/chr/hgdp_wgs.20190516.full.subset.chr${1}.bed"


### Run scadmix
singularity exec --bind $BIND $SIF bash scadmix_variants.sh \
                                                -b $BAM \
                                                -c "chr${1}" \
                                                -e $BED \
                                                -f $FASTA \
                                                -o $OUT



