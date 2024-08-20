#!/bin/bash
#SBATCH --mem=50G
#SBATCH -p himem
#SBATCH -c 32
#SBATCH -N 1
#SBATCH -t 7-00:00:00
#SBATCH --output=%x.out

module load samtools
BAM=/cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/${2}/${1}_possorted_genome_bam.bam
samtools index $BAM

