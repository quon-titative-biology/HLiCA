#!/bin/bash
#SBATCH --mem=10G
#SBATCH -p all
#SBATCH -c 32
#SBATCH -N 1
#SBATCH -t 1-00:00:00
#SBATCH --output=%x.out



module load singularity/3.5.2

echo "$1"


#### Ancestral Admixture Inference
## pruned for admixtures hgdp
#https://www.dropbox.com/scl/fi/wd6hwr3ygdd5mmdio7hcr/continental_w_middle_east-GRCh38.tar.gz?rlkey=kvr3ztknm2m7w5b45sd6wa4q3

SIF=/cluster/projects/macparland/RE/ancestry_calling/scadmix/scadmix.sif
N_IND_FILE=/cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/continental_w_middle_east-GRCh38/nInd_hgdp_wgs.20190516.full.subset.txt
FREQ_FILE=/cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/continental_w_middle_east-GRCh38/refPanel_hgdp_wgs.20190516.full.subset.txt
FILE_NAME=${1}_ancestry
OUT=/cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/${2}/${1}
BIND=/cluster/projects/macparland/RE


singularity exec --bind $BIND $SIF bash scadmix_ancestry.sh \
                                    -r "True" \
                                    -o $OUT \
                                    -n $FILE_NAME \
                                    -f $FREQ_FILE \
                                    -i $N_IND_FILE