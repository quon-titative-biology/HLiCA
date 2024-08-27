


##################################
### GuilliamsScott
##################################
salloc -c 1 -t 1:0:0 --mem 1G
cd /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA
mkdir GuilliamsScott

scp /media/redgar/Seagate\ Portable\ Drive/HLiCA/GuilliamsScott/*_possorted_genome_bam.bam t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/GuilliamsScott


###################
#### Variant calling
###################


cd /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/GuilliamsScott
#mkdir MA_2022_5


# files from:
#http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p13/hg38.p13.fa.gz
#https://www.dropbox.com/s/gt9a8n6thzwaysk/GRCh38_chr_beds.tar.gz


###PERSAMPLE
sample=( 
    ABU8
    CISE06
    CISE07
    CISE08
    CISE09
    CS101
    CS108
    CS109
    CS110
    CS111
    CS112
    CS126
    CS127
    CS162
    CS164
    CS166
    CS167
    CS169
    CS170
    CS171
    CS31
    CS32
    CS33
    CS34
    CS37
    CS38
    CS41
    CS42
    CS43
    CS44
    CS46
    CS71
    CS73
    CS81
    CS83
    CS85
    CS87
)


####################
## Index BAM
####################
for i in "${sample[@]}"; do
    sbatch --job-name "index_$i" /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/index_bams.sh "$i" "GuilliamsScott"
done


####################
## Variant calling (all chromosomes one all samples)
####################
for i in "${sample[@]}"; do
    for chrom in $(seq 1 22); do
        job_name="${i}_${chrom}"
        sbatch --job-name "$job_name" /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/variant_call_one_sample_h4h.sh "$chrom" "$i" "GuilliamsScott"
    done
done



####################
#### Ancestral Admixture Inference
####################
## pruned for admixtures hgdp
#https://www.dropbox.com/scl/fi/wd6hwr3ygdd5mmdio7hcr/continental_w_middle_east-GRCh38.tar.gz?rlkey=kvr3ztknm2m7w5b45sd6wa4q3

for i in "${sample[@]}"; do
    sbatch --job-name "ansinf_$i" /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/ancestral_inference_one_sample_h4h.sh "$i" "GuilliamsScott"
done


## transfer output to make plot locally
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/GuilliamsScott/*/*_ancestry.qopt /home/redgar/Documents/ancestry_calling/scadmix/HLiCA











##################################
### Toronto
##################################
salloc -c 1 -t 1:0:0 --mem 1G
cd /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA
mkdir toronto

scp /media/redgar/Seagate\ Portable\ Drive/transplant_bams/SRR*_possorted_genome_bam.bam t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/toronto
scp /media/redgar/Seagate\ Portable\ Drive/transplant_bams/Human_Caudate_C72_3pr_v2_possorted_genome_bam.bam t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/toronto
scp /media/redgar/Seagate\ Portable\ Drive/transplant_bams/TLH_18_09_20_CITE_3_possorted_genome_bam.bam t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/toronto
scp /media/redgar/Seagate\ Portable\ Drive/transplant_bams/C70_Caudate_3pr_v2_possorted_genome_bam.bam t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/toronto




###################
#### Variant calling
###################


cd /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/toronto


# files from:
#http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p13/hg38.p13.fa.gz
#https://www.dropbox.com/s/gt9a8n6thzwaysk/GRCh38_chr_beds.tar.gz


###PERSAMPLE
sample=( 
    C70_Caudate_3pr_v2
    Human_Caudate_C72_3pr_v2
    TLH_18_09_20_CITE_3
    SRR16227558
    SRR16227559
    SRR16227560
    SRR16227570
    SRR16227577
    SRR16227584
    SRR7276476
    )



####################
## Index BAM
####################
for i in "${sample[@]}"; do
    sbatch --job-name "index_$i" /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/index_bams.sh "$i" "toronto"
done


####################
## Variant calling (all chromosomes one all samples)
####################
for i in "${sample[@]}"; do
    for chrom in $(seq 1 22); do
        job_name="${i}_${chrom}"
        sbatch --job-name "$job_name" /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/variant_call_one_sample_h4h.sh "$chrom" "$i" "toronto"
    done
done



####################
#### Ancestral Admixture Inference
####################
## pruned for admixtures hgdp
#https://www.dropbox.com/scl/fi/wd6hwr3ygdd5mmdio7hcr/continental_w_middle_east-GRCh38.tar.gz?rlkey=kvr3ztknm2m7w5b45sd6wa4q3

for i in "${sample[@]}"; do
    sbatch --job-name "ansinf_$i" /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/ancestral_inference_one_sample_h4h.sh "$i" "toronto"
done



## transfer output to make plot locally
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/toronto/*/*_ancestry.qopt /home/redgar/Documents/ancestry_calling/scadmix/HLiCA











##################################
### Henderson
##################################
salloc -c 1 -t 1:0:0 --mem 1G
cd /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA
mkdir Henderson

scp /media/redgar/Seagate\ Portable\ Drive/HLiCA/Henderson/*_possorted_genome_bam.bam t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/Henderson

###################
#### Variant calling
###################

cd /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/Henderson


# files from:
#http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p13/hg38.p13.fa.gz
#https://www.dropbox.com/s/gt9a8n6thzwaysk/GRCh38_chr_beds.tar.gz


###PERSAMPLE
sample=( 
  EDI001_End
  EDI001_Mes
  EDI002_End
  EDI002_Mes
  EDI003_End
  EDI003_Mes
  EDI004_End
  EDI004_Mes
  EDI005_End
  EDI005_Mes
  EDI006_End
  EDI006_Mes
  EDI007_End
  EDI007_Mes
  EDI008_End
  EDI008_Mes
  EDI009_End
  EDI009_Mes
  EDI010_End
  EDI010_Mes
  EDI011_End
  EDI011_Mes
  EDI012_End
  EDI012_Mes
  EDI013_End
  EDI013_Mes
  EDI014_EndMes
  EDI015_End
  EDI015_Mes
  EDI016_EndMes
  EDI017_End
)

####################
## Index BAM
####################
for i in "${sample[@]}"; do
    sbatch --job-name "index_$i" /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/index_bams.sh "$i" "Henderson"
done


####################
## Variant calling (all chromosomes one all samples)
####################
for i in "${sample[@]}"; do
    for chrom in $(seq 1 22); do
        job_name="${i}_${chrom}"
        sbatch --job-name "$job_name" /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/variant_call_one_sample_h4h.sh "$chrom" "$i" "Henderson"
    done
done



####################
#### Ancestral Admixture Inference
####################
## pruned for admixtures hgdp
#https://www.dropbox.com/scl/fi/wd6hwr3ygdd5mmdio7hcr/continental_w_middle_east-GRCh38.tar.gz?rlkey=kvr3ztknm2m7w5b45sd6wa4q3

for i in "${sample[@]}"; do
    sbatch --job-name "ansinf_$i" /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/ancestral_inference_one_sample_h4h.sh "$i" "Henderson"
done



## transfer output to make plot locally
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/Henderson/*/*_ancestry.qopt /home/redgar/Documents/ancestry_calling/scadmix/HLiCA









##################################
### Andrews
##################################
salloc -c 1 -t 1:0:0 --mem 1G
cd /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA
mkdir Andrews_2022

scp /media/redgar/Seagate\ Portable\ Drive/HLiCA/Andrews_2022/C*_possorted_genome_bam.bam t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/Andrews_2022


###################
#### Variant calling
###################

cd /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/Andrews_2022


# files from:
#http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p13/hg38.p13.fa.gz
#https://www.dropbox.com/s/gt9a8n6thzwaysk/GRCh38_chr_beds.tar.gz


###PERSAMPLE
sample=( 
  C41_NPC_SC_3pr_CBC03ANXX_bamtofastq
  C42_NPC_SC_3pr_CBC03ANXX_bamtofastq
  C43_NPC_SC_3pr_bamtofastq
  C46_SC_3pr_CC00WANXX_bamtofastq
  C48_SC_3pr_CBVE3ANXX_bamtofastq
  C49_SC_3pr_CC9M4ANXX_bamtofastq
  C50_SC_3pr_CC00HANXX_bamtofastq
  C51_flush_SC_3pr_CCFGAANXX_bamtofastq
  C51_SC_3pr_CCFP0ANXX_bamtofastq
  C52_SC_3pr_CC00HANXX_bamtofastq
  C53_SC_3pr_CCH7HANXX_bamtofastq
  C54_SC_3pr_CCJ7PANXX_bamtofastq
  C56_SC_3pr_CCKKGANXX_bamtofastq
  C58_SC_5pr_CD1TGANXX_bamtofastq
  C59_SC_3pr_CD10YANXX_bamtofastq
  C59_SC_5pr_CD1TGANXX_bamtofastq
  C61_SC_3pr_CD5ELANXX_bamtofastq
  C61_SC_5pr_CD3P5ANXX_bamtofastq
  C63_SC_3pr_CDP21ANXXbamtofastq
  C63_SC_5pr_CDP21ANXX_bamtofastq
  C64_SC_3pr_CDPFDANXX_bamtofastq
  C64_SC_5pr_bamtofastq
  C66_SC_3pr_CDPYWANXX_bamtofastq
  C68_SC_3pr_HTGHKDMXX_bamtofastq
  C69_SC_3pr_bamtofastq
  C70_SC_5pr_bamtofastq
)




####################
## Index BAM
####################
for i in "${sample[@]}"; do
    sbatch --job-name "index_$i" /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/index_bams.sh "$i" "Andrews_2022"
done


####################
## Variant calling (all chromosomes one all samples)
####################
for i in "${sample[@]}"; do
    for chrom in $(seq 1 22); do
        job_name="${i}_${chrom}"
        sbatch --job-name "$job_name" /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/variant_call_one_sample_h4h.sh "$chrom" "$i" "Andrews_2022"
    done
done



####################
#### Ancestral Admixture Inference
####################
## pruned for admixtures hgdp
#https://www.dropbox.com/scl/fi/wd6hwr3ygdd5mmdio7hcr/continental_w_middle_east-GRCh38.tar.gz?rlkey=kvr3ztknm2m7w5b45sd6wa4q3

for i in "${sample[@]}"; do
    sbatch --job-name "ansinf_$i" /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/ancestral_inference_one_sample_h4h.sh "$i" "Andrews_2022"
done



## transfer output to make plot locally
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/Andrews_2022/*/*_ancestry.qopt /home/redgar/Documents/ancestry_calling/scadmix/HLiCA



##################################
### MacParland
##################################
salloc -c 1 -t 1:0:0 --mem 1G
cd /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA
mkdir MacParland

scp /media/redgar/Seagate\ Portable\ Drive/HLiCA/MacParland/P*_bamtofastq_possorted_genome_bam.bam t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/MacParland


###################
#### Variant calling
###################


cd /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/MacParland
#mkdir MA_2022_5


# files from:
#http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p13/hg38.p13.fa.gz
#https://www.dropbox.com/s/gt9a8n6thzwaysk/GRCh38_chr_beds.tar.gz


###PERSAMPLE
sample=( 
    P1TLH_bamtofastq
    P2TLH_bamtofastq
    P4TLH_bamtofastq
    P5TLH_bamtofastq
)


####################
## Index BAM
####################
for i in "${sample[@]}"; do
    sbatch --job-name "index_$i" /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/index_bams.sh "$i" "MacParland"
done


####################
## Variant calling (all chromosomes one all samples)
####################
for i in "${sample[@]}"; do
    for chrom in $(seq 1 22); do
        job_name="${i}_${chrom}"
        sbatch --job-name "$job_name" /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/variant_call_one_sample_h4h.sh "$chrom" "$i" "MacParland"
    done
done



####################
#### Ancestral Admixture Inference
####################
## pruned for admixtures hgdp
#https://www.dropbox.com/scl/fi/wd6hwr3ygdd5mmdio7hcr/continental_w_middle_east-GRCh38.tar.gz?rlkey=kvr3ztknm2m7w5b45sd6wa4q3

for i in "${sample[@]}"; do
    sbatch --job-name "ansinf_$i" /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/ancestral_inference_one_sample_h4h.sh "$i" "MacParland"
done



## transfer output to make plot locally
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/MacParland/*/*_ancestry.qopt /home/redgar/Documents/ancestry_calling/scadmix/HLiCA




















#####################################
### DasGupta
#####################################
salloc -c 1 -t 1:0:0 --mem 1G
cd /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA
mkdir DasGupta

scp /media/redgar/Seagate\ Portable\ Drive/HLiCA/DasGupta/XH*possorted_genome_bam.bam t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/DasGupta



###################
#### Variant calling
###################


cd /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/DasGupta
#mkdir MA_2022_5


# files from:
#http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p13/hg38.p13.fa.gz
#https://www.dropbox.com/s/gt9a8n6thzwaysk/GRCh38_chr_beds.tar.gz


###PERSAMPLE
sample=( 
    XHH012
    XHL101
    XHL103
    XHL107
    XHL109
    XHL118
    XHL121
    XHL123
    XHL125
    XHL127
    XHL129
    XHL133
    XHL318
    XHL319
    XHL322
    XHL323
    XHL326
    XHL327
    XHL330
    XHL338
    XHL339
    XHL341
    XHL343
)


####################
## Index BAM
####################
for i in "${sample[@]}"; do
    sbatch --job-name "index_$i" /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/index_bams.sh "$i" "DasGupta"
done


####################
## Variant calling (all chromosomes one all samples)
####################
for i in "${sample[@]}"; do
    for chrom in $(seq 1 22); do
        job_name="${i}_${chrom}"
        sbatch --job-name "$job_name" /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/variant_call_one_sample_h4h.sh "$chrom" "$i" "DasGupta"
    done
done


####################
#### Ancestral Admixture Inference
####################
## pruned for admixtures hgdp
#https://www.dropbox.com/scl/fi/wd6hwr3ygdd5mmdio7hcr/continental_w_middle_east-GRCh38.tar.gz?rlkey=kvr3ztknm2m7w5b45sd6wa4q3

for i in "${sample[@]}"; do
    sbatch --job-name "ansinf_$i" /cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/ancestral_inference_one_sample_h4h.sh "$i" "DasGupta"
done


## transfer output to make plot locally
scp t117652uhn@h4huhndata1.uhnresearch.ca:/cluster/projects/macparland/RE/ancestry_calling/scadmix/HLiCA/DasGupta/*/*_ancestry.qopt /home/redgar/Documents/ancestry_calling/scadmix/HLiCA





