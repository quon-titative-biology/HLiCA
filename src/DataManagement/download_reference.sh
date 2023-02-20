"""

Download reference genome and create ref file for 10x
Reference genome (T2T) from ncbi (refseq)

(https://www.ncbi.nlm.nih.gov/assembly/GCA_009914755.4)
(https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/)

https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count#cr-mkfastq

"""

module load cellranger/7.0.0

"""
T2T from ncbi (refseq)
(https://www.ncbi.nlm.nih.gov/assembly/GCA_009914755.4)
(https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/)
"""

rsync_path="rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0"
rsync --copy-links --recursive --times --verbose "${rsync_path}/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz" ./
rsync --copy-links --recursive --times --verbose "${rsync_path}/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz" ./

file_list=("GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz" \
           "GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz" \
           "GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz")

for file in "${file_list[@]}"
  do

  echo $file
  rsync --copy-links --recursive --times --verbose "${rsync_path}/${file}" ./
  gzip -d ${file}
done

# convert to cell ranger format
cellranger mkref --genome=T2T \
                 --fasta=GCF_009914755.1_T2T-CHM13v2.0_genomic.fna \
                 --genes=GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf \
                 --nthreads=32

"""
T2T-CHM13v2.0 from github
https://github.com/marbl/CHM13
"""

cd T2T
mkdir T2T_chm_github

curl -O https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
curl -O https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_GENCODEv35_CAT_Liftoff.vep.gff3.gz
curl -O https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_GENCODEv35_CAT_Liftoff.vep.gff3.gz.tbi

# run gff2gtf.py

gzip -d chm13v2.0_GENCODEv35_CAT_Liftoff.vep.gtf.gz
gzip -d chm13v2.0.fa.gz

cellranger mkref --genome=T2T_chm_github \
                 --fasta=chm13v2.0.fa \
                 --genes=chm13v2.0_GENCODEv35_CAT_Liftoff.vep.gtf \
                 --nthreads=32


"""
Download reference genome, grch38 from 10x
"""

curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -xzvf refdata-gex-GRCh38-2020-A.tar.gz

"""
Download from gencode
"""
curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/GRCh38.p13.genome.fa.gz
curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gtf.gz

gzip -d GRCh38.p13.genome.fa.gz
gzip -d gencode.v35.annotation.gtf.gz

cellranger mkref --genome=GRCh38p13_gencode_v35 \
                 --fasta=GRCh38.p13.genome.fa \
                 --genes=gencode.v35.annotation.gtf \
                 --nthreads=32

"""
Download from gencode
genome/GRCh38/GRCh38_gencode/v42
"""
ver=42

curl -O "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${ver}/GRCh38.p13.genome.fa.gz"
curl -O "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${ver}/gencode.v42.annotation.gtf.gz"

gzip -d GRCh38.p13.genome.fa.gz
gzip -d gencode.v42.annotation.gtf.gz

cellranger mkref --genome=GRCh38p13_gencode_v42 \
                --fasta=GRCh38.p13.genome.fa \
                --genes="gencode.v${ver}.annotation.gtf" \
                --nthreads=32
