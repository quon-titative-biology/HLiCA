"""
Cell ranger version 7.0.0
"""

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

"""
Run cellranger count
"""

transcriptome="GRCh38p13_gencode_v42"
sample=
fastqs=

cellranger count --id="${sample}" \
                 --transcriptome=${transcriptome} \
                 --fastqs=${fastqs} \
                 --localcores=32 \
                 --localmem=256 \
                 --sample=${sample} \
                 --no-bam
