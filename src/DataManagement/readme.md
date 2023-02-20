# Download the data from aws bucket and perform basic sequence alignment
# with cell ranger

# Pipeline
bash pipeline.sh

# Data structure
genome
fastqfiles
alignments


## Download reference genome
There are multiple options, either GRCh38 or T2T. Additionally multiple versions of gencode can be used.

bash download_reference.sh


## Download fastq files from aws

### Download fastqfiles metadata and use aws
First download the metadata of fastq files from google drive in form of csv (fastqfiles/GoogleSheetMetadata.csv).

bash download_fastqfiles.sh

### Prepare the metadata.csv
python prep_samples.py

This creates GoogleSheetMetadata_file.csv and GoogleSheetMetadata_sample.csv. GoogleSheetMetadata_file.csv adds 3 columns, [SAMPLE, FASTQS_DIR, FILE_DIR]. SAMPLE and FASTQS_DIR are arguments to be used for cellranger count.

s3 uri
"s3://submissions-czi004liv/andrews_2021/SRR7276476/McGilvery_Sonya__TLH_June_29_MissingLibrary_1_CBC03ANXX/SRR7276476_S1_L007_R2_001.fastq.gz"
FILE_DIR
"submissions-czi004liv/andrews_2021/SRR7276476/McGilvery_Sonya__TLH_June_29_MissingLibrary_1_CBC03ANXX/SRR7276476_S1_L007_R2_001.fastq.gz"
FASTQS_DIR
"submissions-czi004liv/andrews_2021/SRR7276476/McGilvery_Sonya__TLH_June_29_MissingLibrary_1_CBC03ANXX"
SAMPLE
"SRR7276476"

_file.csv is metadata per file with updated columns. _sample.csv removes duplicate in columns FASTQS_DIR and SAMPLE so that each row is sample for cellranger.


## Align to ref with cellranger
bash run_cellranger.sh
