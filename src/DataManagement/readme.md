# Download the data from aws bucket and perform basic sequence alignment with cell ranger



# Pipeline
General pipeline including shell script for all components are included in the pipeline.sh

* Download reference genome
* Download fastqfiles
* Prepare metdata
* Align to reference genome

```
bash pipeline.sh
```

# Data structure
genome: directory where reference genome is stored
fastqfiles: directory where fastqfiles are stored
alignments: directory where cellranger counts outputs are stored

# Download reference genome
There are multiple options for reference genome, either GRCh38 or T2T. Additionally multiple versions of gencode can be used. Currently, using GENCODE v42.

```
bash download_reference.sh
```

# Download fastq files

## Download fastq files metadata and use aws
First download the metadata of fastq files from google drive in form of csv (fastqfiles/GoogleSheetMetadata.csv).

```
bash download_fastqfiles/download_with_lattice.sh
```

## Download fastqfiles with SRA from GEO

For datasets that are not yet available from lattice use SRA toolkit for download.

```
bash download_fastqfiles/download_with_sra.sh
```

# Prepare the metadata

Prepare the metadata.csv for lattice datasets

```
python prep_samples.py
```

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


# Align to ref with cellranger

```
bash run_cellranger.sh
```

# Reformat for sharing

Aggregate the results across samples per dataset and create a single h5ad file for sharing.

## 1. Create csv for cellranger aggr function

```
bash reformat_for_sharing/01_cellranger_aggr_create_csv.py
```

## 2. Run the cellranger aggr

```
02_cellranger_aggr.sh
```
