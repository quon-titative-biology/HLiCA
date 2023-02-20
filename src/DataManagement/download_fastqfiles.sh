"""
1. Download data

  a. Fastq files from aws cli
  (https://docs.google.com/spreadsheets/d/1l7XM7wxihdm5JqSRj0Il5piNMGoiINk0KpqyLf2U8Wo/edit#gid=0)

"""

# a. download fastq files

# metadata file
wget --output-file="logs.csv" "https://docs.google.com/spreadsheets/d/1l7XM7wxihdm5JqSRj0Il5piNMGoiINk0KpqyLf2U8Wo/export?format=csv&gid=0" -O "GoogleSheetMetadata.csv"
head -n 1 "GoogleSheetMetadata.csv"

# https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-envvars.html#envvars-set
module load awscli/2.7.12

export AWS_ACCESS_KEY_ID=AKIAYQQPUH52SXUYGWXA
export AWS_SECRET_ACCESS_KEY=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
export AWS_DEFAULT_REGION=us-west-2

# Verify accesss
aws configure list
aws sts get-caller-identity

aws s3 ls s3://submissions-czi004liv
aws s3 ls s3://submissions-czi004liv/andrews_2021 --recursive

# Download all the files from GoogleSheetMetadata.csv
# This assumes that the s3_uri is located in the first column and csv has a header

sed 1d "GoogleSheetMetadata.csv" |
while IFS=, read -r s3_uri rest; do

  # Set download variables
  echo "$s3_uri"

  # Set
  filename=${s3_uri:5} # strip s3://
  download_dir="$(dirname ${filename})" # get dirname
  echo $download_dir

  mkdir -p $download_dir

  # Copy from aws
  aws s3 cp $s3_uri $download_dir

done
