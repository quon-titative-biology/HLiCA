"""
1. Download data

  a. Fastq files from aws cli
  (https://docs.google.com/spreadsheets/d/1l7XM7wxihdm5JqSRj0Il5piNMGoiINk0KpqyLf2U8Wo/edit#gid=0)
  Henderson
  (https://docs.google.com/spreadsheets/d/1Qh3HZIxBXnZg9B9vqPmHUOn1_-DJgGs5RZTB7JNVPTw/export?gid=2032049985&format=csv)
  Gruen
  (https://docs.google.com/spreadsheets/d/1Qh3HZIxBXnZg9B9vqPmHUOn1_-DJgGs5RZTB7JNVPTw/export?gid=0&format=csv)
"""

# Specify information
fastqfiles="fastqfiles" # relative
fastqfiles="/share/quonlab/workspaces/czi_liver_atlas/data/LiverNetworkData/fastqfiles" #absolute

GroupName="${fastqfiles}/Henderson"
googlesheet_link="https://docs.google.com/spreadsheets/d/1Qh3HZIxBXnZg9B9vqPmHUOn1_-DJgGs5RZTB7JNVPTw/export?gid=2032049985&format=csv"
GoogleSheetFile="$GroupName/GoogleSheetMetadata.csv"

GroupName="${fastqfiles}/Gruen"
googlesheet_link="https://docs.google.com/spreadsheets/d/1Qh3HZIxBXnZg9B9vqPmHUOn1_-DJgGs5RZTB7JNVPTw/export?gid=0&format=csv"
GoogleSheetFile="$GroupName/GoogleSheetMetadata.csv"

GroupName="${fastqfiles}/Guiliams"
googlesheet_link="https://docs.google.com/spreadsheets/d/1l7XM7wxihdm5JqSRj0Il5piNMGoiINk0KpqyLf2U8Wo/export?gid=1028823323&format=csv"
GoogleSheetFile="$GroupName/GoogleSheetMetadata.csv"

# Set directories
# cd fastqfiles
mkdir -p $GroupName
# cd $GroupName

# a. download fastq files

# metadata file
wget --output-file="logs.csv" $googlesheet_link -O ${GoogleSheetFile}
head -n 5 ${GoogleSheetFile}

# sed  -i '1i s3_uri' ${GoogleSheetFile} # insert s3_uri as first line
# head -n 5 ${GoogleSheetFile}

#### Verifyy access to aws
# https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-envvars.html#envvars-set
module load awscli/2.7.12

# Load ACCESS_KEY_ID
export AWS_ACCESS_KEY_ID=
export AWS_SECRET_ACCESS_KEY=
export AWS_DEFAULT_REGION=

# Verify accesss
aws configure list
aws sts get-caller-identity

# Show the top directory of lattice bucket
aws s3 ls s3://submissions-czi004liv

## Sync the entire aws directory
# aws s3 sync "s3://submissions-czi004liv/gruen_2023/" "fastqfiles/submissions-czi004liv/gruen_2023"
# aws s3 sync "s3://submissions-czi004liv/henderson_2023/" "fastqfiles/submissions-czi004liv/henderson_2023"
#
# aws s3 ls "s3://submissions-czi004liv/henderson_2023/" --recursive --human-readable --summarize | grep "fastq.gz"

# Download all the files from GoogleSheetMetadata.csv
# This assumes that the s3_uri is located in the first column and csv has a header


dos2unix $GoogleSheetFile

sed 1d ${GoogleSheetFile} |
while IFS=, read -r s3_uri rest; do

  # Set download variables
  echo "s3_uri: ${s3_uri}"

  # Set download file names
  filename=${s3_uri:5} # strip s3://
  download_dir_file="${fastqfiles}/${filename}"
  download_dir="${fastqfiles}/$(dirname ${filename})" # get dirname
  echo $download_dir_file
  # touch $download_dir_file
  if [ -f $download_dir_file ]
  then
      echo "File exist, skip download"
  else
      echo "File does not exist, download"

      mkdir -p $download_dir

      # Copy from aws
      echo "$s3_uri"
      aws s3 cp "$s3_uri" "$download_dir"
  fi


  echo "Finished"
  echo ""

done
