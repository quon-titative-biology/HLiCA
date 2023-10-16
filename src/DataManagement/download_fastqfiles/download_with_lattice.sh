"""
1. Download data

  a. Fastq files from aws cli
  (https://docs.google.com/spreadsheets/d/1l7XM7wxihdm5JqSRj0Il5piNMGoiINk0KpqyLf2U8Wo/edit#gid=0)
  Henderson
  (https://docs.google.com/spreadsheets/d/1Qh3HZIxBXnZg9B9vqPmHUOn1_-DJgGs5RZTB7JNVPTw/export?gid=2032049985&format=csv)
  Gruen
  (https://docs.google.com/spreadsheets/d/1Qh3HZIxBXnZg9B9vqPmHUOn1_-DJgGs5RZTB7JNVPTw/export?gid=0&format=csv)
"""

# ==============================================================================
# Download metadata from google sheet
# ==============================================================================

# Set
# fastqfiles: directory to store fastqfiles
fastqfiles="fastqfiles" # relative
# fastqfiles="/share/quonlab/workspaces/czi_liver_atlas/data/LiverNetworkData/fastqfiles" #absolute

# Set
# Groupname: output path
# googlesheet_link: url
# GoogleSheetFile: output directory/name
GroupName="${fastqfiles}/Henderson"
googlesheet_link="https://docs.google.com/spreadsheets/d/1l7XM7wxihdm5JqSRj0Il5piNMGoiINk0KpqyLf2U8Wo/export?gid=1189575288&format=csv"
GoogleSheetFile="$GroupName/GoogleSheetMetadata.csv"

GroupName="${fastqfiles}/Gruen"
googlesheet_link="https://docs.google.com/spreadsheets/d/1l7XM7wxihdm5JqSRj0Il5piNMGoiINk0KpqyLf2U8Wo/export?gid=654291382&format=csv"
GoogleSheetFile="$GroupName/GoogleSheetMetadata.csv"

GroupName="${fastqfiles}/Guiliams"
googlesheet_link="https://docs.google.com/spreadsheets/d/1l7XM7wxihdm5JqSRj0Il5piNMGoiINk0KpqyLf2U8Wo/export?gid=1028823323&format=csv"
GoogleSheetFile="$GroupName/GoogleSheetMetadata.csv"

GroupName="${fastqfiles}/Toronto"
googlesheet_link="https://docs.google.com/spreadsheets/d/1l7XM7wxihdm5JqSRj0Il5piNMGoiINk0KpqyLf2U8Wo/export?gid=0&format=csv"
GoogleSheetFile="$GroupName/GoogleSheetMetadata.csv"

GroupName="${fastqfiles}/DasGupta"
googlesheet_link="https://docs.google.com/spreadsheets/d/1l7XM7wxihdm5JqSRj0Il5piNMGoiINk0KpqyLf2U8Wo/export?gid=872171588&format=csv"
GoogleSheetFile="$GroupName/GoogleSheetMetadata.csv"
# For DasGupta run modify_filename.py before downloading

# Create directories
mkdir -p $GroupName

# metadata file
wget --output-file="logs.csv" $googlesheet_link -O ${GoogleSheetFile}
head -n 5 ${GoogleSheetFile}


# ==============================================================================
# Download fastqfiles
# ==============================================================================

kinit -l 28d
aklog

#### Verifyy access to aws
# https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-envvars.html#envvars-set
module load awscli/2.7.12

# Load ACCESS_KEY_ID
# (recommended alternative: set up login with aws config files)
# export AWS_ACCESS_KEY_ID=
# export AWS_SECRET_ACCESS_KEY=
# export AWS_DEFAULT_REGION=

# Verify accesss
aws configure list
aws sts get-caller-identity

# Show the top directory of lattice bucket
aws s3 ls s3://submissions-czi004liv

# Download all the files from GoogleSheetMetadata.csv
# This assumes that the s3_uri is located in the first column and csv has a heade
dos2unix $GoogleSheetFile

# Skip first line (header) then download fastqfiles
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
