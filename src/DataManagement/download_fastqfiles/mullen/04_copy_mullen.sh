#### Download fastq data from mullen lab

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
# This assumes that the s3_uri is located in the first column and csv has a header

aws s3 ls s3://submissions-czi004liv/mullen_2023/6854_1/

src="s3://submissions-czi004liv/mullen_2023"
dest="alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13/share/Mullen/"
mkdir -p $dest

aws s3 cp $src $dest --recursive
