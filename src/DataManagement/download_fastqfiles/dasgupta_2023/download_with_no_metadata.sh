"""
Copy the entire sumbissions-czi004liv/dasgupta_2023
"""

s3_uri="s3://submissions-czi004liv/dasgupta_2023"
aws s3 ls "$s3_uri"


download_dir="fastqfiles/submissions-czi004liv/"
mkdir -p $download_dir
aws s3 cp --recursive "$s3_uri" "$download_dir"
