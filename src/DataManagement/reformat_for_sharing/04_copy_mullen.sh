#### Download fastq data from mullen lab

src="s3://submissions-czi004liv/mullen_2023"
dest="alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13/share/"
aws s3 cp $src $dest --recursive
