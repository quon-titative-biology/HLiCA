kinit -l 28d
aklog

module load sratoolkit/2.11.2

'''
Configure vdb-config
Enable local file-caching and set it to fastqfiles/sra_downloads as absolute path
Can increase RAM used to 100 MB
'''

vdb-config -i

"""

"""
cd fastqfiles/sra_downloads

for srr_id in $(cat SRR_Acc_List.txt)
do
  echo $srr_id
  if [ -f "${srr_id}_2.fastq.gz" ]
  then
    echo 'File exists'
  else
    echo 'File does not exist, start fastq-dump'
    fastq-dump --gzip --split-files ${srr_id}
  fi
done


for srr_id in $(cat SRR_Acc_List.txt | tac)
do
  echo $srr_id
  if [ -f "${srr_id}_2.fastq.gz" ]
  then
    echo 'File exists'
  else
    echo 'File does not exist, start fastq-dump'
    # fastq-dump --gzip --split-files ${srr_id}
    time fasterq-dump --split-files ${srr_id}
    time gzip "${srr_id}*.fastq"
  fi
done
