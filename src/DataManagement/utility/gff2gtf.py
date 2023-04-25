import os
import gzip

gff_file_dir = "chm13v2.0_GENCODEv35_CAT_Liftoff.vep.gff3.gz"

gtf_file_dir = "chm13v2.0_GENCODEv35_CAT_Liftoff.vep.gtf.gz"

if gtf_file_dir is None:
    gtf_file_dir = gff_file_dir.replace('gff3','gtf')

def gff2gtf_perline(line):
    # line = line.decode("utf-8")
    # line = line.strip() # remove \n
    line = line.replace(b";",b"; ").replace(b"=",b" ")
    return line
    # seqid,source,feature,start,end,score,strand,phase,attributes = line.decode("utf-8") .split('\t')
    #
    # attr_dict = {attr.split('=')[0]:attr.split('=')[1] for attr in attributes.split(';')}


with gzip.open(gtf_file_dir,'wb') as gtf_file:
    gtf_file.write(b"##description: converted from chm13v2.0_GENCODEv35_CAT_Liftoff.vep.gff3\n")
    gtf_file.write(b"##provider: UCSC GENCODEv35 CAT/Liftoff v2\n")
    gtf_file.write(b"##contact: NA\n")
    gtf_file.write(b"##format: gtf\n")
    gtf_file.write(b"##date: 2023-02-08\n")
    with gzip.open(gff_file_dir,'r') as gff_file:
        for i,line in enumerate(gff_file):
            line_new = gff2gtf_perline(line)
            gtf_file.write(line_new)

"""
Example gtf

##description: evidence-based annotation of the human genome (GRCh38), version 35 (Ensembl 101)
##provider: GENCODE
##contact: gencode-help@ebi.ac.uk
##format: gtf
##date: 2020-06-03
chr1    HAVANA  gene    11869   14409   .       +       .       gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; hgnc_id "HGNC:37102"; havana_gene "OTTHUMG00000000961.2";
chr1    HAVANA  transcript      11869   14409   .       +       .       gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1    HAVANA  exon    11869   12227   .       +       .       gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; exon_number 1; exon_id "ENSE00002234944.1"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1    HAVANA  exon    12613   12721   .       +       .       gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; exon_number 2; exon_id "ENSE00003582793.1"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
chr1    HAVANA  exon    13221   14409   .       +       .       gene_id "ENSG00000223972.5"; transcript_id "ENST00000456328.2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; transcript_type "processed_transcript"; transcript_name "DDX11L1-202"; exon_number 3; exon_id "ENSE00002312635.1"; level 2; transcript_support_level "1"; hgnc_id "HGNC:37102"; tag "basic"; havana_gene "OTTHUMG00000000961.2"; havana_transcript "OTTHUMT00000362751.1";
"""
