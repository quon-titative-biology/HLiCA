"""
3. Perform sequence alignment with cellranger
https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count#cr-mkfastq
"""
module load cellranger/7.0.0

fastqs="fastqfiles/submissions-czi004liv/andrews_2021/PRJNA769141/SRR16227562"
sample="TLH_18_09_20_CITE_3"

fastqs="fastqfiles/submissions-czi004liv/andrews_2021/SRR7276476/McGilvery_Sonya__TLH_June_29_MissingLibrary_1_CBC03ANXX"
sample="SRR7276476"

fastqs="fastqfiles/submissions-czi004liv/andrews_2021/SRR16227585/MacParland_Catia__C73_D1_VIS_0_1_HTKCMDRXX"
sample="SRR16227585"

fastqs="fastqfiles/submissions-czi004liv/andrews_2021/SRR16227570/Count_Results58_0_1_H5LJKCCX2"
sample="SRR16227570"

fastqs="fastqfiles/submissions-czi004liv/andrews_2021/SRR16227585/MacParland_Catia__C73_D1_VIS_0_1_HTKCMDRXX"
sample="SRR16227585"


ref="T2T_chm_github"
transcriptome="genome/T2T/T2T_chm_github/T2T_chm_github"

cellranger count --id="${sample}_${ref}" \
                 --transcriptome=${transcriptome} \
                 --fastqs=${fastqs} \
                 --localcores=64 \
                 --localmem=128 \
                 --sample=${sample}


ref="GRCh38p13_gencode_v35"
transcriptome="genome/GRCh38/GRCh38_gencode/GRCh38p13_gencode_v35"

cellranger count --id="${sample}_${ref}" \
                --transcriptome=${transcriptome} \
                --fastqs=${fastqs} \
                --localcores=64 \
                --localmem=128 \
                --sample=${sample}


ref="T2T_refseq"
transcriptome="genome/T2T/GCF_009914755.1_T2T-CHM13v2.0_ncbi/T2T"

cellranger count --id="${sample}_${ref}" \
                 --transcriptome=${transcriptome} \
                 --fastqs=${fastqs} \
                 --localcores=64 \
                 --localmem=128 \
                 --sample=${sample}


ref="GRCh38p13_hg38_10x"
transcriptome="genome/GRCh38/refdata-gex-GRCh38-2020-A_10x"

cellranger count --id="${sample}_${ref}" \
                 --transcriptome=${transcriptome} \
                 --fastqs=${fastqs} \
                 --localcores=64 \
                 --localmem=128 \
                 --sample=${sample}
