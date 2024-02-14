# ESB
.libPaths()
library(Matrix)
library(data.table)
library(Seurat)
library(SeuratDisk)

path <- '/share/quonlab/HLiCA/czi_liver_atlas/integration/data'
meta_path <- '/share/quonlab//HLiCA/czi_liver_atlas/integration/metadata'
setwd(path)

datasets <- c('Mullen', 'DasGupta', 'Gruen',
              'Henderson','GuilliamsScott',
              'MacParland', 'Toronto', 'Andrews')
covariates <- c('STUDY', 'suspension_type','assay', 'donor_uuid', 'library_alias')
meta_col <- c('STUDY', 'SAMPLE', 'tissue', 
               'suspension_type','assay', 'sequencing_run',
                 'library_alias','library_uuid', 
                   'donor_uuid','donor_sex', 'donor_age', 'donor_ethnicity')
sample_num <- list()
sequencing_run_num <- list()
donor_num <- list()
library_num <- list()
cell_num <- list()
for (study in datasets){
    print(study)
    study_path = file.path(path, study)
    raw_metadata = read.csv(file.path(meta_path,study,'GoogleSheetMetadata_sample_qced_data.csv'), header=TRUE, sep=",")
    raw_metadata$AlignmentCellRanger <- toupper(raw_metadata$AlignmentCellRanger)
    raw_metadata$QC_Pipeline <- toupper(raw_metadata$QC_Pipeline)
    metadata = raw_metadata[raw_metadata$AlignmentCellRanger=='TRUE',]
    metadata = metadata[metadata$QC_Pipeline=='TRUE',]
    metadata$STUDY = study
    metadata = metadata[meta_col]
    # Toronto to Toronto_2
    # Andrews to Toronto_3
    # MacParland to Toronto_1
    # Gruen to Grun
    # GuilliamsScott to Scott
    if (study == 'Toronto'){
        metadata$STUDY = 'Toronto_2'
    } else if (study == 'Andrews'){
        metadata$STUDY = 'Toronto_3'
    } else if (study == 'MacParland'){
        metadata$STUDY = 'Toronto_1'
    } else if (study == 'Gruen'){
        metadata$STUDY = 'Grun'
    } else if (study == 'GuilliamsScott'){
        metadata$STUDY = 'Scott'
    }
    print(head(metadata))
    print(table(metadata$SAMPLE))
    if (sum(table(metadata$SAMPLE)>1)){
        # break the loop
        print('Error: Duplicate sample names')
        break
    }
    # obtain samples to select
    SAMPLE_ID = metadata$SAMPLE
    rownames(metadata) = metadata$SAMPLE
    study_list = list()
    for (sample in SAMPLE_ID){
        print(sample)
        filename = Sys.glob(file.path(path,study,paste0(sample, '_RNA.rds')))
        processed_RNA = readRDS(filename)
        print(dim(processed_RNA))
        sample_metadata = processed_RNA@meta.data
        sample_metadata = sample_metadata[c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mt',
                                            'cell_barcode', 'seurat_clusters', 'S.Score', 'G2M.Score', 'Phase', 
                                            'scmap_cluster_anno', 'scmap_cell_anno', 'general_labs', 'consistent_labs', 'marker_labs', 'marker_general_labs','SAMPLE')]
        metadata_ = metadata[names(table(sample_metadata$SAMPLE)),]
        for (meta_col_ in meta_col){
            # print(meta_col_)
            sample_metadata[meta_col_] = metadata_[meta_col_]
        }
        processed_RNA@meta.data = sample_metadata
        study_list[[sample]] = processed_RNA
    }
    # length(study_list)
    # study_list
    merged_seurat <- merge(x = study_list[[1]],
            y = study_list[2:length(study_list)],
            merge.data = TRUE)
    # remove duplicated samples in Andrews dataset (Toronto_3)
    # C46_SC_3pr_HLV55DRXX_bamtofastq
    # C61_SC_3pr_HTGHKDMXX_bamtofastq
    # C63_SC_3pr_HLV55DRXX_bamtofastq
    # C63_SC_5pr_HLV55DRXX_bamtofastq
    # C68_SC_3pr_CE0GVANXX_bamtofastq
    if (study == 'Andrews'){
        rm_list <- c('C46_SC_3pr_HLV55DRXX_bamtofastq', 
        'C61_SC_3pr_HTGHKDMXX_bamtofastq', 
        'C63_SC_3pr_HLV55DRXX_bamtofastq', 
        'C63_SC_5pr_HLV55DRXX_bamtofastq', 
        'C68_SC_3pr_CE0GVANXX_bamtofastq')
        # remove the duplicated sample
        merged_seurat = merged_seurat[,-which(merged_seurat@meta.data$SAMPLE%in%rm_list)]
    }
    metadata = merged_seurat@meta.data
    sample_num[study] = dim(table(metadata$SAMPLE))
    sequencing_run_num[study] = dim(table(metadata$sequencing_run))
    donor_num[study] = dim(table(metadata$donor_uuid))
    library_num[study] = dim(table(metadata$library_alias))
    saveRDS(merged_seurat, file.path(path, paste0(study, '_RNA_merged.rds')))
    fwrite(as.data.frame(merged_seurat@meta.data), file.path(meta_path, paste0(study, '_RNA_merged_metadata.csv')), row.names=TRUE, col.names=TRUE, sep=',')
    cell_num[study] = dim(merged_seurat)[2]
}

# cell number of all duplicated samples in the Andrews dataset (Toronto_3)
# C46_SC_3pr_CC00WANXX_bamtofastq C46_SC_3pr_HLV55DRXX_bamtofastq
#                            1612.                           3760
# C61_SC_3pr_CD5ELANXX_bamtofastq C61_SC_3pr_HTGHKDMXX_bamtofastq
#                           10719.                           5036
#  C63_SC_3pr_CDP21ANXXbamtofastq C63_SC_3pr_HLV55DRXX_bamtofastq
#                           12760.                           3830
# C63_SC_5pr_CDP21ANXX_bamtofastq C63_SC_5pr_HLV55DRXX_bamtofastq
#                           15798.                          17318
# C68_SC_3pr_CE0GVANXX_bamtofastq C68_SC_3pr_HTGHKDMXX_bamtofastq
#                           13414                           13782.

# to be removed object in the Andrews dataset (Toronto_3) [confirmed by Rachel]
# table(merged_seurat[,which(merged_seurat@meta.data$SAMPLE%in%rm_list)]@meta.data$SAMPLE)
# C46_SC_3pr_HLV55DRXX_bamtofastq C61_SC_3pr_HTGHKDMXX_bamtofastq
#                            3760                            5036
# C63_SC_3pr_HLV55DRXX_bamtofastq C63_SC_5pr_HLV55DRXX_bamtofastq
#                            3830                           17318
# C68_SC_3pr_CE0GVANXX_bamtofastq
#                           13414


sample_num
sequencing_run_num
donor_num
library_num
cell_num
# make the four lists into a dataframe
sample_num = as.data.frame(sample_num)
rownames(sample_num) = 'sample'
sequencing_run_num = as.data.frame(sequencing_run_num)
rownames(sequencing_run_num) = 'sequencing_run'
donor_num = as.data.frame(donor_num)
rownames(donor_num) = 'donor'
library_num = as.data.frame(library_num)
rownames(library_num) = 'library'
cell_num = as.data.frame(cell_num)
rownames(cell_num) = 'cell'
# merge the four dataframes
stats = Reduce('rbind', list(sample_num, sequencing_run_num, donor_num, library_num, cell_num))
colnames(stats) = c('Mullen', 'DasGupta', 'Grün', 'Henderson','Scott', 'MacParland2018/Toronto1', 'Andrews2022/Toronto2', 'Andrews2024/Toronto3')
fwrite(stats, file.path(meta_path, 'sample_stats.csv'), row.names=TRUE, col.names=TRUE, sep=',')
t(stats)
#                         sample sequencing_run donor library   cell
# Mullen                      23             23    23      23 137360
# DasGupta                    23             23    23      23  49921
# Grün                        4              4     4       4  29242
# Henderson                   29             29    17      29 105525
# Scott                       37             37    19      37  59201
# MacParland2018/Toronto1      4              4     4       4   6603
# Andrews2022/Toronto2        10             10     4      10  80130
# Andrews2024/Toronto3        26             26    21      26 176698
rowSums(stats)
# sample sequencing_run          donor        library           cell
#    156            156            115            156         644680




######################################################################
######################################################################
######################################################################
######################################################################
study_list = list()
common_gene_list = list()
for (study in datasets){
    print(study)
    obj = readRDS(file.path(path, paste0(study, '_RNA_merged.rds')))
    counts <- GetAssayData(object=obj, slot = "counts")
    nonzero <- counts > 0
    keep_genes <- Matrix::rowSums(nonzero) >= min(100, 0.01*dim(obj)[2]) # keep genes with at least expressed over 1% or 100 cells, whichever is smaller
    filtered_counts <- counts[keep_genes, ]
    genes_kept <- rownames(filtered_counts)
    print("Number of genes kept")
    print(length(genes_kept))
    common_gene_list[[study]] = genes_kept
    study_list[[study]] = obj
}
# [1] "Mullen"
# [1] "Number of genes kept"
# [1] 30480
# [1] "DasGupta"
# [1] "Number of genes kept"
# [1] 24770
# [1] "Gruen"
# [1] "Number of genes kept"
# [1] 20461
# [1] "Henderson"
# [1] "Number of genes kept"
# [1] 28871
# [1] "GuilliamsScott"
# [1] "Number of genes kept"
# [1] 24706
# [1] "MacParland"
# [1] "Number of genes kept"
# [1] 14787
# [1] "Toronto"
# [1] "Number of genes kept"
# [1] 25518
# [1] "Andrews"
# [1] "Number of genes kept"
# [1] 24034
saveRDS(common_gene_list, file.path(path, 'common_gene_list.rds'))
common_genes = Reduce(union, common_gene_list)
length(common_genes)
# [1] 33670

merged_seurat <- merge(x = study_list[[1]],
            y = study_list[2:length(study_list)],
            merge.data = TRUE)
print(dim(merged_seurat))
# [1]  62696 644680

# change Unknown to unknown
merged_seurat@meta.data[which(merged_seurat@meta.data$donor_ethnicity == 'Unknown'), ]$donor_ethnicity = 'unknown'


# Toronto 1 dataset metadata modification
# SRR7276476: C41/P3TLH --> change the study to Toronto_1
merged_seurat@meta.data[which(merged_seurat@meta.data$SAMPLE == 'SRR7276476'), ]$SAMPLE = 'P3TLH'
merged_seurat@meta.data[which(merged_seurat@meta.data$SAMPLE == 'P3TLH'), ]$STUDY = 'Toronto_1'
merged_seurat@meta.data[which(merged_seurat@meta.data$SAMPLE == 'P1TLH_bamtofastq'), ]$SAMPLE = 'P1TLH'
merged_seurat@meta.data[which(merged_seurat@meta.data$SAMPLE == 'P2TLH_bamtofastq'), ]$SAMPLE = 'P2TLH'
merged_seurat@meta.data[which(merged_seurat@meta.data$SAMPLE == 'P4TLH_bamtofastq'), ]$SAMPLE = 'P4TLH'
merged_seurat@meta.data[which(merged_seurat@meta.data$SAMPLE == 'P5TLH_bamtofastq'), ]$SAMPLE = 'P5TLH'
# toronto_data_2018 <- merged_seurat[,which(merged_seurat@meta.data$STUDY=='Toronto_1')]
# table(toronto_data_2018@meta.data$SAMPLE)
# P1TLH P2TLH P3TLH P4TLH P5TLH
#  1206  1437  6875  1224  2736

# Toronto 2 dataset metadata modification
# change the names to:
# C70_Caudate_3pr_v2: C70_RESEQ
# Human_Caudate_C72_3pr_v2:  C72_RESEQ
# SRR16227558:  C41_CST
# SRR16227559:  C41_NST
# SRR16227560:  C41_TST
# SRR16227570:  C58_TST
# SRR16227577:  C70_TST
# SRR16227584:  C72_TST
# TLH_18_09_20_CITE_3: C58_RESEQ
merged_seurat@meta.data[which(merged_seurat@meta.data$SAMPLE == 'C70_Caudate_3pr_v2'), ]$SAMPLE = 'C70_RESEQ'
merged_seurat@meta.data[which(merged_seurat@meta.data$SAMPLE == 'Human_Caudate_C72_3pr_v2'), ]$SAMPLE = 'C72_RESEQ'
merged_seurat@meta.data[which(merged_seurat@meta.data$SAMPLE == 'SRR16227558'), ]$SAMPLE = 'C41_CST'
merged_seurat@meta.data[which(merged_seurat@meta.data$SAMPLE == 'SRR16227559'), ]$SAMPLE = 'C41_NST'
merged_seurat@meta.data[which(merged_seurat@meta.data$SAMPLE == 'SRR16227560'), ]$SAMPLE = 'C41_TST'
merged_seurat@meta.data[which(merged_seurat@meta.data$SAMPLE == 'SRR16227570'), ]$SAMPLE = 'C58_TST'
merged_seurat@meta.data[which(merged_seurat@meta.data$SAMPLE == 'SRR16227577'), ]$SAMPLE = 'C70_TST'
merged_seurat@meta.data[which(merged_seurat@meta.data$SAMPLE == 'SRR16227584'), ]$SAMPLE = 'C72_TST'
merged_seurat@meta.data[which(merged_seurat@meta.data$SAMPLE == 'TLH_18_09_20_CITE_3'), ]$SAMPLE = 'C58_RESEQ'
toronto_data_2022 <- merged_seurat[,which(merged_seurat@meta.data$STUDY=='Toronto_2')]
table(toronto_data_2022@meta.data$SAMPLE)
#   C41_CST   C41_NST   C41_TST C58_RESEQ   C58_TST C70_RESEQ   C70_TST C72_RESEQ  C72_TST
#      6715     13017      8214      4123      3326      6519     13357      8988     8996


# Toronto 3 dataset metadata modification
# change the names to:
# C41_NPC_SC_3pr_CBC03ANXX_bamtofastq: C41_NPC
# C42_NPC_SC_3pr_CBC03ANXX_bamtofastq: C42_NPC
# C43_NPC_SC_3pr_bamtofastq: C43_NPC
# C46_SC_3pr_CC00WANXX_bamtofastq: ?
# C46_SC_3pr_HLV55DRXX_bamtofastq: ? x
# C48_SC_3pr_CBVE3ANXX_bamtofastq: C48
# C49_SC_3pr_CC9M4ANXX_bamtofastq: C49
# C50_SC_3pr_CC00HANXX_bamtofastq: C50
# C51_flush_SC_3pr_CCFGAANXX_bamtofastq: C51_Flush
# C51_SC_3pr_CCFP0ANXX_bamtofastq: C51
# C52_SC_3pr_CC00HANXX_bamtofastq: C52
# C53_SC_3pr_CCH7HANXX_bamtofastq: C53_RESEQ
# C54_SC_3pr_CCJ7PANXX_bamtofastq: C54
# C56_SC_3pr_CCKKGANXX_bamtofastq: C56_RESEQ
# C58_SC_5pr_CD1TGANXX_bamtofastq: C58_5pr
# C59_SC_3pr_CD10YANXX_bamtofastq: C59
# C59_SC_5pr_CD1TGANXX_bamtofastq: C59_5pr
# C61_SC_3pr_CD5ELANXX_bamtofastq: ?
# C61_SC_3pr_HTGHKDMXX_bamtofastq: ? x
# C61_SC_5pr_CD3P5ANXX_bamtofastq: C61_5pr
# C63_SC_3pr_CDP21ANXXbamtofastq: ?
# C63_SC_3pr_HLV55DRXX_bamtofastq: ? x
# C63_SC_5pr_CDP21ANXX_bamtofastq: ?
# C63_SC_5pr_HLV55DRXX_bamtofastq: ? x
# C64_SC_3pr_CDPFDANXX_bamtofastq: C64_RESEQ
# C64_SC_5pr_bamtofastq: C64_5pr
# C66_SC_3pr_CDPYWANXX_bamtofastq: C66_RESEQ
# C68_SC_3pr_CE0GVANXX_bamtofastq: ? x
# C68_SC_3pr_HTGHKDMXX_bamtofastq: ?
# C69_SC_3pr_bamtofastq: C69
# C70_SC_5pr_bamtofastq: C70_5pr_RESEQ







head(merged_seurat@meta.data)
covariates <- c('STUDY', 'suspension_type','assay', 'donor_uuid', 'library_alias')

table(merged_seurat@meta.data$STUDY)
#  DasGupta      Grun Henderson    Mullen     Scott Toronto_1 Toronto_2 Toronto_3
#     49921     29242    105525    137360     59201     13478     73255    176698
table(merged_seurat@meta.data$suspension_type)
#    cell nucleus
#  450377  194303
table(merged_seurat@meta.data$assay)
# 10x 3' v2 10x 3' v3 10x 5' v1
#    202221    384986     57473
table(merged_seurat@meta.data$donor_sex)
# female   male
# 332803 311877
table(merged_seurat@meta.data$donor_ethnicity)
            #        African                    Chinese
            #           6740                      12569
            #       European Hispanic or Latin American
            #         262033                       1612
            #      Malaysian        Singaporean Chinese
            #           1918                      35253
            #    South Asian           South East Asian
            #          13782                      25342
            #           Thai                    unknown
            #           3381                     282050
saveRDS(merged_seurat, file.path(path, 'healthy_RNA_merged.rds'))
merged_seurat <- merged_seurat[common_genes, ]
saveRDS(merged_seurat, file.path(path, 'healthy_RNA_merged_common_genes.rds'))
merged_seurat # 33670 genes x 644680 cells
# 33670 features across 644680 samples within 1 assay


merged_seurat <- readRDS(file.path(path, 'healthy_RNA_merged_common_genes.rds'))
library(dplyr)
library(harmony)
merged_seurat <- merged_seurat %>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>% ScaleData()
merged_seurat <- RunPCA(merged_seurat, assay="RNA", features=rownames(merged_seurat), npcs=40)
covariates
harmonized_seurat <- RunHarmony(merged_seurat, group.by.vars = covariates, reduction = "pca", assay.use = "RNA", reduction.save = "harmony")
harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony", assay = "RNA", dims = 1:40)
harmonized_seurat

HARMONY = Embeddings(harmonized_seurat, reduction = "harmony")
UMAP = Embeddings(harmonized_seurat, reduction = "umap")
META = harmonized_seurat@meta.data
fwrite(as.data.frame(HARMONY), file.path(path, 'harmony_embedding.csv'), sep = ',', quote=F, row.names = T, col.names = T)
fwrite(as.data.frame(UMAP), file.path(path, 'harmony_umap.csv'), sep = ',', quote=F, row.names = T, col.names = T)
fwrite(as.data.frame(META), file.path(path, 'healthy_metadata.csv'), sep = ',', quote=F, row.names = T, col.names = T)
saveRDS(harmonized_seurat, file.path(path, 'healthy_RNA_merged_harmonized.rds'))
