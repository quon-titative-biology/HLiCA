# ESB
# (r_env) hrhu@esb:/share/quonlab/workspaces/hrhu/HLiCA/healthy_gamma$

library(Matrix)
library(data.table)
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(harmony)

merged_seurat <- readRDS('/share/quonlab/workspaces/hrhu/HLiCA/healthy_gamma/healthy_RNA_merged_harmonized_gamma.rds')

table(merged_seurat@meta.data$STUDY)
#  DasGupta      Grun Henderson    Mullen     Scott Toronto_1 Toronto_2 Toronto_3
#     43318     26404     92193    112703     51554     10750     56925    130852
table(merged_seurat@meta.data$donor_uuid) # dim = 110
table(merged_seurat@meta.data$library_uuid) # dim = 156
table(merged_seurat@meta.data$assay)
# 10x 3' v2 10x 3' v3 10x 5' v1
#    157663    324241     42795
table(merged_seurat@meta.data$suspension_type)
#    cell nucleus
#  367328  157371
hvgs <- VariableFeatures(merged_seurat)


merged_seurat@assays$RNA@scale.data[1:5,1:5]
#          Mullen_6854_1_AAACCCAAGACACACG-1 Mullen_6854_1_AAACCCAAGCGGATCA-1
# MTND1P23                       -0.1960194                       -0.1960194
# MTCO1P12                       -0.2158544                       -0.2158544
# HES4                           -0.2483795                       -0.2483795
# ISG15                          -0.3792342                       -0.3792342
# TNFRSF18                       -0.1131502                       -0.1131502
#          Mullen_6854_1_AAACCCACACACACGC-1 Mullen_6854_1_AAACCCACAGCGTACC-1
# MTND1P23                       -0.1960194                       -0.1960194
# MTCO1P12                       -0.2158544                       -0.2158544
# HES4                           -0.2483795                       -0.2483795
# ISG15                          -0.3792342                        2.9960194
# TNFRSF18                       -0.1131502                       -0.1131502
#          Mullen_6854_1_AAACCCACAGCTGTCG-1
# MTND1P23                       -0.1960194
# MTCO1P12                       -0.2158544
# HES4                           -0.2483795
# ISG15                          -0.3792342
# TNFRSF18                       -0.1131502
merged_seurat <- merged_seurat[hvgs, ]


disease_merged <- readRDS('/share/quonlab/workspaces/hrhu/HLiCA/disease_atlas/disease_merged_Oct2024.rds')
disease_merged <- disease_merged[hvgs, ]



disease_merged
# 4000 features across 319754 samples within 1 assay
merged_seurat
# 4000 features across 524699 samples within 1 assay



table(disease_merged@meta.data$STUDY)
        #   Amit         Aronow            Cao           Chen       DasGupta
        #   4956           3100          10459          17762          32744
        #    Fan           Fong   Heikenwalder      Henderson Henderson_2024
        #   6976           5847          80617          19295          22747
        #     Li           Lleo             Qu        Schramm        Schwabe
        #  22387          22472           1051           7205            896
        #    Shi            Sun            Yan
        #  16119           4588          40533
colnames(disease_merged@meta.data)
#  [1] "orig.ident"          "nCount_RNA"          "nFeature_RNA"
#  [4] "percent.mt"          "cell_barcode"        "donor"
#  [7] "S.Score"             "G2M.Score"           "Phase"
# [10] "old.ident"           "general_labs"        "consistent_labs"
# [13] "marker_labs"         "marker_general_labs" "Dataset ID"
# [16] "Run ID"              "Output Name"         "Assay Type"
# [19] "BioProject"          "BioSample"           "Center Name"
# [22] "Experiment"          "Platform"            "Sample Name"
# [25] "SRA Study"           "Disease Type"        "STUDY"
# [28] "PMID"                "Sequencing"          "FACS sorting"
# [31] "Pass_CellRanger"     "Pass_QC"             "Chemistry"
# [34] "Sex"

table(disease_merged@meta.data$'Sample Name') # dim = 140
disease_merged@meta.data$donor_uuid <- disease_merged@meta.data$'Sample Name'
dim(table(disease_merged@meta.data$donor_uuid)) # dim = 140


table(disease_merged@meta.data$'Output Name') # dim = 324
disease_merged@meta.data$library_uuid <- disease_merged@meta.data$'Output Name'
dim(table(disease_merged@meta.data$library_uuid)) # dim = 324


table(disease_merged@meta.data$'Chemistry') # dim = 140
#      Single Cell 3' v2      Single Cell 3' v3      Single Cell 5' PE
#                 113510                 182706                   5847
# Single Cell 5' R2-only
#                  17691
disease_merged@meta.data$assay <- disease_merged@meta.data$'Chemistry'
table(disease_merged@meta.data$assay)  # 10x 3' v2 10x 3' v3 10x 5' v1
#      Single Cell 3' v2      Single Cell 3' v3      Single Cell 5' PE
#                 113510                 182706                   5847
# Single Cell 5' R2-only
#                  17691
# Single Cell 3' v2 --> 10x 3' v2
# Single Cell 3' v3 --> 10x 3' v3
# Single Cell 5' PE --> 10x 5' PE
# Single Cell 5' R2-only --> 10x 5' R2
disease_merged@meta.data <- disease_merged@meta.data %>%
  mutate(assay = case_when(
    assay == "Single Cell 3' v2" ~ "10x 3' v2",
    assay == "Single Cell 3' v3" ~ "10x 3' v3",
    assay == "Single Cell 5' PE" ~ "10x 5' PE",
    assay == "Single Cell 5' R2-only" ~ "10x 5' R2",
    TRUE ~ assay
  ))

# Verify the changes
table(disease_merged@meta.data$assay)
# 10x 3' v2 10x 3' v3 10x 5' PE 10x 5' R2
#    113510    182706      5847     17691



table(disease_merged@meta.data$'Sequencing') # dim = 140
# scRNA-seq snRNA-seq
#    292487     27267
disease_merged@meta.data$suspension_type <- disease_merged@meta.data$'Sequencing'
table(disease_merged@meta.data$suspension_type) # cell nucleus
# scRNA-seq snRNA-seq
#    292487     27267
# scRNA-seq --> cell
# snRNA-seq --> nucleus
disease_merged@meta.data <- disease_merged@meta.data %>%
  mutate(suspension_type = case_when(
    suspension_type == "scRNA-seq" ~ "cell",
    suspension_type == "snRNA-seq" ~ "nucleus",
    TRUE ~ suspension_type
  ))

# Verify the changes
table(disease_merged@meta.data$suspension_type)
#    cell nucleus
#  292487   27267


merged_seurat@meta.data$disease_status <- 'healthy'
disease_merged@meta.data$disease_status <- 'disease'


merged_all <- merge(x = merged_seurat,
                    y = disease_merged,
                    merge.data = TRUE)
covariates <- c('STUDY', 'suspension_type','assay', 'donor_uuid', 'library_uuid', 'disease_status')
merged_all
# 4000 features across 844453 samples within 1 assay

table(merged_all@meta.data$disease_status)
# disease healthy
#  319754  524699
table(merged_all@meta.data$STUDY)
#           Amit         Aronow            Cao           Chen       DasGupta
#           4956           3100          10459          17762          76062
#            Fan           Fong           Grun   Heikenwalder      Henderson
#           6976           5847          26404          80617         111488
# Henderson_2024             Li           Lleo         Mullen             Qu
#          22747          22387          22472         112703           1051
#        Schramm        Schwabe          Scott            Shi            Sun
#           7205            896          51554          16119           4588
#      Toronto_1      Toronto_2      Toronto_3            Yan
#          10750          56925         130852          40533
table(merged_all@meta.data$assay)
# 10x 3' v2 10x 3' v3 10x 5' PE 10x 5' R2 10x 5' v1
#    271173    506947      5847     17691     42795
table(merged_all@meta.data$suspension_type)
#    cell nucleus
#  659815  184638
saveRDS(merged_all, '/share/quonlab/workspaces/hrhu/HLiCA/healthy_gamma_healthy_disease_merged.rds')


# scale data
merged_all <- ScaleData(merged_all)
merged_all <- RunPCA(merged_all, assay="RNA", features=rownames(merged_all), npcs=50)
covariates
harmonized_seurat <- RunHarmony(merged_all, group.by.vars = covariates, reduction = "pca", assay.use = "RNA", reduction.save = "harmony")
harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony", assay = "RNA", dims = 1:50)
harmonized_seurat


path <- '/share/quonlab/workspaces/hrhu/HLiCA/healthy_gamma/'
HARMONY = Embeddings(harmonized_seurat, reduction = "harmony")
UMAP = Embeddings(harmonized_seurat, reduction = "umap")

META = as.data.frame(harmonized_seurat@meta.data)
colnames(META)
# fill <NA> to "unknown" for the entire dataframe
META[is.na(META)] <- "unknown"
head(META)

META = META[,c('Phase', 'SAMPLE', 'STUDY',
        'tissue', 'suspension_type', 'assay', 'sequencing_run', 'library_alias',
        'library_uuid', 'donor_uuid', 'donor_sex', 'donor_age', 'donor_ethnicity', 'condition',
        'Gamma.Annotation', 'disease_status', 'donor', 'old.ident', 'general_labs', 'consistent_labs',
        'Output Name','Sample Name','Sequencing','Chemistry','Sex','Disease Type')]



fwrite(META, 'all_metadata.csv', sep = ',', quote=F, row.names = T, col.names = T)

fwrite(as.data.frame(HARMONY), file.path(path, 'all_harmony_embedding.csv'), sep = ',', quote=F, row.names = T, col.names = T)
fwrite(as.data.frame(UMAP), file.path(path, 'all_harmony_umap.csv'), sep = ',', quote=F, row.names = T, col.names = T)
saveRDS(harmonized_seurat, file.path(path, 'all_RNA_merged_harmonized.rds'))

