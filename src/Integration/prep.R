
# ESB
.libPaths(.libPaths()[2])
.libPaths()
library(Matrix)
library(data.table)
library(Seurat)
library(SeuratDisk)


dataset = c('Gruen', 'GuilliamsScott', 'Henderson', 'Mullen', 'Toronto') 
# 'DasGupta' has not extra meta data # 
metatype = c('SAMPLE', 'suspension_type', 'assay', 'tissue', 'donor_sex', 'donor_age', 'donor_ethnicity')
path = "/share/quonlab/workspaces/czi_liver_atlas/data/LiverNetworkData/alignment/ref_GRCh38p13_gencode_v42/share_2023_05_13/share_qc"
output_path = '/share/quonlab/workspaces/hrhu/HLiCA'
for (study in dataset){
    print(study)
    metadata = read.csv(file.path(path,study,'GoogleSheetMetadata_sample_qced.csv'), header=TRUE, sep=",")
    metadata = metadata[metadata$AlignmentCellRanger=='True',]
    metadata = metadata[metadata$QC_Pipeline=='True',]
    SAMPLE_ID = metadata$SAMPLE
    rownames(metadata) = metadata$SAMPLE
    for (sample in SAMPLE_ID){
        print(sample)
        filename = Sys.glob(file.path(path,study,sample,'raw_feature_bc_matrix/qc_output/Cleaned_output_EmptyOnly.rds'))
        processed_RNA = readRDS(filename)
        print(dim(processed_RNA))
        meta_name = metadata[rownames(metadata)==sample,]
        for (meta_col in metatype){
            print(meta_col)
            processed_RNA@meta.data[meta_col] = meta_name[meta_col]
        }
        prefix = paste0(study, '_', sample)
        processed_RNA = RenameCells(processed_RNA, add.cell.id=prefix)
        processed_RNA@meta.data$orig.ident = prefix
        processed_RNA@meta.data$donor = NULL
        processed_RNA@meta.data$cell_ID = NULL
        processed_RNA@meta.data$scmap_cluster_anno = processed_RNA@meta.data$scmap_cluster_anno$lm1
        updated_metadata = processed_RNA@meta.data
        counts = processed_RNA@assays$RNA@counts
        RNA = CreateSeuratObject(counts=counts, meta.data=updated_metadata)
        print('saving files...')
        saveRDS(RNA, paste0(file.path(output_path, study), '/', sample, "_RNA.rds"))
        SaveH5Seurat(RNA, filename = paste0(file.path(output_path, study), '/', sample, "_RNA.h5Seurat"), overwrite = TRUE)
        Convert(paste0(file.path(output_path, study), '/', sample, "_RNA.h5Seurat"), dest = "h5ad", overwrite = TRUE)
        fwrite(updated_metadata, file = paste0(file.path(output_path, study), '/', sample, "_metadata.csv"), quote=F, row.names=T, col.names=T)
        print('files saved')
    }
}








study = 'DasGupta'
print(study)
metadata = read.csv(file.path(path,study,'GoogleSheetMetadata_sample_qced.csv'), header=TRUE, sep=",")
extra_metatype = c('suspension_type', 'assay', 'tissue', 'donor_sex', 'donor_age', 'donor_ethnicity')
for (extra_meta_col in extra_metatype){
     metadata[extra_meta_col] = "unknown_"
}
metadata = metadata[metadata$AlignmentCellRanger=='True',]
metadata = metadata[metadata$QC_Pipeline=='True',]
SAMPLE_ID = metadata$SAMPLE
rownames(metadata) = metadata$SAMPLE
for (sample in SAMPLE_ID){
    print(sample)
    filename = Sys.glob(file.path(path,study,sample,'raw_feature_bc_matrix/qc_output/Cleaned_output_EmptyOnly.rds'))
    processed_RNA = readRDS(filename)
    print(dim(processed_RNA))
    meta_name = metadata[rownames(metadata)==sample,]
    for (meta_col in metatype){
        print(meta_col)
        processed_RNA@meta.data[meta_col] = meta_name[meta_col]
    }
    prefix = paste0(study, '_', sample)
    processed_RNA = RenameCells(processed_RNA, add.cell.id=prefix)
    processed_RNA@meta.data$orig.ident = prefix
    processed_RNA@meta.data$donor = NULL
    processed_RNA@meta.data$cell_ID = NULL
    processed_RNA@meta.data$scmap_cluster_anno = processed_RNA@meta.data$scmap_cluster_anno$lm1
    updated_metadata = processed_RNA@meta.data
    counts = processed_RNA@assays$RNA@counts
    RNA = CreateSeuratObject(counts=counts, meta.data=updated_metadata)
    print('saving files...')
    saveRDS(RNA, paste0(file.path(output_path, study), '/', sample, "_RNA.rds"))
    SaveH5Seurat(RNA, filename = paste0(file.path(output_path, study), '/', sample, "_RNA.h5Seurat"), overwrite = TRUE)
    Convert(paste0(file.path(output_path, study), '/', sample, "_RNA.h5Seurat"), dest = "h5ad", overwrite = TRUE)
    fwrite(updated_metadata, file = paste0(file.path(output_path, study), '/', sample, "_metadata.csv"), quote=F, row.names=T, col.names=T)
    print('files saved')
}
