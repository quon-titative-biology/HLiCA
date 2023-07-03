# ESB
.libPaths(.libPaths()[2])
.libPaths()
library(Matrix)
library(data.table)
library(Seurat)
library(SeuratDisk)
library(harmony)


dataset = c('Gruen', 'GuilliamsScott', 'Henderson', 'Mullen', 'Toronto', 'DasGupta')
# 'DasGupta' has not extra meta data # 
metatype = c('SAMPLE', 'suspension_type', 'assay', 'tissue', 'donor_sex', 'donor_age', 'donor_ethnicity')
path = '/share/quonlab/workspaces/hrhu/HLiCA'

# for (study in dataset){
#     print(study)
#     study_list = list()
#     for (i in (1:length(Sys.glob(file.path(path,study,'*.rds'))))){
#         print(i)
#         print(Sys.glob(file.path(path,study,'*.rds'))[i])
#         seurat = readRDS(Sys.glob(file.path(path,study,'*.rds'))[i])
#         study_list[[i]] = seurat
#     }
#     merged_seurat <- merge(x = study_list[[1]],
#                 y = study_list[2:length(study_list)],
#                 merge.data = TRUE)
#     merged_seurat <- merged_seurat %>%
#         NormalizeData() %>%
#         FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
#         ScaleData()
#     merged_seurat <- RunPCA(merged_seurat, assay = "RNA", npcs = 50)
#     metatype_ = list()
#     for (label in metatype){
#         print(label)
#         print(table(merged_seurat@meta.data[[label]]))
#         if (length(table((merged_seurat@meta.data[[label]]))) > 1){
#             metatype_ = c(metatype_, label)
#         }
#     }
#     metatype_ = unlist(metatype_)
#     harmonized_seurat <- RunHarmony(merged_seurat, 
#                     group.by.vars = metatype_, 
#                     reduction = "pca", assay.use = "RNA", reduction.save = "harmony")
#     harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony", assay = "RNA", dims = 1:40)
#     # Make plots.
#     plot_list = list()
#     for (i in 1:length(metatype_)){
#         print(i)
#         p = DimPlot(harmonized_seurat, reduction = "umap", group.by = metatype_[i], label = TRUE) + NoLegend()
#         plot_list[[i]] = p
#     }
#     for (i in 1:length(metatype_)){
#         print(i)
#         filename = file.path(path, paste0(study, '_', as.character(dim(merged_seurat)[2]), '', 'cells_', metatype_[i], '.png'))
#         png(filename, width=1000, height=1000, res=300)
#         print(plot_list[[i]])
#         dev.off()
#         print("done")
#     }
#     metatype_ = c("Phase", "scmap_cluster_anno", "scmap_cell_anno", "general_labs", "consistent_labs", "marker_labs" , "marker_general_labs")
#     plot_list = list()
#     for (i in 1:length(metatype_)){
#         print(i)
#         p = DimPlot(harmonized_seurat, reduction = "umap", group.by = metatype_[i], label = TRUE) + NoLegend()
#         plot_list[[i]] = p
#     }
#     for (i in 1:length(metatype_)){
#         print(i)
#         filename = file.path(path, paste0(study, '_', as.character(dim(merged_seurat)[2]), '', 'cells_', metatype_[i], '.png'))
#         png(filename, width=1000, height=1000, res=300)
#         print(plot_list[[i]])
#         dev.off()
#         print("done")
#     }
# }





# # higher resolution???


# for (study in dataset){
#     print(study)
#     study_list = list()
#     for (i in (1:length(Sys.glob(file.path(path,study,'*.rds'))))){
#         print(i)
#         print(Sys.glob(file.path(path,study,'*.rds'))[i])
#         seurat = readRDS(Sys.glob(file.path(path,study,'*.rds'))[i])
#         study_list[[i]] = seurat
#     }
#     merged_seurat <- merge(x = study_list[[1]],
#                 y = study_list[2:length(study_list)],
#                 merge.data = TRUE)
#     merged_seurat <- merged_seurat %>%
#         NormalizeData() %>%
#         FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
#         ScaleData()
#     merged_seurat <- RunPCA(merged_seurat, assay = "RNA", npcs = 50)
#     metatype_ = list()
#     for (label in metatype){
#         print(label)
#         print(table(merged_seurat@meta.data[[label]]))
#         if (length(table((merged_seurat@meta.data[[label]]))) > 1){
#             metatype_ = c(metatype_, label)
#         }
#     }
#     metatype_ = unlist(metatype_)
#     harmonized_seurat <- RunHarmony(merged_seurat, 
#                     group.by.vars = metatype_, 
#                     reduction = "pca", assay.use = "RNA", reduction.save = "harmony")
#     harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony", assay = "RNA", dims = 1:40)
#     # Make plots.
#     plot_list = list()
#     for (i in 1:length(metatype_)){
#         print(i)
#         p = DimPlot(harmonized_seurat, reduction = "umap", group.by = metatype_[i], label = TRUE) + NoLegend()
#         plot_list[[i]] = p
#     }
#     for (i in 1:length(metatype_)){
#         print(i)
#         filename = file.path(path, paste0(study, '_', as.character(dim(merged_seurat)[2]), '', 'cells_', metatype_[i], '_.png'))
#         png(filename, width=4000, height=4000, res=300)
#         print(plot_list[[i]])
#         dev.off()
#         print("done")
#     }
#     metatype_ = c("Phase", "scmap_cluster_anno", "scmap_cell_anno", "general_labs", "consistent_labs", "marker_labs" , "marker_general_labs")
#     plot_list = list()
#     for (i in 1:length(metatype_)){
#         print(i)
#         p = DimPlot(harmonized_seurat, reduction = "umap", group.by = metatype_[i], label = TRUE) + NoLegend()
#         plot_list[[i]] = p
#     }
#     for (i in 1:length(metatype_)){
#         print(i)
#         filename = file.path(path, paste0(study, '_', as.character(dim(merged_seurat)[2]), '', 'cells_', metatype_[i], '_.png'))
#         png(filename, width=4000, height=4000, res=300)
#         print(plot_list[[i]])
#         dev.off()
#         print("done")
#     }
# }






# metatype = c('study', 'SAMPLE', 'suspension_type', 'assay', 'tissue', 'donor_sex', 'donor_age', 'donor_ethnicity')
metatype = c('study', 'SAMPLE')#, 'suspension_type', 'assay')

study_list = list()
for (i in (1:length(Sys.glob(file.path(path,'*','*.rds'))))){
    print(i)
    print(Sys.glob(file.path(path,'*','*.rds'))[i])
    seurat = readRDS(Sys.glob(file.path(path,'*','*.rds'))[i])
    seurat@meta.data$study = basename(dirname(Sys.glob(file.path(path,'*','*.rds'))[i]))
    study_list[[i]] = seurat
}
merged_seurat <- merge(x = study_list[[1]],
            y = study_list[2:length(study_list)],
            merge.data = TRUE)
merged_seurat <- merged_seurat %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData()
merged_seurat <- RunPCA(merged_seurat, assay = "RNA", npcs = 50)
metatype_ = list()
for (label in metatype){
    print(label)
    print(table(merged_seurat@meta.data[[label]]))
    if (length(table((merged_seurat@meta.data[[label]]))) > 1){
        metatype_ = c(metatype_, label)
    }
}
metatype_ = unlist(metatype_)
metatype_ = c('study', 'SAMPLE')
harmonized_seurat <- RunHarmony(merged_seurat, 
                group.by.vars = metatype_, 
                reduction = "pca", assay.use = "RNA", reduction.save = "harmony")
harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony", assay = "RNA", dims = 1:40)
harmonized_seurat
# # Make plots.
# plot_list = list()
# for (i in 1:length(metatype_)){
#     print(i)
#     p = DimPlot(harmonized_seurat, reduction = "umap", group.by = metatype_[i], label = TRUE) + NoLegend()
#     plot_list[[i]] = p
# }
# for (i in 1:length(metatype_)){
#     print(i)
#     filename = file.path(path, paste0('all_', as.character(dim(merged_seurat)[2]), '', 'cells_', metatype_[i], '_.png'))
#     png(filename, width=4000, height=4000, res=300)
#     print(plot_list[[i]])
#     dev.off()
#     print("done")
# }
# metatype_ = c("Phase", "scmap_cluster_anno", "scmap_cell_anno", "general_labs", "consistent_labs", "marker_labs" , "marker_general_labs")
# plot_list = list()
# for (i in 1:length(metatype_)){
#     print(i)
#     p = DimPlot(harmonized_seurat, reduction = "umap", group.by = metatype_[i], label = TRUE) + NoLegend()
#     plot_list[[i]] = p
# }
# for (i in 1:length(metatype_)){
#     print(i)
#     filename = file.path(path, paste0('all_', as.character(dim(merged_seurat)[2]), '', 'cells_', metatype_[i], '_.png'))
#     png(filename, width=4000, height=4000, res=300)
#     print(plot_list[[i]])
#     dev.off()
#     print("done")
# }

UMAP = Embeddings(harmonized_seurat, reduction = "umap")
META = harmonized_seurat@meta.data

fwrite(as.data.frame(UMAP), file.path(path, 'all_umap.csv'), sep = ',', quote=F, row.names = T, col.names = T)
fwrite(as.data.frame(META), file.path(path, 'all_meta.csv'), sep = ',', quote=F, row.names = T, col.names = T)

merged_seurat     # 62696 features across 595369 samples within 1 assay
harmonized_seurat # 62696 features across 595369 samples within 1 assay
table(merged_seurat@meta.data$study)
    #   DasGupta          Gruen GuilliamsScott      Henderson         Mullen
    #     183100          29242          59201         106336         137360
    #    Toronto
    #      80130
table(harmonized_seurat@meta.data$study)




saveRDS(harmonized_seurat, 'all_harmonized_seurat.rds')

saveRDS(merged_seurat, 'all_merged_seurat.rds')
SaveH5Seurat(harmonized_seurat, filename = "all_harmonized_seurat.h5Seurat")
Convert("all_harmonized_seurat.h5Seurat", dest = "h5ad")
SaveH5Seurat(merged_seurat, filename = "all_merged_seurat.h5Seurat")
Convert("all_merged_seurat.h5Seurat", dest = "h5ad")



marker_genes = read.csv('/share/quonlab/workspaces/hrhu/HLiCA/marker_genes.csv', header = T, sep = ',')
head(marker_genes)
#    Gene Specific.Type General.Type
# 1 CD79A        B cell       B cell
# 2 CD79B        B cell       B cell
# 3  IGHM  Naive B cell       B cell
# 4 IGHA1  Plama B cell       B cell
# 5 IGHA2  Plama B cell       B cell
# 6 IGHG1  Plama B cell       B cell
marker_genes[] <- lapply(marker_genes, function(x) gsub("/", "_", x))
marker_genes[] <- lapply(marker_genes, function(x) gsub(" ", "_", x))

table(marker_genes$Specific.Type)
table(marker_genes$General.Type)
#               B cell            Cellcycle        Cholangiocyte
#                   14                    7                   30
#          Endothelial           Hepatocyte           Lymphocyte
#                   34                   27                   31
#   Macrophage/Myeloid                  RBC Stellate/Mesenchymal
#                   73                    4                   22

#               B_cell            Cellcycle        Cholangiocyte
#                   14                    7                   30
#          Endothelial           Hepatocyte           Lymphocyte
#                   34                   27                   31
#   Macrophage_Myeloid                  RBC Stellate_Mesenchymal
#                   73                    4                   22

spec_type = unique(marker_genes$Specific.Type)
general_type = unique(marker_genes$General.Type)
# feature plot marker genes of each type for all
feature_plots = list()
for (i in 1:length(spec_type)){
    print(spec_type[i])
    marker_genes_ = marker_genes[marker_genes$Specific.Type == spec_type[i],]$Gene
    p = FeaturePlot(harmonized_seurat, features = marker_genes_, ncol=1, raster = FALSE)
    feature_plots[[i]] = p
    print("generating figure")
    png(paste0('all_', spec_type[i], '_Spec.png'), width=1000, height=1000*length(marker_genes_),res=300)
    print(feature_plots[[i]])
    dev.off()
    print("done")
}

# feature_plots_ = list()
# for (i in 1:length(general_type)){
#     print(general_type[i])
#     marker_genes_ = marker_genes[marker_genes$General.Type == general_type[i],]$Gene
#     p = FeaturePlot(harmonized_seurat, features = marker_genes_,  ncol=1, raster = FALSE)
#     feature_plots_[[i]] = p
#     print("generating figure")
#     png(paste0('all_', general_type[i], '_General.png'), width=1000, height=1000*length(marker_genes_), res=300)
#     print(feature_plots_[[i]])
#     dev.off()
#     print("done")
# }



# feature plot marker genes of each type for each study
for (study in dataset){
    obj = harmonized_seurat[, which(harmonized_seurat@meta.data$study == study)]
    print(study)
    for (i in 1:length(spec_type)){
        print(spec_type[i])
        marker_genes_ = marker_genes[marker_genes$Specific.Type == spec_type[i],]$Gene
        p = FeaturePlot(obj, features = marker_genes_, ncol=1, raster = FALSE)
        print("generating figure...")
        png(paste0(study, '_', spec_type[i], '_Spec.png'), width=1000, height=1000*length(marker_genes_),res=300)
        print(p)
        dev.off()
        print("done")
    }
}




# load libraries and functions
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# get cell-type-specific gene sets from our in-built database (DB)
gs_list = gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Liver") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Liver" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
# prepare gene sets
gs_list_ = gene_sets_prepare(db_, tissue)


harmonized_seurat <- FindNeighbors(harmonized_seurat, dims = 1:40)
harmonized_seurat <- FindClusters(harmonized_seurat, resolution = 0.5)

table(harmonized_seurat@meta.data$seurat_clusters)


png('cluster_umap.png', width=4000, height=4000, res=300)
DimPlot(harmonized_seurat, reduction = "umap", group.by = "seurat_clusters", label = TRUE, raster = FALSE)
dev.off()





###############################################
###############################################
    # scType #
###############################################
###############################################
head(harmonized_seurat@meta.data)
harmonized_seurat
es.max = sctype_score(scRNAseqData=harmonized_seurat[["RNA"]]@scale.data, scaled=TRUE, gs=gs_list$gs_positive, gs2=gs_list$gs_negative)
es.max_ = sctype_score(scRNAseqData=harmonized_seurat[["RNA"]]@scale.data, scaled=TRUE, gs=gs_list_$gs_positive, gs2=gs_list_$gs_negative)
# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either obj[["RNA"]]@scale.data (default), obj[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or obj[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(harmonized_seurat@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(harmonized_seurat@meta.data[harmonized_seurat@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(harmonized_seurat@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])


harmonized_seurat@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  harmonized_seurat@meta.data$customclassif[harmonized_seurat@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
png('cluster_umap_customclassif_scType.png', width=4000, height=4000, res=300)
DimPlot(harmonized_seurat, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif', raster = FALSE)    
dev.off()

fwrite(as.data.frame(harmonized_seurat@meta.data$customclassif), file.path(path, 'scType.csv'), sep = ',', quote=F, row.names = T, col.names = T)


png('cluster_umap_ct.png', width=6000, height=1000, res=300)
DimPlot(harmonized_seurat, reduction = "umap", group.by = 'study', split.by='study', ncol=6, raster = FALSE)
dev.off()

metatype_ = c("Phase", "scmap_cluster_anno", "scmap_cell_anno", "general_labs", "consistent_labs", "marker_labs" , "marker_general_labs", "customclassif")

for (mt_ in metatype_){
    print(mt_)
    png(paste0(mt_,'_cluster_umap_ct.png'), width=12000, height=4000, res=400)
    print(DimPlot(harmonized_seurat, reduction = "umap", group.by = mt_, split.by='study', label=T, repel=T, ncol=6, raster = FALSE))
    dev.off()
    print('done')
}
